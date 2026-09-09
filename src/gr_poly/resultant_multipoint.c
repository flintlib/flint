/*
    Copyright (C) 2026 Mael Hostettler
    Copyright (C) 2026 Antoine Bak
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"
#include "gr_poly.h"
#include "thread_support.h"
#if FLINT_HAVE_FFT_SMALL
# include "fft_small.h"
#endif

/* Target size, in words, of the tables holding the evaluated coefficients.
   The evaluation points are processed in blocks chosen so that these tables
   stay around this size; the memory usage is then proportional to the input
   and output sizes instead of to their product. */
#define RESULTANT_MULTIPOINT_BLOCK_WORDS (WORD(1) << 18)

/* When the modulus supports a radix-2 transform long enough for the number of
   points, the coefficients are evaluated at roots of unity with DFTs instead
   of the Bluestein product of the geometric method; measured on this machine,
   that phase becomes six to twelve times faster. */

/* The evaluation of the coefficients in y, and the univariate resultants at
   the points, are both loops of independent iterations, and both are split
   over several workers when there is enough to do. The work estimates are in
   units of word operations; under this many, the thread pool costs more than
   it saves. */
#define RESULTANT_MULTIPOINT_PARALLEL_CUTOFF (WORD(1) << 16)

static slong
_num_workers(double work)
{
    slong nworkers = flint_get_num_threads();

    if (nworkers <= 1 || work < (double) RESULTANT_MULTIPOINT_PARALLEL_CUTOFF)
        return 1;

    return nworkers;
}

/* Returns 1 if none of q, q^2, ..., q^(len-1) is one, i.e. if the len
   points 1, q, q^2, ..., q^(len-1) are pairwise distinct. */
static int
_nmod_geometric_points_distinct(ulong q, slong len, nmod_t mod)
{
    ulong qi = 1;
    slong i;

    for (i = 1; i < len; i++)
    {
        qi = nmod_mul(qi, q, mod);

        if (qi == 1)
            return 0;
    }

    return 1;
}

/* Returns r such that q = r^2 is of order at least len, or 0 if there is no
   such r. The group of squares is cyclic of order (n - 1) / 2, so one exists
   exactly when 2 len <= n - 1, and trying 2, 3, 4, ... in turn then finds one
   after a couple of candidates: the squares of order exactly (n - 1) / 2 by
   themselves have density phi(m) / m with m = (n - 1) / 2, which stays above
   1/8 for every m that fits in a word. Going through the candidates in order
   rather than at random keeps the algorithm deterministic. */
static ulong
_nmod_find_geometric_ratio(slong len, nmod_t mod)
{
    ulong r;

    if (mod.n < 2 * (ulong) len + 1)
        return 0;

    for (r = 2; r < mod.n; r++)
        if (_nmod_geometric_points_distinct(nmod_mul(r, r, mod), len, mod))
            return r;

    return 0;
}

/* Sets res to res_y(A, B) evaluated at one point, given in f and g the values
   there of the lenA and lenB coefficients in y. A specialisation may drop in
   degree; res_y(A, B) is nevertheless the determinant of the Sylvester matrix
   built from the degrees lenA - 1 and lenB - 1, which is what the corrections
   below reduce to a resultant of the specialisations. */
static void
_resultant_at_point(ulong * res, nn_ptr f, slong lenA,
                    nn_ptr g, slong lenB, nmod_t mod)
{
    slong len1 = lenA, len2 = lenB;
    ulong t;

    MPN_NORM(f, len1);
    MPN_NORM(g, len2);

    /* one of them vanishes identically: the Sylvester matrix has a zero row */
    if (len1 == 0 || len2 == 0)
    {
        *res = 0;
        return;
    }

    if (len1 == lenA)
    {
        /* deg_y(A) is preserved:
           res(A, B) = lc(A)^(deg_y(B) - deg(g)) res(f, g) */
        t = _nmod_poly_resultant(f, len1, g, len2, mod);

        if (len2 < lenB)
            t = nmod_mul(t, nmod_pow_ui(f[len1 - 1], lenB - len2, mod), mod);
    }
    else if (len2 == lenB)
    {
        /* deg_y(B) is preserved instead; the same identity applied to
           res(B, A) = (-1)^(deg_y(A) deg_y(B)) res(A, B) gives
           res(A, B) = (-1)^(deg_y(B) (deg_y(A) + deg(f)))
                           lc(B)^(deg_y(A) - deg(f)) res(f, g) */
        if (len1 >= len2)
        {
            t = _nmod_poly_resultant(f, len1, g, len2, mod);
        }
        else
        {
            t = _nmod_poly_resultant(g, len2, f, len1, mod);

            if (((len1 - 1) & 1) && ((len2 - 1) & 1))
                t = nmod_neg(t, mod);
        }

        t = nmod_mul(t, nmod_pow_ui(g[len2 - 1], lenA - len1, mod), mod);

        if (((lenB - 1) & 1) && ((lenA + len1) & 1))
            t = nmod_neg(t, mod);
    }
    else
    {
        /* both drop: the first column of that Sylvester matrix holds only the
           two leading coefficients, so it is zero and so is the determinant */
        t = 0;
    }

    *res = t;
}


/* One block of univariate resultants. The values of the lenA + lenB
   coefficients in y at the point held in position s are val[i * stride + s].
   Where each result goes depends on how the caller interpolates: at the
   point index base + step * s, or, when the caller is going to run an
   inverse transform, at position s of its buffer. */
typedef struct
{
    nn_srcptr valA;
    nn_srcptr valB;
    slong lenA;
    slong lenB;
    slong stride;
    slong blen;
    nn_ptr res;
    double * dres;
    const ulong * sidx;     /* NULL, or the point index of position s */
    slong base;
    slong step;
    slong npoints;
    nn_ptr fbuf;
    slong nworkers;
    nmod_t mod;
}
_multipoint_res_args_t;

static void
_multipoint_res_worker(slong widx, void * argsv)
{
    _multipoint_res_args_t * A = argsv;
    nn_ptr f = A->fbuf + widx * (A->lenA + A->lenB);
    nn_ptr g = f + A->lenA;
    slong i, s;

    for (s = widx; s < A->blen; s += A->nworkers)
    {
        slong j = A->base + A->step * (A->sidx == NULL ? s : (slong) A->sidx[s]);

        if (A->sidx != NULL && j >= A->npoints)
            continue;

        for (i = 0; i < A->lenA; i++)
            f[i] = A->valA[i * A->stride + s];

        for (i = 0; i < A->lenB; i++)
            g[i] = A->valB[i * A->stride + s];

        if (A->dres != NULL)
        {
            ulong res;

            _resultant_at_point(&res, f, A->lenA, g, A->lenB, A->mod);
            A->dres[s] = (double) res;
        }
        else
        {
            _resultant_at_point(A->res + j, f, A->lenA, g, A->lenB, A->mod);
        }
    }
}

/* One block of the geometric evaluation: the i-th coefficient in y, scaled
   by spow, evaluated at blen points into val + i * stride. One worker takes
   every nworkers-th coefficient, with its own scratch; the progression G is
   only read, so all of them can share it. */
typedef struct
{
    nn_ptr valA;
    nn_ptr valB;
    const gr_poly_struct * Ax;
    const gr_poly_struct * Bx;
    slong lenA;
    slong lenB;
    slong stride;
    slong blen;
    slong maxlen;
    nn_srcptr spow;
    nn_ptr tmpbuf;
    slong nworkers;
    nmod_geometric_progression_struct * G;
    nmod_t mod;
}
_geom_eval_args_t;

static void
_geom_eval_worker(slong widx, void * argsv)
{
    _geom_eval_args_t * A = argsv;
    nn_ptr tmp = A->tmpbuf + widx * A->maxlen;
    slong i, k;

    for (i = widx; i < A->lenA + A->lenB; i += A->nworkers)
    {
        const gr_poly_struct * P;
        nn_ptr out;

        if (i < A->lenA)
        {
            P = A->Ax + i;
            out = A->valA + i * A->stride;
        }
        else
        {
            P = A->Bx + (i - A->lenA);
            out = A->valB + (i - A->lenA) * A->stride;
        }

        for (k = 0; k < P->length; k++)
            tmp[k] = nmod_mul(((nn_srcptr) P->coeffs)[k], A->spow[k], A->mod);

        _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(out, tmp,
            P->length, A->G, A->blen, A->mod);
    }
}

#if FLINT_HAVE_FFT_SMALL

/* sd_fft leaves its outputs in (-2p, 2p), as exact integers in a double. */
static ulong
_sd_fft_get_nmod(double a, ulong p)
{
    slong s = (slong) a;

    if (s < 0) s += (slong) p;
    if (s < 0) s += (slong) p;
    if ((ulong) s >= p) s -= (slong) p;

    return (ulong) s;
}

/* Evaluates (poly, plen), with its k-th coefficient scaled by spow[k], at the
   M-th roots of unity, leaving the values in the raw order sd_fft_trunc
   produces them in: the value at the s-th root lands at the position
   sd_fft_ctx_trunc_index(m, n_revbin(s, m)). plen is at most M, so the
   coefficients need no reduction modulo x^M - 1 first. */
static void
_dft_evaluate(nn_ptr vs, nn_srcptr poly, slong plen, nn_srcptr spow,
              slong M, ulong m, double * dbuf, sd_fft_ctx_t Q, nmod_t mod)
{
    slong k;

    if (plen == 0)
    {
        _nmod_vec_zero(vs, M);
        return;
    }

    for (k = 0; k < plen; k++)
        dbuf[k] = (double) nmod_mul(poly[k], spow[k], mod);

    sd_fft_trunc(Q, dbuf, m, plen, M);

    for (k = 0; k < M; k++)
        vs[k] = _sd_fft_get_nmod(dbuf[k], mod.n);
}

/* One block of the transform evaluation, one worker per every nworkers-th
   coefficient in y. Each has its own transform buffer; the context Q is only
   grown under its own lock, so all of them can share it. */
typedef struct
{
    nn_ptr valA;
    nn_ptr valB;
    const gr_poly_struct * Ax;
    const gr_poly_struct * Bx;
    slong lenA;
    slong lenB;
    slong M;
    ulong m;
    nn_srcptr spow;
    double * dbuf;
    slong nworkers;
    sd_fft_ctx_struct * Q;
    nmod_t mod;
}
_dft_eval_args_t;

static void
_dft_eval_worker(slong widx, void * argsv)
{
    _dft_eval_args_t * A = argsv;
    double * dbuf = A->dbuf + widx * A->M;
    slong i;

    for (i = widx; i < A->lenA + A->lenB; i += A->nworkers)
    {
        const gr_poly_struct * P;
        nn_ptr out;

        if (i < A->lenA)
        {
            P = A->Ax + i;
            out = A->valA + i * A->M;
        }
        else
        {
            P = A->Bx + (i - A->lenA);
            out = A->valB + (i - A->lenA) * A->M;
        }

        _dft_evaluate(out, (nn_srcptr) P->coeffs, P->length, A->spow,
            A->M, A->m, dbuf, A->Q, A->mod);
    }
}

/* res_y(A, B) by evaluation at roots of unity. Returns GR_UNABLE when the
   modulus does not support a transform long enough for the number of points.

   The values are interpolated back either with one inverse transform, or
   with the geometric interpolation, whichever is cheaper; see below. */
static int
_gr_poly_resultant_multipoint_dft(gr_poly_struct * resx,
                                  const gr_poly_struct * Ax, slong lenA,
                                  const gr_poly_struct * Bx, slong lenB,
                                  slong npoints, slong maxlen, nmod_t mod,
                                  gr_ctx_t cctx)
{
    sd_fft_ctx_t Q;
    nmod_geometric_progression_t G;
    double * dbuf, * d = NULL;
    nn_ptr valA, valB, valres = NULL, resc = NULL, spow, wpow = NULL, f;
    ulong * sidx = NULL;
    slong k, t, u, N, M, lgnb, nblocks, nblocks_needed, ptrunc, nworkers;
    _dft_eval_args_t eargs;
    _multipoint_res_args_t rargs;
    ulong depth, m, v, w, r = 0;
    int use_ifft;

    depth = FLINT_MAX(n_clog2((ulong) npoints), 4);
    N = n_pow2(depth);

    if (!fft_small_mulmod_satisfies_bounds(mod.n))
        return GR_UNABLE;

    /* A primitive N-th root of unity w must exist. The geometric
       interpolation additionally needs a square root r of it, which is of
       order 2N and so costs one more power of two; the inverse transform
       needs no such thing. */
    v = n_trailing_zeros(mod.n - 1);

    if (v < depth || depth > SD_FFT_CTX_W2TAB_SIZE)
        return GR_UNABLE;

    /* The inverse transform reads back a whole number of blocks of BLK_SZ
       values, so it needs ptrunc >= npoints of them where the geometric
       interpolation needs exactly npoints. Those few extra resultants are
       all that it costs, and they pay for themselves only when they are a
       small part of the whole; below BLK_SZ points it would have to invert
       the full length N, which is never worth it. */
    if (depth > LG_BLK_SZ)
        ptrunc = FLINT_MIN(N, (slong) (n_cdiv((ulong) npoints, BLK_SZ) * BLK_SZ));
    else
        ptrunc = N;

    use_ifft = ((ptrunc - npoints) * 64 < npoints);

    /* with no square root of w available the transform is the only option */
    if (v < depth + 1)
        use_ifft = 1;

    sd_fft_ctx_init_prime(Q, mod.n);

    /* The context provides a primitive 2^v-th root, of which w and r are
       powers, so binary powering gives both. Reading them off the table of
       roots instead would mean building it up to depth, which is
       proportional to the number of points, while the transforms below only
       ever run at depth m and fit the table to that themselves. */
    w = nmod_pow_ui(Q->primitive_2power_root, n_pow2(v - depth), mod);

    if (!use_ifft)
        r = nmod_pow_ui(Q->primitive_2power_root, n_pow2(v - depth - 1), mod);

    /* The N points are split into nblocks = N / M classes modulo nblocks, the
       block t holding the points w^(t + nblocks s) for 0 <= s < M. Those are
       the M-th roots of unity scaled by w^t, so the block is evaluated by the
       substitution x -> w^t x followed by a transform of length M only. The
       block is a power of two at least maxlen, so that a coefficient always
       fits in one, and otherwise as large as the memory target allows; the
       tables of values are then proportional to the input and output sizes
       instead of to their product, as in the geometric evaluation. */
    m = FLINT_MAX(n_clog2((ulong) maxlen), 4);

    while (m < depth && (double) (lenA + lenB) * n_pow2(m + 1)
                            <= (double) RESULTANT_MULTIPOINT_BLOCK_WORDS)
        m++;

    M = n_pow2(m);
    lgnb = depth - m;
    nblocks = n_pow2(lgnb);

    nworkers = flint_get_num_threads();

    dbuf = flint_aligned_alloc(32, FLINT_MAX(32, nworkers * M * sizeof(double)));
    valA = _nmod_vec_init(lenA * M);
    valB = _nmod_vec_init(lenB * M);
    spow = _nmod_vec_init(maxlen);
    f = _nmod_vec_init(nworkers * (lenA + lenB));

    if (use_ifft)
    {
        /* Writing j = t + nblocks s, so that t is the low lgnb bits of j and
           s the high m bits, n_revbin(j, depth) = n_revbin(s, m)
           + M n_revbin(t, lgnb); and since M >= 16, trunc_index only shuffles
           the low four bits, so
             trunc_index(depth, revbin(j)) = trunc_index(m, revbin(s))
                                               + M revbin(t).
           The first term is the position the length-M transform already put
           that value in, so visiting the blocks in the order t = revbin(u)
           and copying each block's values verbatim to d[u M ..] builds
           exactly the buffer the inverse transform expects, with no
           permutation anywhere. */
        d = flint_aligned_alloc(32, FLINT_MAX(32, N * sizeof(double)));
        /* block u fills d[u M .. u M + M), so only the blocks below
           ptrunc are needed */
        nblocks_needed = (ptrunc + M - 1) / M;
    }
    else
    {
        /* sidx[s] is the power of the M-th root whose value the transform
           leaves in position s, that is the inverse of the map above */
        sidx = flint_malloc(M * sizeof(ulong));
        for (k = 0; k < M; k++)
            sidx[k] = n_revbin(sd_fft_ctx_trunc_index(m, k), m);

        valres = _nmod_vec_init(npoints);
        resc = _nmod_vec_init(npoints);
        wpow = _nmod_vec_init(maxlen);

        /* spow[k] = (w^t)^k for the current block; multiplying by
           wpow[k] = w^k advances it by one block */
        for (k = 0; k < maxlen; k++)
            spow[k] = 1;

        wpow[0] = 1;
        for (k = 1; k < maxlen; k++)
            wpow[k] = nmod_mul(wpow[k - 1], w, mod);

        /* the points of a block are spread over the whole range, so every
           block holds some below npoints and all of them are needed */
        nblocks_needed = nblocks;
    }

    eargs.Ax = Ax; eargs.Bx = Bx;
    eargs.lenA = lenA; eargs.lenB = lenB;
    eargs.valA = valA; eargs.valB = valB;
    eargs.M = M; eargs.m = m;
    eargs.spow = spow; eargs.dbuf = dbuf;
    eargs.Q = Q; eargs.mod = mod;

    rargs.valA = valA; rargs.valB = valB;
    rargs.lenA = lenA; rargs.lenB = lenB;
    rargs.stride = M;
    rargs.res = valres; rargs.dres = NULL;
    rargs.sidx = sidx; rargs.step = nblocks;
    rargs.npoints = npoints; rargs.base = 0;
    rargs.fbuf = f; rargs.mod = mod;

    for (u = 0; u < nblocks_needed; u++)
    {
        if (use_ifft)
        {
            t = n_revbin(u, lgnb);

            /* the blocks are not visited in increasing order of t, so the
               scaling is recomputed rather than advanced */
            spow[0] = 1;

            if (maxlen > 1)
            {
                ulong wt = nmod_pow_ui(w, t, mod);

                spow[1] = wt;
                for (k = 2; k < maxlen; k++)
                    spow[k] = nmod_mul(spow[k - 1], wt, mod);
            }
        }
        else
        {
            t = u;
        }

        eargs.nworkers = FLINT_MIN(nworkers,
            _num_workers((double) (lenA + lenB) * M * m));
        flint_parallel_do(_dft_eval_worker, &eargs, eargs.nworkers, 0, 0);

        if (use_ifft)
        {
            rargs.blen = FLINT_MIN(M, ptrunc - u * M);
            rargs.dres = d + u * M;
        }
        else
        {
            rargs.blen = M;
            rargs.base = t;
            rargs.dres = NULL;
        }

        rargs.nworkers = FLINT_MIN(nworkers,
            _num_workers((double) rargs.blen * lenA * lenB));
        flint_parallel_do(_multipoint_res_worker, &rargs, rargs.nworkers, 0, 0);

        if (!use_ifft)
            for (k = 1; k < maxlen; k++)
                spow[k] = nmod_mul(spow[k], wpow[k], mod);
    }

    gr_poly_fit_length(resx, npoints, cctx);
    _gr_poly_set_length(resx, npoints, cctx);

    if (use_ifft)
    {
        /* the inverse transform leaves the result scaled by 2^depth */
        ulong Ninv = nmod_inv(nmod_pow_ui(UWORD(2), depth, mod), mod);

        sd_ifft_trunc(Q, d, depth, ptrunc);

        for (k = 0; k < npoints; k++)
            ((nn_ptr) resx->coeffs)[k] =
                nmod_mul(_sd_fft_get_nmod(d[k], mod.n), Ninv, mod);
    }
    else
    {
        _nmod_geometric_progression_init_function(G, r, npoints, mod, UWORD(2));
        _nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(resc, valres, G, npoints, mod);
        nmod_geometric_progression_clear(G);

        _nmod_vec_set((nn_ptr) resx->coeffs, resc, npoints);
    }

    _gr_poly_normalise(resx, cctx);

    sd_fft_ctx_clear(Q);
    flint_aligned_free(dbuf);
    _nmod_vec_clear(valA);
    _nmod_vec_clear(valB);
    _nmod_vec_clear(spow);
    _nmod_vec_clear(f);

    if (use_ifft)
    {
        flint_aligned_free(d);
    }
    else
    {
        flint_free(sidx);
        _nmod_vec_clear(valres);
        _nmod_vec_clear(resc);
        _nmod_vec_clear(wpow);
    }

    return GR_SUCCESS;
}

#endif

int
_gr_poly_resultant_multipoint(gr_ptr res, gr_srcptr A, slong lenA,
                              gr_srcptr B, slong lenB, gr_ctx_t ctx)
{
    const gr_poly_struct * Ax = A;
    const gr_poly_struct * Bx = B;
    gr_poly_struct * resx = res;
    gr_ctx_struct * cctx;
    nmod_geometric_progression_t G;
    nmod_t mod;
    nn_ptr valA, valB, valres, resc, tmp, spow, w, f;
    slong i, j, k, blenA, blenB, npoints, maxlen, batch, nblocks, block, nworkers;
    ulong r, q, t;
    _geom_eval_args_t eargs;
    _multipoint_res_args_t rargs;

    if (ctx->which_ring != GR_CTX_GR_POLY)
        return GR_UNABLE;

    cctx = POLYNOMIAL_ELEM_CTX(ctx);

    if (cctx->which_ring != GR_CTX_NMOD)
        return GR_UNABLE;

    if (lenB <= 1)
        return _gr_poly_resultant_small(res, A, lenA, B, lenB, ctx);

    /* the degree bound and the specialisation property below both read the
       degrees in y off lenA and lenB */
    FLINT_ASSERT(Ax[lenA - 1].length != 0);
    FLINT_ASSERT(Bx[lenB - 1].length != 0);

    /* blenA - 1 and blenB - 1 bound the degrees in x */
    blenA = 0;
    for (i = 0; i < lenA; i++)
        blenA = FLINT_MAX(blenA, Ax[i].length);

    blenB = 0;
    for (i = 0; i < lenB; i++)
        blenB = FLINT_MAX(blenB, Bx[i].length);

    /* deg_x(res) <= deg_y(B) deg_x(A) + deg_y(A) deg_x(B); note that this
       makes npoints >= maxlen, since lenA, lenB >= 2 */
    npoints = (lenB - 1) * (blenA - 1) + (lenA - 1) * (blenB - 1) + 1;
    maxlen = FLINT_MAX(blenA, blenB);

    mod = NMOD_CTX(cctx);

    /* The precomputations assume a prime modulus. */
    if (gr_ctx_is_field(cctx) != T_TRUE)
        return GR_UNABLE;

#if FLINT_HAVE_FFT_SMALL
    if (_gr_poly_resultant_multipoint_dft(resx, Ax, lenA, Bx, lenB,
            npoints, maxlen, mod, cctx) == GR_SUCCESS)
        return GR_SUCCESS;
#endif

    /* The ratio of the progression is a square, so its order divides
       (n - 1) / 2 and there is nothing to be done when that is less than
       npoints. */
    r = _nmod_find_geometric_ratio(npoints, mod);

    if (r == 0)
        return GR_UNABLE;

    q = nmod_mul(r, r, mod);

    /* The points are handled in blocks of `batch` of them. Evaluating a
       coefficient of length maxlen at fewer than maxlen points saves nothing,
       so that is the smallest useful block; above that, the block is as large
       as the memory target allows. */
    batch = FLINT_MAX(maxlen, RESULTANT_MULTIPOINT_BLOCK_WORDS / (lenA + lenB));
    batch = FLINT_MIN(batch, npoints);
    nblocks = (npoints + batch - 1) / batch;

    nworkers = flint_get_num_threads();

    valA = _nmod_vec_init(lenA * batch);
    valB = _nmod_vec_init(lenB * batch);
    valres = _nmod_vec_init(npoints);
    resc = _nmod_vec_init(npoints);
    tmp = _nmod_vec_init(nworkers * maxlen);
    spow = _nmod_vec_init(maxlen);
    w = _nmod_vec_init(maxlen);
    f = _nmod_vec_init(nworkers * (lenA + lenB));

    /* The evaluation points are 1, q, q^2, ..., q^(npoints-1), where q = r^2.
       The block starting at index k*batch consists of the points s q^i for
       0 <= i < batch, with s = q^(k*batch); the substitution x -> s x on the
       input polynomials turns those into the powers of q that the geometric
       progression provides, so a single precomputation of length batch serves
       for all blocks. */
    _nmod_geometric_progression_init_function(G, r, batch, mod, UWORD(1));

    /* spow[j] = s^j for the current block, starting with s = 1; multiplying
       by w[j] = (q^batch)^j advances s by one block */
    t = nmod_pow_ui(q, batch, mod);
    spow[0] = 1;
    w[0] = 1;
    for (j = 1; j < maxlen; j++)
    {
        spow[j] = 1;
        w[j] = nmod_mul(w[j - 1], t, mod);
    }

    eargs.Ax = Ax; eargs.Bx = Bx;
    eargs.lenA = lenA; eargs.lenB = lenB;
    eargs.valA = valA; eargs.valB = valB;
    eargs.stride = batch; eargs.maxlen = maxlen;
    eargs.spow = spow; eargs.tmpbuf = tmp;
    eargs.G = G; eargs.mod = mod;

    rargs.valA = valA; rargs.valB = valB;
    rargs.lenA = lenA; rargs.lenB = lenB;
    rargs.stride = batch;
    rargs.res = valres; rargs.dres = NULL; rargs.sidx = NULL;
    rargs.step = 1; rargs.npoints = npoints;
    rargs.fbuf = f; rargs.mod = mod;

    for (block = 0; block < nblocks; block++)
    {
        slong blen = FLINT_MIN(batch, npoints - block * batch);

        eargs.blen = blen;
        eargs.nworkers = FLINT_MIN(nworkers, _num_workers((double) (lenA + lenB)
                                * (maxlen + blen) * FLINT_BIT_COUNT(blen)));
        flint_parallel_do(_geom_eval_worker, &eargs, eargs.nworkers, 0, 0);

        rargs.blen = blen;
        rargs.base = block * batch;
        rargs.nworkers = FLINT_MIN(nworkers,
            _num_workers((double) blen * lenA * lenB));
        flint_parallel_do(_multipoint_res_worker, &rargs, rargs.nworkers, 0, 0);

        for (k = 1; k < maxlen; k++)
            spow[k] = nmod_mul(spow[k], w[k], mod);
    }

    nmod_geometric_progression_clear(G);

    _nmod_geometric_progression_init_function(G, r, npoints, mod, UWORD(2));
    _nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(resc, valres, G, npoints, mod);
    nmod_geometric_progression_clear(G);

    gr_poly_fit_length(resx, npoints, cctx);
    _gr_poly_set_length(resx, npoints, cctx);
    _nmod_vec_set((nn_ptr) resx->coeffs, resc, npoints);
    _gr_poly_normalise(resx, cctx);

    _nmod_vec_clear(valA);
    _nmod_vec_clear(valB);
    _nmod_vec_clear(valres);
    _nmod_vec_clear(resc);
    _nmod_vec_clear(tmp);
    _nmod_vec_clear(spow);
    _nmod_vec_clear(w);
    _nmod_vec_clear(f);

    return GR_SUCCESS;
}

int
gr_poly_resultant_multipoint(gr_ptr r, const gr_poly_t f,
                             const gr_poly_t g, gr_ctx_t ctx)
{
    slong len1 = f->length;
    slong len2 = g->length;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (len1 == 0 || len2 == 0)
    {
        return gr_zero(r, ctx);
    }

    if (gr_is_zero(GR_ENTRY(f->coeffs, len1 - 1, sz), ctx) != T_FALSE ||
        gr_is_zero(GR_ENTRY(g->coeffs, len2 - 1, sz), ctx) != T_FALSE)
    {
        return GR_UNABLE;
    }

    if (len1 >= len2)
    {
        status |= _gr_poly_resultant_multipoint(r, f->coeffs, len1, g->coeffs, len2, ctx);
    }
    else
    {
        status |= _gr_poly_resultant_multipoint(r, g->coeffs, len2, f->coeffs, len1, ctx);

        if (((len1 | len2) & 1) == 0)
            status |= gr_neg(r, r, ctx);
    }

    return status;
}
