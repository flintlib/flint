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
   M-th roots of unity, writing the values to vs in the order of increasing
   powers. The transform enumerates them in bit-reversed order, so rev[] maps
   the power to the position holding its value. plen is at most M, so the
   coefficients need no reduction modulo x^M - 1 first. */
static void
_dft_evaluate(nn_ptr vs, nn_srcptr poly, slong plen, nn_srcptr spow,
              const ulong * rev, slong M, ulong m,
              double * dbuf, sd_fft_ctx_t Q, nmod_t mod)
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
        vs[k] = _sd_fft_get_nmod(dbuf[rev[k]], mod.n);
}

/* res_y(A, B) by evaluation at roots of unity. Returns GR_UNABLE when the
   modulus does not support a transform long enough for the number of points. */
static int
_gr_poly_resultant_multipoint_dft(gr_poly_struct * resx,
                                  const gr_poly_struct * Ax, slong lenA,
                                  const gr_poly_struct * Bx, slong lenB,
                                  slong npoints, slong maxlen, nmod_t mod,
                                  gr_ctx_t cctx)
{
    sd_fft_ctx_t Q;
    nmod_geometric_progression_t G;
    ulong * rev;
    double * dbuf;
    nn_ptr valA, valB, valres, resc, spow, wpow, f, g;
    slong i, j, k, s, t, N, M, nblocks;
    ulong depth, m, w, r;

    depth = FLINT_MAX(n_clog2((ulong) npoints), 4);
    N = n_pow2(depth);

    if (!fft_small_mulmod_satisfies_bounds(mod.n))
        return GR_UNABLE;

    /* A primitive N-th root of unity must exist, and so must a square root r
       of it. That r is of order 2N, since its order divides 2N but not N, so
       this needs one power of two more than the transform itself does. */
    if (n_trailing_zeros(mod.n - 1) < depth + 1)
        return GR_UNABLE;

    sd_fft_ctx_init_prime(Q, mod.n);
    sd_fft_ctx_fit_depth(Q, depth);

    /* the transform evaluates at powers of w; the geometric interpolation
       below works with the points q^i for q = r^2, so r is a square root of w,
       which the test above makes sure exists */
    w = _sd_fft_get_nmod(sd_fft_ctx_w(Q, N / 2), mod.n);
    r = n_sqrtmod(w, mod.n);

    FLINT_ASSERT(r != 0);

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
    nblocks = N / M;

    rev = flint_malloc(M * sizeof(ulong));
    for (k = 0; k < M; k++)
        rev[k] = sd_fft_ctx_trunc_index(m, n_revbin(k, m));

    dbuf = flint_aligned_alloc(32, FLINT_MAX(32, M * sizeof(double)));
    valA = _nmod_vec_init(lenA * M);
    valB = _nmod_vec_init(lenB * M);
    valres = _nmod_vec_init(npoints);
    resc = _nmod_vec_init(npoints);
    spow = _nmod_vec_init(maxlen);
    wpow = _nmod_vec_init(maxlen);
    f = _nmod_vec_init(lenA + lenB);
    g = f + lenA;

    /* spow[k] = (w^t)^k for the current block; multiplying by wpow[k] = w^k
       advances it by one block */
    for (k = 0; k < maxlen; k++)
        spow[k] = 1;

    wpow[0] = 1;
    for (k = 1; k < maxlen; k++)
        wpow[k] = nmod_mul(wpow[k - 1], w, mod);

    for (t = 0; t < nblocks; t++)
    {
        for (i = 0; i < lenA; i++)
            _dft_evaluate(valA + i * M, (nn_srcptr) Ax[i].coeffs,
                Ax[i].length, spow, rev, M, m, dbuf, Q, mod);

        for (i = 0; i < lenB; i++)
            _dft_evaluate(valB + i * M, (nn_srcptr) Bx[i].coeffs,
                Bx[i].length, spow, rev, M, m, dbuf, Q, mod);

        for (s = 0; s < M && (j = t + nblocks * s) < npoints; s++)
        {
            for (i = 0; i < lenA; i++)
                f[i] = valA[i * M + s];

            for (i = 0; i < lenB; i++)
                g[i] = valB[i * M + s];

            _resultant_at_point(valres + j, f, lenA, g, lenB, mod);
        }

        for (k = 1; k < maxlen; k++)
            spow[k] = nmod_mul(spow[k], wpow[k], mod);
    }

    _nmod_geometric_progression_init_function(G, r, npoints, mod, UWORD(2));
    _nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(resc, valres, G, npoints, mod);
    nmod_geometric_progression_clear(G);

    gr_poly_fit_length(resx, npoints, cctx);
    _gr_poly_set_length(resx, npoints, cctx);
    _nmod_vec_set((nn_ptr) resx->coeffs, resc, npoints);
    _gr_poly_normalise(resx, cctx);

    sd_fft_ctx_clear(Q);
    flint_aligned_free(dbuf);
    flint_free(rev);
    _nmod_vec_clear(valA);
    _nmod_vec_clear(valB);
    _nmod_vec_clear(valres);
    _nmod_vec_clear(resc);
    _nmod_vec_clear(spow);
    _nmod_vec_clear(wpow);
    _nmod_vec_clear(f);

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
    nn_ptr valA, valB, valres, resc, tmp, spow, w, f, g;
    slong i, j, k, blenA, blenB, npoints, maxlen, batch, nblocks, block;
    ulong r, q, t;

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

    valA = _nmod_vec_init(lenA * batch);
    valB = _nmod_vec_init(lenB * batch);
    valres = _nmod_vec_init(npoints);
    resc = _nmod_vec_init(npoints);
    tmp = _nmod_vec_init(maxlen);
    spow = _nmod_vec_init(maxlen);
    w = _nmod_vec_init(maxlen);
    f = _nmod_vec_init(lenA + lenB);
    g = f + lenA;

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

    for (block = 0; block < nblocks; block++)
    {
        slong blen = FLINT_MIN(batch, npoints - block * batch);

        for (i = 0; i < lenA; i++)
        {
            for (k = 0; k < Ax[i].length; k++)
                tmp[k] = nmod_mul(((nn_srcptr) Ax[i].coeffs)[k], spow[k], mod);

            _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(valA + i * batch,
                tmp, Ax[i].length, G, blen, mod);
        }

        for (i = 0; i < lenB; i++)
        {
            for (k = 0; k < Bx[i].length; k++)
                tmp[k] = nmod_mul(((nn_srcptr) Bx[i].coeffs)[k], spow[k], mod);

            _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(valB + i * batch,
                tmp, Bx[i].length, G, blen, mod);
        }

        for (j = 0; j < blen; j++)
        {
            for (i = 0; i < lenA; i++)
                f[i] = valA[i * batch + j];

            for (i = 0; i < lenB; i++)
                g[i] = valB[i * batch + j];

            _resultant_at_point(valres + block * batch + j, f, lenA, g, lenB, mod);
        }

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
