/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_poly.h"
#include "nmod_mat.h"
#include "mpn_mod.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_poly.h"
#include "padic_radix.h"
#include "impl.h"

/*
    Matrix-vector products with a matrix packed once in the limb-aligned
    form (impl.h), for algorithms that multiply the same matrix by many
    vectors, such as the Berkowitz algorithm. Each output is computed with
    the semantics of padic_radix_dot_strided_delayed.
*/

typedef struct
{
    slong r;
    slong c;
    slong L;            /* limbs per residue */
    slong * vv;
    slong * ve;
    slong * NN;
    slong * ww;
    slong * sz;
    nn_ptr U;
}
_packed_mat_struct;

/* pack A (all entries) with residues of L limbs; returns 0 if an exponent
   is not moderate */
static int
_packed_mat_init(_packed_mat_struct * P, const gr_mat_t A, slong L, gr_ctx_t ctx)
{
    const radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    padic_radix_pack_meta_struct meta;
    slong i, j, q, len = A->r * A->c;

    P->r = A->r;
    P->c = A->c;
    P->L = L;
    P->vv = flint_malloc(sizeof(slong) * 5 * FLINT_MAX(len, 1));
    P->ve = P->vv + len;
    P->NN = P->ve + len;
    P->ww = P->NN + len;
    P->sz = P->ww + len;
    P->U = flint_malloc(sizeof(ulong) * L * FLINT_MAX(len, 1));

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
        {
            const padic_radix_struct * x = gr_mat_entry_srcptr(A, i, j, ctx);
            q = i * A->c + j;
            if (!_padic_radix_pack_meta(&meta, x, radix))
                return 0;
            P->vv[q] = meta.vv; P->ve[q] = meta.ve; P->NN[q] = meta.NN; P->ww[q] = meta.ww;
            P->sz[q] = (meta.vv >= PADIC_RADIX_PACK_NONE) ? 0 :
                _padic_radix_pack_unit(P->U + L * q, x, meta.ww, L, radix);
            if (P->sz[q] == 0)
                flint_mpn_zero(P->U + L * q, L);
        }

    return 1;
}

static void
_packed_mat_clear(_packed_mat_struct * P)
{
    flint_free(P->vv);
    flint_free(P->U);
}

/*
    res[i - r0] = sum_{j < len} P_ij vec[j] for r0 <= i < r1. res must not
    alias vec. Returns GR_UNABLE (res untouched) if some output needs more
    than P->L limbs or an exponent of vec is not moderate.
*/
static int
_packed_mat_vec(padic_radix_struct * res, const _packed_mat_struct * P, slong r0, slong r1,
    const padic_radix_struct * vec, slong len, gr_ctx_t ctx)
{
    const radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong e = radix->exp, L = P->L, i, j;
    slong * xv, * xe, * xN, * xw, * xs, * oN, * oW, * os;
    nn_ptr XU, acc;
    padic_radix_pack_meta_struct meta;
    int status = GR_SUCCESS;
    TMP_INIT;

    TMP_START;
    xv = TMP_ALLOC(sizeof(slong) * (5 * len + 3 * (r1 - r0)));
    xe = xv + len; xN = xe + len; xw = xN + len; xs = xw + len;
    oN = xs + len; oW = oN + (r1 - r0); os = oW + (r1 - r0);

    for (j = 0; j < len; j++)
    {
        if (!_padic_radix_pack_meta(&meta, vec + j, radix))
            goto unable;
        xv[j] = meta.vv; xe[j] = meta.ve; xN[j] = meta.NN; xw[j] = meta.ww;
    }

    for (i = r0; i < r1; i++)
    {
        const slong * avv = P->vv + i * P->c, * ave = P->ve + i * P->c;
        const slong * aNN = P->NN + i * P->c, * aww = P->ww + i * P->c;
        slong vmin = 2 * PADIC_RADIX_PACK_BIG, Nt = 2 * PADIC_RADIX_PACK_BIG;
        slong W = 2 * PADIC_RADIX_PACK_BIG, N, s;

        for (j = 0; j < len; j++)
        {
            slong x1 = avv[j] + xv[j];
            slong e1 = aNN[j] + xe[j];
            slong e2 = ave[j] + xN[j];
            slong f = aww[j] + xw[j];
            vmin = FLINT_MIN(vmin, x1);
            e1 = FLINT_MIN(e1, e2);
            Nt = FLINT_MIN(Nt, e1);
            W = FLINT_MIN(W, f);
        }

        if (vmin >= PADIC_RADIX_PACK_NONE)
        {
            N = (Nt >= PADIC_RADIX_PACK_NONE) ? PADIC_RADIX_EXACT : Nt;
            s = 0;
        }
        else
        {
            N = _padic_radix_pack_prec(vmin, Nt, ctx);
            s = (N <= vmin) ? 0 : (slong) radix_div_exp((ulong) (N - e * W) + e - 1, radix);
            if (s > L)
                goto unable;
        }

        oN[i - r0] = N;
        oW[i - r0] = W;
        os[i - r0] = s;
    }

    XU = TMP_ALLOC(sizeof(ulong) * L * FLINT_MAX(len, 1));
    for (j = 0; j < len; j++)
    {
        xs[j] = (xv[j] >= PADIC_RADIX_PACK_NONE) ? 0 :
            _padic_radix_pack_unit(XU + L * j, vec + j, xw[j], L, radix);
        if (xs[j] == 0)
            flint_mpn_zero(XU + L * j, L);
    }

    acc = TMP_ALLOC(sizeof(ulong) * 3 * L);

    for (i = r0; i < r1; i++)
    {
        const slong * aww = P->ww + i * P->c, * asz = P->sz + i * P->c;
        nn_srcptr ua = P->U + L * i * P->c;
        slong N = oN[i - r0], W = oW[i - r0], nslots = os[i - r0];
        padic_radix_struct * c = res + (i - r0);
        nn_ptr d;

        if (nslots == 0)
        {
            radix_integer_zero(&c->u, radix);
            c->v = 0;
            c->N = N;
            status |= _padic_radix_finalize(c, ctx);
            continue;
        }

        flint_mpn_zero(acc, 3 * nslots);

        if (L == 1)
        {
            ulong c0 = 0, c1 = 0, c2 = 0, ph, pl;
            for (j = 0; j < len; j++)
            {
                if (aww[j] + xw[j] == W)
                {
                    umul_ppmm(ph, pl, ua[j], XU[j]);
                    add_sssaaaaaa(c2, c1, c0, c2, c1, c0, 0, ph, pl);
                }
            }
            acc[0] = c0; acc[1] = c1; acc[2] = c2;
        }
        else if (L == 2)
        {
            ulong c0 = 0, c1 = 0, c2 = 0, d0 = 0, d1 = 0, d2 = 0, ph, pl;
            for (j = 0; j < len; j++)
            {
                slong o = aww[j] + xw[j] - W;
                ulong a0 = ua[2 * j], a1 = ua[2 * j + 1], b0 = XU[2 * j], b1 = XU[2 * j + 1];

                if (o == 0)
                {
                    umul_ppmm(ph, pl, a0, b0);
                    add_sssaaaaaa(c2, c1, c0, c2, c1, c0, 0, ph, pl);
                    if (nslots == 2)
                    {
                        umul_ppmm(ph, pl, a0, b1);
                        add_sssaaaaaa(d2, d1, d0, d2, d1, d0, 0, ph, pl);
                        umul_ppmm(ph, pl, a1, b0);
                        add_sssaaaaaa(d2, d1, d0, d2, d1, d0, 0, ph, pl);
                    }
                }
                else if (o == 1 && nslots == 2)
                {
                    umul_ppmm(ph, pl, a0, b0);
                    add_sssaaaaaa(d2, d1, d0, d2, d1, d0, 0, ph, pl);
                }
            }
            acc[0] = c0; acc[1] = c1; acc[2] = c2;
            acc[3] = d0; acc[4] = d1; acc[5] = d2;
        }
        else
        {
            for (j = 0; j < len; j++)
            {
                slong an = asz[j], bn = xs[j], o;

                o = aww[j] + xw[j] - W;
                if (o >= nslots || an == 0 || bn == 0)
                    continue;

                _padic_radix_pack_addmul(acc, ua + L * j, an, XU + L * j, bn, o, nslots);
            }
        }

        d = radix_integer_fit_limbs(&c->u, nslots + 3, radix);
        c->u.size = _padic_radix_pack_normalise(d, acc, nslots, radix);
        status |= _padic_radix_pack_finish(c, W, N, ctx);
    }

    TMP_END;
    return status;

unable:
    TMP_END;
    return GR_UNABLE;
}

/*
    The Berkowitz algorithm (as _gr_mat_charpoly_berkowitz), with the
    matrix-vector products a_k = M_t a_(k-1) done by _packed_mat_vec on a
    matrix packed once. Falls back to the generic algorithm for contexts
    without a relative precision (which bounds the windows a priori), for
    matrices with exact nonzero entries, and for extreme exponents.
*/
int
_padic_radix_mat_charpoly_berkowitz(padic_radix_struct * cp, const gr_mat_t mat, gr_ctx_t ctx)
{
    const radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong n = mat->r, e = radix->exp, r = PADIC_RADIX_CTX_PREC_REL(ctx), L, Lmax;
    slong i, k, t;
    padic_radix_struct * a, * A;
    _packed_mat_struct P;
    int status = GR_SUCCESS;

    if (n <= 5 || r == PADIC_RADIX_PREC_INF || r > PADIC_RADIX_PACK_MODERATE)
        return _gr_mat_charpoly_berkowitz(cp, mat, ctx);

    /* exact nonzero entries can give dot products with only exact terms,
       which keep exact sums exact: leave those to the generic algorithm */
    for (i = 0; i < n; i++)
        for (k = 0; k < n; k++)
        {
            const padic_radix_struct * x = gr_mat_entry_srcptr(mat, i, k, ctx);
            if (x->u.size != 0 && x->N == PADIC_RADIX_EXACT)
                return _gr_mat_charpoly_berkowitz(cp, mat, ctx);
        }

    /* residues of ceil(r / e) limbs suffice when the lowest term of an
       output is aligned with its frame (e.g. all valuations 0); in general
       an output window is at most
       r + 2e - 2 digits (see impl.h), and the matrix is repacked with that
       many limbs the first time a product needs them */
    L = FLINT_MAX(1, (slong) radix_div_exp((ulong) (r + e - 1), radix));
    Lmax = (slong) radix_div_exp((ulong) (r + 3 * e - 3), radix);

    if (!_packed_mat_init(&P, mat, L, ctx))
    {
        _packed_mat_clear(&P);
        return _gr_mat_charpoly_berkowitz(cp, mat, ctx);
    }

#define MAT(ii, jj) gr_mat_entry_srcptr(mat, ii, jj, ctx)

    a = gr_heap_init_vec(n * n, ctx);
    A = a + (n - 1) * n;

    status |= _gr_vec_zero(cp, n + 1, ctx);
    status |= padic_radix_neg(cp, MAT(0, 0), ctx);

    for (t = 1; t < n && status == GR_SUCCESS; t++)
    {
        for (i = 0; i <= t; i++)
            status |= padic_radix_set(a + i, MAT(i, t), ctx);

        status |= padic_radix_set(A, MAT(t, t), ctx);

        for (k = 1; k < t; k++)
        {
            int s = _packed_mat_vec(a + k * n, &P, 0, t + 1, a + (k - 1) * n, t + 1, ctx);
            if (s == GR_UNABLE && P.L < Lmax)
            {
                _packed_mat_clear(&P);
                if (!_packed_mat_init(&P, mat, Lmax, ctx))
                    flint_throw(FLINT_ERROR, "(%s)\n", __func__);
                s = _packed_mat_vec(a + k * n, &P, 0, t + 1, a + (k - 1) * n, t + 1, ctx);
            }
            if (s == GR_UNABLE)
            {
                s = GR_SUCCESS;
                for (i = 0; i <= t; i++)
                    s |= _gr_vec_dot(a + k * n + i, NULL, 0, MAT(i, 0), a + (k - 1) * n, t + 1, ctx);
            }
            status |= s;
            status |= padic_radix_set(A + k, a + k * n + t, ctx);
        }

        {
            int s = _packed_mat_vec(A + t, &P, t, t + 1, a + (t - 1) * n, t + 1, ctx);
            if (s == GR_UNABLE && P.L < Lmax)
            {
                _packed_mat_clear(&P);
                if (!_packed_mat_init(&P, mat, Lmax, ctx))
                    flint_throw(FLINT_ERROR, "(%s)\n", __func__);
                s = _packed_mat_vec(A + t, &P, t, t + 1, a + (t - 1) * n, t + 1, ctx);
            }
            if (s == GR_UNABLE)
                s = _gr_vec_dot(A + t, NULL, 0, MAT(t, 0), a + (t - 1) * n, t + 1, ctx);
            status |= s;
        }

        for (k = 0; k <= t; k++)
        {
            status |= _gr_vec_dot_rev(cp + k, cp + k, 1, A, cp, k, ctx);
            status |= padic_radix_sub(cp + k, cp + k, A + k, ctx);
        }
    }

    for (i = n; i > 0; i--)
        padic_radix_swap(cp + i, cp + i - 1, ctx);

    status |= padic_radix_one(cp, ctx);
    status |= _gr_poly_reverse(cp, cp, n + 1, n + 1, ctx);

    gr_heap_clear_vec(a, n * n, ctx);
    _packed_mat_clear(&P);

#undef MAT

    return status;
}

/*
    Reduction to Hessenberg form by similarity transformations, with the
    pivot of minimal valuation in each column (so that all multipliers have
    nonnegative valuation and the transformations are in GL_n(Z_p)), and
    the characteristic polynomial of the Hessenberg matrix by the usual
    O(n^3) recurrence. Returns GR_UNABLE if a column below the subdiagonal
    has no nonzero entry but some inexact zero.
*/
int
padic_radix_mat_charpoly_hessenberg(padic_radix_struct * cp, const gr_mat_t mat, gr_ctx_t ctx)
{
    slong n = mat->r, m, i, j, r;
    gr_mat_t H;
    padic_radix_t inv, f;
    int status = GR_SUCCESS;

    if (n != mat->c)
        return GR_DOMAIN;

    if (n <= 2)
        return _gr_mat_charpoly_berkowitz(cp, mat, ctx);

    gr_mat_init(H, n, n, ctx);
    status |= gr_mat_set(H, mat, ctx);
    padic_radix_init(inv, ctx);
    padic_radix_init(f, ctx);

#define HE(ii, jj) ((padic_radix_struct *) gr_mat_entry_ptr(H, ii, jj, ctx))

    for (m = 0; m + 2 < n && status == GR_SUCCESS; m++)
    {
        slong bv = WORD_MAX;
        int unknown = 0;

        /* pivot in column m, rows m + 1, ..., n - 1 */
        r = -1;
        for (i = m + 1; i < n; i++)
        {
            const padic_radix_struct * x = HE(i, m);
            if (x->u.size != 0 && x->v < bv)
            {
                bv = x->v;
                r = i;
            }
            else if (x->u.size == 0 && x->N != PADIC_RADIX_EXACT)
                unknown = 1;
        }

        if (r < 0)
        {
            if (unknown)
                status = GR_UNABLE;
            continue;       /* already reduced */
        }

        if (r != m + 1)
        {
            _gr_vec_swap(HE(r, 0), HE(m + 1, 0), n, ctx);
            for (i = 0; i < n; i++)
                padic_radix_swap(HE(i, r), HE(i, m + 1), ctx);
        }

        status |= padic_radix_inv(inv, HE(m + 1, m), ctx);

        for (i = m + 2; i < n; i++)
        {
            padic_radix_struct * x = HE(i, m);

            if (x->u.size == 0 && x->N == PADIC_RADIX_EXACT)
                continue;

            /* row_i -= f row_(m+1), then col_(m+1) += f col_i */
            status |= padic_radix_mul(f, x, inv, ctx);
            status |= _padic_radix_vec_submul_scalar(HE(i, m), HE(m + 1, m), n - m, f, ctx);
            for (j = 0; j < n; j++)
                status |= padic_radix_addmul(HE(j, m + 1), f, HE(j, i), ctx);
        }
    }

#undef HE

    if (status == GR_SUCCESS)
        status |= _gr_mat_charpoly_from_hessenberg(cp, H, ctx);

    padic_radix_clear(inv, ctx);
    padic_radix_clear(f, ctx);
    gr_mat_clear(H, ctx);
    return status;
}

/* Z/PZ with the fastest available representation */
static void
_mod_ctx_init(gr_ctx_t mctx, const fmpz_t P)
{
    if (fmpz_abs_fits_ui(P))
        GR_MUST_SUCCEED(gr_ctx_init_nmod(mctx, fmpz_get_ui(P)));
    else if (fmpz_size(P) <= MPN_MOD_MAX_LIMBS)
        GR_MUST_SUCCEED(gr_ctx_init_mpn_mod(mctx, P));
    else
        gr_ctx_init_fmpz_mod(mctx, P);
}

/* the valuation of the integer x (destroyed), capped at cap; q = p */
static slong
_fmpz_val_capped(fmpz_t x, const fmpz_t q, slong cap)
{
    slong v;

    if (fmpz_is_zero(x))
        return cap;

    if (!COEFF_IS_MPZ(*x))
    {
        ulong a = FLINT_ABS(*x), pp = *q;
        for (v = 0; v < cap && a % pp == 0; v++)
            a /= pp;
        return v;
    }

    v = fmpz_remove(x, x, q);
    return FLINT_MIN(v, cap);
}

/*
    Characteristic polynomial with the optimal (first-order) precision.

    Scale A by p^s so that all entries (and the error balls of the inexact
    entries) are integral, and let M be the exact lift of the scaled matrix
    (inexact zeros lift to 0). For a perturbation E with v(E_ij) >= N_ij
    (E_ij = 0 for exact entries), the coefficients of chi(M + E) are those
    of chi(M) plus the first-order term -tr(B_k E), where
    adj(xI - M) = sum_k B_k x^k, plus terms of order at least two in E,
    whose coefficients (minors of xI - M) are integral: those have valuation
    at least 2 Nmin. Hence c_k(A) = c_k(M) + O(p^(N_k)) with

        N_k = min(2 Nmin, min_ij (N_ij + v((B_k)_ji))),

    which is the optimal precision whenever the first-order term dominates
    (it is the precision of the differential of chi). We compute chi(M)
    modulo p^(2 Nmin) by Berkowitz's algorithm, and the B_k modulo p^Nmin
    (enough for the valuations that matter) by the recurrence
    B_(n-1) = I, B_(k-1) = M B_k + c_k I, which costs n matrix products.
    For an exact matrix chi(M) is computed exactly (multimodular).
*/
static int
_charpoly_lift(padic_radix_struct * cp, const gr_mat_t mat, int * capped, gr_ctx_t ctx)
{
    const radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    ulong p = GR_PADIC_RADIX_CTX(ctx)->p;
    slong n = mat->r, i, j, k, sc = 0, Nmin = WORD_MAX, Wtop;
    fmpz_mat_t M;
    fmpz_poly_t c;
    fmpz_t t, P1, P2;
    slong * Nk, * NN;
    int status = GR_SUCCESS;

    if (n != mat->c)
        return GR_DOMAIN;

    if (n == 0)
        return padic_radix_one(cp, ctx);

    /* scaling and precisions */
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            const padic_radix_struct * x = gr_mat_entry_srcptr(mat, i, j, ctx);
            slong lo;

            if (FLINT_ABS(x->v) > PADIC_RADIX_PACK_MODERATE / (4 * n)
                || (x->N != PADIC_RADIX_EXACT && FLINT_ABS(x->N) > PADIC_RADIX_PACK_MODERATE / (4 * n)))
                return GR_UNABLE;

            if (x->u.size != 0)
                lo = x->v;
            else if (x->N != PADIC_RADIX_EXACT)
                lo = x->N;
            else
                continue;

            sc = FLINT_MAX(sc, -lo);
            if (x->N != PADIC_RADIX_EXACT)
                Nmin = FLINT_MIN(Nmin, x->N);
        }

    if (Nmin != WORD_MAX && Nmin + sc <= 0)
        return GR_UNABLE;       /* nothing known modulo p: no gain */

    /*
        Sizes. With an inexact entry, everything is computed modulo
        p^Wtop (Wtop the largest working precision below), so the lifted
        entries are reduced modulo p^Wtop; the cost is then governed by
        K = Nmin + s, which is large only when the entries have valuations
        (or error terms) spread far beyond their relative precisions R. In
        that case the precision is that of a much smaller computation, and
        the lift is not worth it (UNABLE: the caller uses the algorithms that
        work with the precisions of the entries). For an exact matrix, the
        lift computes the exact characteristic polynomial, of about n D
        digits where D is the largest size of an entry of M: only done when
        the context keeps that many digits (or D is moderate compared with
        the context precision).
    */
    {
        slong R = 0, D = 0, Du = 0, vhi = 0, vlo = WORD_MAX;

        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
            {
                const padic_radix_struct * x = gr_mat_entry_srcptr(mat, i, j, ctx);

                if (x->u.size != 0)
                {
                    vlo = FLINT_MIN(vlo, x->v);
                    vhi = FLINT_MAX(vhi, x->v + sc);
                    Du = FLINT_MAX(Du, (slong) radix->exp * FLINT_ABS(x->u.size));
                    D = FLINT_MAX(D, x->v + sc + (slong) radix->exp * FLINT_ABS(x->u.size));
                    if (x->N != PADIC_RADIX_EXACT)
                        R = FLINT_MAX(R, x->N - x->v);
                }
            }

        if (Nmin != WORD_MAX)
        {
            slong K = Nmin + sc;

            if ((double) K > 64.0 * ((double) R + n) + 1024)
                return GR_UNABLE;
        }
        else
        {
            slong rel = PADIC_RADIX_CTX_PREC_REL(ctx), pabs = PADIC_RADIX_CTX_PREC_ABS(ctx), C;

            if (D > (WORD_MAX / 2) / n || !_padic_radix_exact_digits_ok(n * D, ctx))
                return GR_UNABLE;

            if (rel != PADIC_RADIX_PREC_INF)
                C = rel;
            else if (pabs != PADIC_RADIX_PREC_INF && vlo != WORD_MAX)
                C = pabs - FLINT_MIN(vlo, 0);
            else
                C = WORD_MAX;

            if (C != WORD_MAX && D > 4 * C + 256)
                return GR_UNABLE;

            /* the integer matrix is dominated by powers of p which the
               p-adic representation (and division-free algorithms in it)
               avoids */
            if (vhi > 4 * Du + 256)
                return GR_UNABLE;
        }
    }

    /* the largest working precision: W = 2 K + ((2 K / 4 + 16) << 2) */
    Wtop = (Nmin == WORD_MAX) ? WORD_MAX : 2 * (Nmin + sc) + ((2 * (Nmin + sc) / 4 + 16) << 2);

    fmpz_mat_init(M, n, n);
    fmpz_poly_init(c);
    fmpz_init(t);
    fmpz_init(P1);
    fmpz_init(P2);
    NN = flint_malloc(sizeof(slong) * (n * n + n + 1));
    Nk = NN + n * n;

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            const padic_radix_struct * x = gr_mat_entry_srcptr(mat, i, j, ctx);
            fmpz * m = fmpz_mat_entry(M, i, j);

            NN[i * n + j] = (x->N == PADIC_RADIX_EXACT) ? WORD_MAX : x->N + sc;

            if (x->u.size != 0)
            {
                if (Nmin == WORD_MAX)
                {
                    radix_integer_get_fmpz(m, &x->u, radix);
                    fmpz_ui_pow_ui(t, p, x->v + sc);
                    fmpz_mul(m, m, t);
                }
                else if (x->v + sc < Wtop)
                {
                    /* modulo p^Wtop */
                    radix_integer_t r;
                    radix_integer_init(r, radix);
                    radix_integer_mod_digits(r, &x->u, Wtop - (x->v + sc), radix);
                    radix_integer_get_fmpz(m, r, radix);
                    radix_integer_clear(r, radix);
                    fmpz_ui_pow_ui(t, p, x->v + sc);
                    fmpz_mul(m, m, t);
                }
            }
        }

    if (Nmin == WORD_MAX)
    {
        /* exact matrix */
        fmpz_mat_charpoly(c, M);
        for (k = 0; k < n; k++)
            Nk[k] = WORD_MAX;
    }
    else
    {
        slong K = Nmin + sc;

        fmpz_ui_pow_ui(P1, p, K);
        fmpz_mul(P2, P1, P1);

        /* chi(M) modulo p^(2 Nmin). M is exact, so any algorithm can be run
           at a working precision W above 2 Nmin, with entries known to
           O(p^W): the Hessenberg reduction with valuation pivots (O(n^3)),
           whose precision loss is tracked and usually small, retrying with
           a larger W if some coefficient falls short of 2 Nmin, and finally
           Berkowitz's division-free algorithm (O(n^4), no loss at all). */
        {
            slong W, attempt;
            int ok = 0;

            fmpz_poly_fit_length(c, n + 1);

            for (attempt = 0; attempt < 4 && !ok; attempt++)
            {
                gr_ctx_t F2;
                gr_mat_t M2;
                padic_radix_struct * c2;
                const radix_struct * radix2;
                int s2;

                W = 2 * K + ((attempt < 3) ? (2 * K / 4 + 16) << attempt : 0);
                gr_ctx_init_padic_radix(F2, p, W, W, 0);
                radix2 = PADIC_RADIX_CTX_RADIX(F2);
                gr_mat_init(M2, n, n, F2);
                for (i = 0; i < n; i++)
                    for (j = 0; j < n; j++)
                    {
                        padic_radix_struct * x = gr_mat_entry_ptr(M2, i, j, F2);
                        status |= padic_radix_exact_set_fmpz(x, fmpz_mat_entry(M, i, j), F2);
                        x->N = W;
                        status |= _padic_radix_finalize(x, F2);
                    }
                c2 = gr_heap_init_vec(n + 1, F2);

                if (attempt < 3)
                    s2 = padic_radix_mat_charpoly_hessenberg(c2, M2, F2);
                else
                    s2 = _padic_radix_mat_charpoly_berkowitz(c2, M2, F2);

                if (s2 == GR_SUCCESS)
                {
                    ok = 1;
                    for (k = 0; k < n; k++)
                        if (c2[k].N != PADIC_RADIX_EXACT && c2[k].N < 2 * K)
                            ok = 0;
                }

                if (ok || attempt == 3)
                {
                    for (k = 0; k <= n; k++)
                    {
                        fmpz * ck = c->coeffs + k;
                        if (c2[k].u.size == 0)
                            fmpz_zero(ck);
                        else
                        {
                            radix_integer_get_fmpz(ck, &c2[k].u, radix2);
                            fmpz_ui_pow_ui(t, p, c2[k].v);
                            fmpz_mul(ck, ck, t);
                        }
                        Nk[k] = (k < n && c2[k].N != PADIC_RADIX_EXACT) ? FLINT_MIN(2 * K, c2[k].N) : 2 * K;
                    }
                    if (s2 != GR_SUCCESS)
                        status |= s2;
                    ok = 1;
                }

                gr_heap_clear_vec(c2, n + 1, F2);
                gr_mat_clear(M2, F2);
                gr_ctx_clear(F2);
            }

            _fmpz_poly_set_length(c, n + 1);
        }

        /* the B_k modulo p^Kmod: first with one word, then (if a valuation
           at the cap matters) modulo p^K */
        {
            slong Kmod = FLINT_MIN(K, (slong) radix_limb_exponent(radix)), pass;
            slong * Nfirst = flint_malloc(sizeof(slong) * n);
            fmpz_t pq;

            fmpz_init_set_ui(pq, p);

            for (k = 0; k < n; k++)
                Nfirst[k] = Nk[k];

            for (pass = 0; status == GR_SUCCESS; pass++)
            {
                gr_ctx_t m1;
                gr_mat_t Mm, Bk, T;
                int more = 0;

                fmpz_ui_pow_ui(P1, p, Kmod);
                _mod_ctx_init(m1, P1);
                gr_mat_init(Mm, n, n, m1);
                gr_mat_init(Bk, n, n, m1);
                gr_mat_init(T, n, n, m1);
                for (i = 0; i < n; i++)
                    for (j = 0; j < n; j++)
                        status |= gr_set_fmpz(gr_mat_entry_ptr(Mm, i, j, m1), fmpz_mat_entry(M, i, j), m1);
                status |= gr_mat_one(Bk, m1);

                for (k = n - 1; k >= 0 && status == GR_SUCCESS; k--)
                {
                    slong best = Nfirst[k];
                    int atcap = 0;

                    for (i = 0; i < n; i++)
                        for (j = 0; j < n; j++)
                        {
                            slong Nij = NN[i * n + j], vb;
                            if (Nij == WORD_MAX || Nij >= best)
                                continue;
                            status |= gr_get_fmpz(t, gr_mat_entry_srcptr(Bk, j, i, m1), m1);
                            vb = _fmpz_val_capped(t, pq, Kmod);
                            if (Nij + vb < best)
                            {
                                best = Nij + vb;
                                atcap = (vb == Kmod);
                            }
                        }

                    if (atcap && Kmod < K)
                        more = 1;

                    Nk[k] = best;
                    if (best == Nfirst[k])
                        *capped = 1;

                    if (k > 0)
                    {
                        status |= gr_mat_mul(T, Mm, Bk, m1);
                        gr_mat_swap(T, Bk, m1);
                        fmpz_mod(t, c->coeffs + k, P1);
                        for (i = 0; i < n; i++)
                            status |= gr_add_fmpz(gr_mat_entry_ptr(Bk, i, i, m1), gr_mat_entry_srcptr(Bk, i, i, m1), t, m1);
                    }
                }

                gr_mat_clear(Mm, m1);
                gr_mat_clear(Bk, m1);
                gr_mat_clear(T, m1);
                gr_ctx_clear(m1);

                if (!more)
                    break;
                Kmod = K;
            }

            flint_free(Nfirst);
            fmpz_clear(pq);
        }
    }

    /* c_k(A) = c_k(M) p^(-s (n - k)) */
    for (k = 0; k < n && status == GR_SUCCESS; k++)
    {
        padic_radix_struct * r = cp + k;
        slong shift = sc * (n - k);

        status |= padic_radix_exact_set_fmpz(r, c->coeffs + k, ctx);
        if (Nk[k] != WORD_MAX)
        {
            fmpz_mod(t, c->coeffs + k, P2);
            status |= padic_radix_exact_set_fmpz(r, t, ctx);
        }
        if (r->u.size != 0)
            r->v -= shift;
        else
            r->v = 0;
        r->N = (Nk[k] == WORD_MAX) ? PADIC_RADIX_EXACT : Nk[k] - shift;
        status |= _padic_radix_finalize(r, ctx);
    }
    status |= padic_radix_one(cp + n, ctx);

    flint_free(NN);
    fmpz_mat_clear(M);
    fmpz_poly_clear(c);
    fmpz_clear(t);
    fmpz_clear(P1);
    fmpz_clear(P2);
    return status;
}

int
padic_radix_mat_charpoly_lift(padic_radix_struct * cp, const gr_mat_t mat, gr_ctx_t ctx)
{
    int capped = 0;
    return _charpoly_lift(cp, mat, &capped, ctx);
}

/* whether the matrix of integers M is cyclic modulo p (a random Krylov
   basis has full rank; may say no for a cyclic matrix with small chance) */
static int
_cyclic_mod_p(const gr_mat_t mat, ulong p, gr_ctx_t ctx)
{
    const radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong n = mat->r, i, j, k, tries;
    nmod_mat_t A, K;
    nn_ptr v, w;
    flint_rand_t state;
    int res = 0;

    nmod_mat_init(A, n, n, p);
    nmod_mat_init(K, n, n, p);
    v = flint_malloc(sizeof(ulong) * 2 * n);
    w = v + n;
    flint_rand_init(state);

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            const padic_radix_struct * x = gr_mat_entry_srcptr(mat, i, j, ctx);
            ulong r = 0;
            if (x->u.size != 0 && x->v == 0)
            {
                r = n_mod2_preinv(x->u.d[0], p, n_preinvert_limb(p));
                if (x->u.size < 0)
                    r = nmod_neg(r, A->mod);
            }
            nmod_mat_entry(A, i, j) = r;
        }
    (void) radix;

    for (tries = 0; tries < 2 && !res; tries++)
    {
        for (i = 0; i < n; i++)
            v[i] = n_randint(state, p);

        for (k = 0; k < n; k++)
        {
            for (i = 0; i < n; i++)
                nmod_mat_entry(K, i, k) = v[i];
            for (i = 0; i < n; i++)
            {
                ulong s = 0;
                for (j = 0; j < n; j++)
                    s = nmod_add(s, nmod_mul(nmod_mat_entry(A, i, j), v[j], A->mod), A->mod);
                w[i] = s;
            }
            flint_mpn_copyi(v, w, n);
        }

        res = (nmod_mat_rank(K) == n);
    }

    flint_rand_clear(state);
    flint_free(v);
    nmod_mat_clear(A);
    nmod_mat_clear(K);
    return res;
}

/*
    Berkowitz's division-free algorithm gives precision N for all
    coefficients of a matrix with integral entries known to O(p^N); this is
    optimal when the matrix is cyclic modulo p (then every B_k has a unit
    entry, B_k being a monic polynomial in M of degree < n). Otherwise (or
    for jagged precisions, negative valuations, exact entries), use the
    lifting method with the precision of the differential. With
    PADIC_RADIX_FAST_LINALG, use the Hessenberg reduction.
*/
int
padic_radix_mat_charpoly(padic_radix_struct * cp, const gr_mat_t mat, gr_ctx_t ctx)
{
    slong n = mat->r, i, j, N0 = WORD_MAX;
    int uniform = 1, status, capped = 0;

    if (n <= 1)
        return _gr_mat_charpoly_berkowitz(cp, mat, ctx);

    for (i = 0; i < n && uniform; i++)
        for (j = 0; j < n && uniform; j++)
        {
            const padic_radix_struct * x = gr_mat_entry_srcptr(mat, i, j, ctx);

            if (x->N == PADIC_RADIX_EXACT || (x->u.size != 0 && x->v < 0) || x->N <= 0)
                uniform = 0;
            else if (N0 == WORD_MAX)
                N0 = x->N;
            else if (x->N != N0)
                uniform = 0;
        }

    /* speed rather than optimal precision: the O(n^3) Hessenberg reduction
       (which loses a few digits on structured matrices) */
    if (PADIC_RADIX_CTX_FLAGS(ctx) & PADIC_RADIX_FAST_LINALG)
    {
        status = padic_radix_mat_charpoly_hessenberg(cp, mat, ctx);
        if (status == GR_UNABLE)
            status = _padic_radix_mat_charpoly_berkowitz(cp, mat, ctx);
        return status;
    }

    if (uniform && _cyclic_mod_p(mat, GR_PADIC_RADIX_CTX(ctx)->p, ctx))
    {
        /* the optimal precision is N0 for every coefficient: the Hessenberg
           reduction usually attains it in O(n^3) operations (it does not
           lose any digits on matrices that are generic modulo p), which is
           then certified; Berkowitz's algorithm (O(n^4)) otherwise */
        if (padic_radix_mat_charpoly_hessenberg(cp, mat, ctx) == GR_SUCCESS)
        {
            int ok = 1;
            for (i = 0; i < n && ok; i++)
                ok = (cp[i].N >= N0);
            if (ok)
                return GR_SUCCESS;
        }
        return _padic_radix_mat_charpoly_berkowitz(cp, mat, ctx);
    }

    status = _charpoly_lift(cp, mat, &capped, ctx);
    if (status == GR_UNABLE)
        return _padic_radix_mat_charpoly_berkowitz(cp, mat, ctx);

    /* where the bound 2 Nmin for the terms of higher order was reached, the
       interval arithmetic of Danilevsky's algorithm (similarity transforms,
       good for structured matrices over large p) may do better: keep the
       more precise of the two (valid) results for each coefficient */
    if (status == GR_SUCCESS && capped)
    {
        padic_radix_struct * tmp = gr_heap_init_vec(n + 1, ctx);
        if (_gr_mat_charpoly_danilevsky(tmp, mat, ctx) == GR_SUCCESS)
        {
            for (i = 0; i < n; i++)
                if (cp[i].N != PADIC_RADIX_EXACT && (tmp[i].N == PADIC_RADIX_EXACT || tmp[i].N > cp[i].N))
                    padic_radix_swap(cp + i, tmp + i, ctx);
        }
        gr_heap_clear_vec(tmp, n + 1, ctx);
    }

    return status;
}
