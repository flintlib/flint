/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef PADIC_RADIX_IMPL_H
#define PADIC_RADIX_IMPL_H

#include "longlong.h"
#include "mpn_extras.h"
#include "padic_radix.h"

/*
    Packed limb-aligned form, used internally by the block kernels (matrix
    and polynomial multiplication).

    An element x = u p^v (u a unit) is stored relative to the limb frame
    w = floor(v / e) as the L-limb residue U = u p^(v - e w) mod B^L, so that
    x = U B^w modulo B^(w + L). All products then land at whole-limb offsets
    and every output is accumulated in a single set of slot accumulators,
    without digit shifts. U is nonnegative: a negative (exact) unit is stored
    as its complement modulo B^L, which is correct as long as L covers every
    limb that any output needs.

    The precision bookkeeping uses four exponents per element, chosen so that
    the per-output reductions are branch-free min-plus products:

        vv  valuation for the minimum valuation of the nonzero terms
            (PADIC_RADIX_PACK_BIG for a zero);
        ve  valuation in the error term of a product: v for a nonzero x, N
            for an inexact zero 0 + O(p^N) (which lies in p^N Z_p), BIG for
            an exact zero (which annihilates exactly);
        NN  absolute precision N (BIG if exact);
        ww  frame floor(v / e) (BIG for a zero).

    The error term of a product a b is min(NN_a + ve_b, ve_a + NN_b), as in
    padic_radix_mul, and a sum of terms has precision the minimum of these.
    The kernels are only used when every exponent has absolute value at most
    PADIC_RADIX_PACK_MODERATE, so that sums of two genuine exponents stay
    below PADIC_RADIX_PACK_NONE <= BIG + (any genuine exponent) and sums never
    overflow; a minimum >= PACK_NONE means that there is no such term.
*/

/*
    In contexts without finite precision, results are limited by the
    size_limit of the context, in digits (unlimited by default): exact
    results beyond it give GR_UNABLE, and inexact results keep at most
    size_limit digits of relative precision. In contexts with a finite
    precision, intermediate results are bounded in terms of the working
    precision instead.
*/
static inline int
_padic_radix_exact_digits_ok(slong digits, gr_ctx_t ctx)
{
    return digits < 0 || (ulong) digits <= ctx->size_limit;
}

/* in a context without finite precision, the absolute precision N of an
   inexact result of valuation v, lowered (soundly) to at most size_limit
   digits of relative precision */
static inline slong
_padic_radix_cap_prec(slong N, slong v, gr_ctx_t ctx)
{
    if (N > v && (ulong) (N - v) > ctx->size_limit)
        return v + (slong) ctx->size_limit;
    return N;
}

#define PADIC_RADIX_PACK_BIG      (WORD_MAX / 4)
#define PADIC_RADIX_PACK_MODERATE (WORD_MAX / 16)
#define PADIC_RADIX_PACK_NONE     (3 * (WORD_MAX / 16))

typedef struct
{
    slong vv;
    slong ve;
    slong NN;
    slong ww;
}
padic_radix_pack_meta_struct;

/* exponent data of an element; returns 0 if an exponent is not moderate */
static inline int
_padic_radix_pack_meta(padic_radix_pack_meta_struct * m, const padic_radix_struct * x,
    const radix_struct * radix)
{
    if (x->N != PADIC_RADIX_EXACT && FLINT_ABS(x->N) > PADIC_RADIX_PACK_MODERATE)
        return 0;

    if (x->u.size != 0)
    {
        slong v = x->v;
        if (FLINT_ABS(v) > PADIC_RADIX_PACK_MODERATE)
            return 0;
        m->vv = v;
        m->ve = v;
        m->NN = (x->N == PADIC_RADIX_EXACT) ? PADIC_RADIX_PACK_BIG : x->N;
        m->ww = (v >= 0) ? (slong) radix_div_exp((ulong) v, radix)
                         : -(slong) radix_div_exp((ulong) (-v) + radix->exp - 1, radix);
    }
    else if (x->N != PADIC_RADIX_EXACT)
    {
        m->vv = PADIC_RADIX_PACK_BIG;
        m->ve = x->N;
        m->NN = x->N;
        m->ww = PADIC_RADIX_PACK_BIG;
    }
    else
    {
        m->vv = PADIC_RADIX_PACK_BIG;
        m->ve = PADIC_RADIX_PACK_BIG;
        m->NN = PADIC_RADIX_PACK_BIG;
        m->ww = PADIC_RADIX_PACK_BIG;
    }

    return 1;
}

/* U[0, L) = u p^(v - e w) mod B^L (complemented if u < 0); returns the
   number of limbs up to the highest nonzero one (0 for a zero element) */
static inline slong
_padic_radix_pack_unit(nn_ptr U, const padic_radix_struct * x, slong w, slong L,
    const radix_struct * radix)
{
    slong n = FLINT_ABS(x->u.size), s, m;
    ulong cy;

    if (n == 0)
    {
        flint_mpn_zero(U, L);
        return 0;
    }

    s = x->v - (slong) radix->exp * w;
    m = FLINT_MIN(n, L);

    if (s == 0)
    {
        flint_mpn_copyi(U, x->u.d, m);
        if (m < L)
            flint_mpn_zero(U + m, L - m);
    }
    else
    {
        cy = radix_mul_1(U, x->u.d, m, radix->bpow[s], radix);
        if (m < L)
        {
            U[m] = cy;
            if (m + 1 < L)
                flint_mpn_zero(U + m + 1, L - m - 1);
        }
    }

    if (x->u.size < 0)
        radix_neg(U, U, L, radix);

    while (L > 0 && U[L - 1] == 0)
        L--;
    return L;
}

/*
    acc[3 o + 3 i, ...] += (a[0, an) * b[0, bn)) mod B^(nslots - o), as an
    unnormalised sum of limb products (three words per slot).
*/
static inline void
_padic_radix_pack_addmul(nn_ptr acc, nn_srcptr a, slong an, nn_srcptr b, slong bn,
    slong o, slong nslots)
{
    slong m, mmax, ii;

    mmax = an + bn - 2;
    if (mmax > nslots - 1 - o)
        mmax = nslots - 1 - o;

    for (m = 0; m <= mmax; m++)
    {
        slong lo = (m >= bn) ? (m - bn + 1) : 0;
        slong hi = (m < an) ? m : (an - 1);
        ulong c0 = acc[3 * (o + m)];
        ulong c1 = acc[3 * (o + m) + 1];
        ulong c2 = acc[3 * (o + m) + 2];

        for (ii = lo; ii <= hi; ii++)
        {
            ulong ph, pl;
            umul_ppmm(ph, pl, a[ii], b[m - ii]);
            add_sssaaaaaa(c2, c1, c0, c2, c1, c0, 0, ph, pl);
        }

        acc[3 * (o + m)] = c0;
        acc[3 * (o + m) + 1] = c1;
        acc[3 * (o + m) + 2] = c2;
    }
}

/* normalise slot accumulators into radix-B digits d (room for nslots + 3
   limbs); returns the normalised length */
static inline slong
_padic_radix_pack_normalise(nn_ptr d, nn_srcptr acc, slong nslots, const radix_struct * radix)
{
    nmod_t B = radix->B;
    ulong cy[3] = { 0, 0, 0 };
    slong s, outn = 0;

    if (B.norm == 0)
    {
        for (s = 0; s < nslots; s++)
        {
            add_sssaaaaaa(cy[2], cy[1], cy[0], cy[2], cy[1], cy[0],
                acc[3 * s + 2], acc[3 * s + 1], acc[3 * s]);
            d[outn++] = flint_mpn_divrem_3_1_preinv_norm(cy, cy, B.n, B.ninv);
        }
        while (cy[0] | cy[1] | cy[2])
            d[outn++] = flint_mpn_divrem_3_1_preinv_norm(cy, cy, B.n, B.ninv);
    }
    else
    {
        for (s = 0; s < nslots; s++)
        {
            add_sssaaaaaa(cy[2], cy[1], cy[0], cy[2], cy[1], cy[0],
                acc[3 * s + 2], acc[3 * s + 1], acc[3 * s]);
            d[outn++] = flint_mpn_divrem_3_1_preinv_unnorm(cy, cy, B.n, B.ninv, B.norm);
        }
        while (cy[0] | cy[1] | cy[2])
            d[outn++] = flint_mpn_divrem_3_1_preinv_unnorm(cy, cy, B.n, B.ninv, B.norm);
    }

    while (outn > 0 && d[outn - 1] == 0)
        outn--;
    return outn;
}

/* whether the nonzero limb x is divisible by p */
static inline int
_padic_radix_pack_divisible_p(ulong x, const radix_struct * radix)
{
    if (radix->b.n == 2)
        return !(x & 1);
    return n_divisible_odd_gm(x, radix->bpow_oddinv[1].a, radix->bpow_oddinv[1].b);
}

/*
    Finish an output c whose unit field holds the normalised value V of the
    accumulated sum (digits relative to p^(e W)), with precision N computed
    by _padic_radix_pack_prec from vmin. When the lowest digit of V is a
    unit (the common case: no cancellation in the lowest digit) the result
    is already canonical with v = e W = vmin, the target precision is at
    least N, and it only remains to truncate to N - v digits.
*/
static inline int
_padic_radix_pack_finish(padic_radix_struct * c, slong W, slong N, gr_ctx_t ctx)
{
    const radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong v = (slong) radix->exp * W;

    c->v = v;
    c->N = N;

    if (c->u.size > 0 && c->u.d[0] != 0 && !_padic_radix_pack_divisible_p(c->u.d[0], radix)
        && N != PADIC_RADIX_EXACT && N > v)
    {
        slong k = N - v, limbs, rem, n = c->u.size;

        limbs = radix_div_exp(k, radix);
        if (n > limbs)
        {
            rem = k - limbs * (slong) radix->exp;
            if (rem == 0)
                n = limbs;
            else
            {
                n = limbs + 1;
                c->u.d[limbs] = n_rem_precomp(c->u.d[limbs], radix->bpow[rem], radix->bpow_div + rem);
            }
            while (n > 0 && c->u.d[n - 1] == 0)
                n--;
            c->u.size = n;
        }
        return GR_SUCCESS;
    }

    return _padic_radix_finalize(c, ctx);
}

/* the precision of a sum of terms with minimum valuation vmin and error
   term Nt, capped by the context; as in padic_radix_dot_strided_delayed */
static inline slong
_padic_radix_pack_prec(slong vmin, slong Nt, gr_ctx_t ctx)
{
    slong N = Nt;
    if (N >= PADIC_RADIX_PACK_NONE)
        N = PADIC_RADIX_EXACT;
    N = FLINT_MIN(N, PADIC_RADIX_CTX_PREC_ABS(ctx));
    if (PADIC_RADIX_CTX_PREC_REL(ctx) != PADIC_RADIX_PREC_INF && vmin < PADIC_RADIX_PACK_NONE)
        N = FLINT_MIN(N, vmin + PADIC_RADIX_CTX_PREC_REL(ctx));
    return N;
}

/* precision marker for an output whose terms are all exact: such an
   output is computed separately with padic_radix_dot_strided_delayed,
   which keeps the exact sum exact when it fits in the context precision */
#define PADIC_RADIX_PACK_EXACT_SUM (WORD_MIN + 1)

/* min-plus convolutions of exponent sequences (minplus.c):
   res[k - klo] = min_i (a_i + b_(k-i)), exactly or as a lower bound */
void _padic_radix_minplus_conv_exact(slong * res, const slong * a, slong alen,
    const slong * b, slong blen, slong klo, slong khi);
void _padic_radix_minplus_conv_lower(slong * res, const slong * a, slong alen,
    const slong * b, slong blen, slong klo, slong khi);
int _padic_radix_minplus_conv(slong * res, const slong * a, slong alen,
    const slong * b, slong blen, slong klo, slong khi, int exact);

#endif
