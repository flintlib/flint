/*
    Copyright (C) 2024, 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "mpn_mod.h"
#include "mpn_mod/impl.h"

/*
    Schoolbook polynomial multiplication and squaring with delayed
    reduction: the coefficients of the product are computed exactly, as
    unreduced multiprecision integers of a fixed size slimbs, with no
    modular reductions. This lets callers such as the cyclotomic
    arithmetic in aprcl perform linear operations (and a single reduction
    per coefficient) on the unreduced coefficients.

    The mathematical kernels overlap with the _mpn_dot_rev functions used
    by the classical and Karatsuba mullow code, and the generic branches
    below reuse those functions directly. For nlimbs == 2 (and squaring
    with nlimbs == 3) the coefficients are instead accumulated by fully
    inlined loops holding the accumulator in registers: at the very short
    dot lengths arising here (at most around 10 terms, often 1 to 4), a
    function call per coefficient plus a memory round trip for doubling
    and square terms measurably costs 5-30%.
*/

slong
_mpn_mod_poly_mul_unreduced_slimbs(slong len, gr_ctx_t ctx)
{
    flint_bitcnt_t sbits;
    slong slimbs;

    sbits = 2 * MPN_MOD_CTX_MODULUS_BITS(ctx) + FLINT_BIT_COUNT(len);
    slimbs = (sbits + FLINT_BITS - 1) / FLINT_BITS;

    /* for nlimbs == 3 the 6-limb register kernels beat the 5-limb
       function-call dot products at short lengths */
    if (MPN_MOD_CTX_NLIMBS(ctx) == 3)
        slimbs = FLINT_MAX(slimbs, 6);

    return slimbs;
}

/* Truncated square: low three limbs of (a1 B + a0)^2 with B = 2^FLINT_BITS,
   valid whenever the square fits in three limbs. */
#define _MPN_MOD_SQR_3P2X2(r2, r1, r0, a1, a0)              \
    do {                                                  \
        ulong __q1, __q0;                                 \
        umul_ppmm(r1, r0, a0, a0);                        \
        (r2) = (a1) * (a1);                               \
        umul_ppmm(__q1, __q0, a0, a1);                    \
        add_ssaaaa(r2, r1, r2, r1, __q1, __q0);           \
        add_ssaaaa(r2, r1, r2, r1, __q1, __q0);           \
    } while (0)

static void
_mpn_mod_unred_mul_acc_2_3(nn_ptr c, nn_srcptr G, nn_srcptr H, slong d)
{
    slong idx, j, jlo, jhi;

    for (idx = 0; idx <= 2 * d - 2; idx++)
    {
        ulong r0 = 0, r1 = 0, r2 = 0;
        ulong u1 = 0, u2 = 0, v2 = 0;
        ulong p0, p1, p2;

        jlo = FLINT_MAX(0, idx - d + 1);
        jhi = FLINT_MIN(idx, d - 1);

        for (j = jlo; j <= jhi; j++)
        {
            ulong A0 = G[2 * j], A1 = G[2 * j + 1];
            ulong B0 = H[2 * (idx - j)], B1 = H[2 * (idx - j) + 1];

            umul_ppmm(p2, p1, A1, B0);
            add_ssaaaa(u2, u1, u2, u1, p2, p1);
            p2 = A1 * B1;
            umul_ppmm(p1, p0, A0, B0);
            add_sssaaaaaa(r2, r1, r0, r2, r1, r0, p2, p1, p0);
            umul_ppmm(p2, p1, A0, B1);
            add_ssaaaa(v2, u1, v2, u1, p2, p1);
        }

        u2 = u2 + v2;
        add_ssaaaa(r2, r1, r2, r1, u2, u1);

        c[3 * idx + 0] = r0;
        c[3 * idx + 1] = r1;
        c[3 * idx + 2] = r2;
    }
}

static void
_mpn_mod_unred_sqr_acc_2_3(nn_ptr c, nn_srcptr G, slong d)
{
    slong idx, j, jlo;

    for (idx = 0; idx <= 2 * d - 2; idx++)
    {
        ulong r0 = 0, r1 = 0, r2 = 0;
        ulong u1 = 0, u2 = 0, v2 = 0;
        ulong p0, p1, p2;

        jlo = FLINT_MAX(0, idx - d + 1);

        /* cross terms G_j G_{idx - j} with j < idx - j, counted twice */
        for (j = jlo; 2 * j < idx; j++)
        {
            ulong A0 = G[2 * j], A1 = G[2 * j + 1];
            ulong B0 = G[2 * (idx - j)], B1 = G[2 * (idx - j) + 1];

            umul_ppmm(p2, p1, A1, B0);
            add_ssaaaa(u2, u1, u2, u1, p2, p1);
            p2 = A1 * B1;
            umul_ppmm(p1, p0, A0, B0);
            add_sssaaaaaa(r2, r1, r0, r2, r1, r0, p2, p1, p0);
            umul_ppmm(p2, p1, A0, B1);
            add_ssaaaa(v2, u1, v2, u1, p2, p1);
        }

        u2 = u2 + v2;
        add_ssaaaa(r2, r1, r2, r1, u2, u1);

        /* double */
        add_sssaaaaaa(r2, r1, r0, r2, r1, r0, r2, r1, r0);

        /* square term */
        if (idx % 2 == 0)
        {
            _MPN_MOD_SQR_3P2X2(p2, p1, p0, G[idx + 1], G[idx]);
            add_sssaaaaaa(r2, r1, r0, r2, r1, r0, p2, p1, p0);
        }

        c[3 * idx + 0] = r0;
        c[3 * idx + 1] = r1;
        c[3 * idx + 2] = r2;
    }
}

static void
_mpn_mod_unred_mul_acc_2_4(nn_ptr c, nn_srcptr G, nn_srcptr H, slong d)
{
    slong idx, j, jlo, jhi;

    for (idx = 0; idx <= 2 * d - 2; idx++)
    {
        ulong r0 = 0, r1 = 0, r2 = 0, r3 = 0;
        ulong p0, p1, p2, p3;

        jlo = FLINT_MAX(0, idx - d + 1);
        jhi = FLINT_MIN(idx, d - 1);

        for (j = jlo; j <= jhi; j++)
        {
            FLINT_MPN_MUL_2X2(p3, p2, p1, p0,
                              G[2 * j + 1], G[2 * j],
                              H[2 * (idx - j) + 1], H[2 * (idx - j)]);
            add_ssssaaaaaaaa(r3, r2, r1, r0, r3, r2, r1, r0, p3, p2, p1, p0);
        }

        c[4 * idx + 0] = r0;
        c[4 * idx + 1] = r1;
        c[4 * idx + 2] = r2;
        c[4 * idx + 3] = r3;
    }
}

static void
_mpn_mod_unred_sqr_acc_2_4(nn_ptr c, nn_srcptr G, slong d)
{
    slong idx, j, jlo;

    for (idx = 0; idx <= 2 * d - 2; idx++)
    {
        ulong r0 = 0, r1 = 0, r2 = 0, r3 = 0;
        ulong p0, p1, p2, p3;

        jlo = FLINT_MAX(0, idx - d + 1);

        for (j = jlo; 2 * j < idx; j++)
        {
            FLINT_MPN_MUL_2X2(p3, p2, p1, p0,
                              G[2 * j + 1], G[2 * j],
                              G[2 * (idx - j) + 1], G[2 * (idx - j)]);
            add_ssssaaaaaaaa(r3, r2, r1, r0, r3, r2, r1, r0, p3, p2, p1, p0);
        }

        add_ssssaaaaaaaa(r3, r2, r1, r0, r3, r2, r1, r0, r3, r2, r1, r0);

        if (idx % 2 == 0)
        {
            FLINT_MPN_SQR_2X2(p3, p2, p1, p0, G[idx + 1], G[idx]);
            add_ssssaaaaaaaa(r3, r2, r1, r0, r3, r2, r1, r0, p3, p2, p1, p0);
        }

        c[4 * idx + 0] = r0;
        c[4 * idx + 1] = r1;
        c[4 * idx + 2] = r2;
        c[4 * idx + 3] = r3;
    }
}

static void
_mpn_mod_unred_mul_acc_2_5(nn_ptr c, nn_srcptr G, nn_srcptr H, slong d)
{
    slong idx, j, jlo, jhi;

    for (idx = 0; idx <= 2 * d - 2; idx++)
    {
        ulong r0 = 0, r1 = 0, r2 = 0, r3 = 0, r4 = 0;
        ulong p0, p1, p2, p3;

        jlo = FLINT_MAX(0, idx - d + 1);
        jhi = FLINT_MIN(idx, d - 1);

        for (j = jlo; j <= jhi; j++)
        {
            FLINT_MPN_MUL_2X2(p3, p2, p1, p0,
                              G[2 * j + 1], G[2 * j],
                              H[2 * (idx - j) + 1], H[2 * (idx - j)]);
            add_sssssaaaaaaaaaa(r4, r3, r2, r1, r0,
                                r4, r3, r2, r1, r0,
                                0, p3, p2, p1, p0);
        }

        c[5 * idx + 0] = r0;
        c[5 * idx + 1] = r1;
        c[5 * idx + 2] = r2;
        c[5 * idx + 3] = r3;
        c[5 * idx + 4] = r4;
    }
}

static void
_mpn_mod_unred_sqr_acc_2_5(nn_ptr c, nn_srcptr G, slong d)
{
    slong idx, j, jlo;

    for (idx = 0; idx <= 2 * d - 2; idx++)
    {
        ulong r0 = 0, r1 = 0, r2 = 0, r3 = 0, r4 = 0;
        ulong p0, p1, p2, p3;

        jlo = FLINT_MAX(0, idx - d + 1);

        for (j = jlo; 2 * j < idx; j++)
        {
            FLINT_MPN_MUL_2X2(p3, p2, p1, p0,
                              G[2 * j + 1], G[2 * j],
                              G[2 * (idx - j) + 1], G[2 * (idx - j)]);
            add_sssssaaaaaaaaaa(r4, r3, r2, r1, r0,
                                r4, r3, r2, r1, r0,
                                0, p3, p2, p1, p0);
        }

        add_sssssaaaaaaaaaa(r4, r3, r2, r1, r0,
                            r4, r3, r2, r1, r0,
                            r4, r3, r2, r1, r0);

        if (idx % 2 == 0)
        {
            FLINT_MPN_SQR_2X2(p3, p2, p1, p0, G[idx + 1], G[idx]);
            add_sssssaaaaaaaaaa(r4, r3, r2, r1, r0,
                                r4, r3, r2, r1, r0,
                                0, p3, p2, p1, p0);
        }

        c[5 * idx + 0] = r0;
        c[5 * idx + 1] = r1;
        c[5 * idx + 2] = r2;
        c[5 * idx + 3] = r3;
        c[5 * idx + 4] = r4;
    }
}

static void
_mpn_mod_unred_mul_acc_3(nn_ptr c, nn_srcptr G, nn_srcptr H, slong d, slong s)
{
    slong idx, j, jlo, jhi;
    ulong t[6];

    for (idx = 0; idx <= 2 * d - 2; idx++)
    {
        ulong r0 = 0, r1 = 0, r2 = 0, r3 = 0, r4 = 0, r5 = 0, r6 = 0;

        jlo = FLINT_MAX(0, idx - d + 1);
        jhi = FLINT_MIN(idx, d - 1);

        if (s == 6)
        {
            for (j = jlo; j <= jhi; j++)
            {
                flint_mpn_mul_n(t, G + 3 * j, H + 3 * (idx - j), 3);
                add_ssssssaaaaaaaaaaaa(r5, r4, r3, r2, r1, r0,
                                       r5, r4, r3, r2, r1, r0,
                                       t[5], t[4], t[3], t[2], t[1], t[0]);
            }
        }
        else
        {
            for (j = jlo; j <= jhi; j++)
            {
                flint_mpn_mul_n(t, G + 3 * j, H + 3 * (idx - j), 3);
                add_sssssssaaaaaaaaaaaaaa(r6, r5, r4, r3, r2, r1, r0,
                                          r6, r5, r4, r3, r2, r1, r0,
                            0, t[5], t[4], t[3], t[2], t[1], t[0]);
            }
        }

        c[s * idx + 0] = r0;
        c[s * idx + 1] = r1;
        c[s * idx + 2] = r2;
        c[s * idx + 3] = r3;
        c[s * idx + 4] = r4;
        c[s * idx + 5] = r5;
        if (s == 7)
            c[s * idx + 6] = r6;
    }
}

static void
_mpn_mod_unred_sqr_acc_3(nn_ptr c, nn_srcptr G, slong d, slong s)
{
    slong idx, j, jlo;
    ulong t[6];

    for (idx = 0; idx <= 2 * d - 2; idx++)
    {
        ulong r0 = 0, r1 = 0, r2 = 0, r3 = 0, r4 = 0, r5 = 0, r6 = 0;

        jlo = FLINT_MAX(0, idx - d + 1);

        if (s == 6)
        {
            for (j = jlo; 2 * j < idx; j++)
            {
                flint_mpn_mul_n(t, G + 3 * j, G + 3 * (idx - j), 3);
                add_ssssssaaaaaaaaaaaa(r5, r4, r3, r2, r1, r0,
                                       r5, r4, r3, r2, r1, r0,
                                       t[5], t[4], t[3], t[2], t[1], t[0]);
            }

            add_ssssssaaaaaaaaaaaa(r5, r4, r3, r2, r1, r0,
                                   r5, r4, r3, r2, r1, r0,
                                   r5, r4, r3, r2, r1, r0);

            if (idx % 2 == 0)
            {
                flint_mpn_sqr(t, G + 3 * (idx / 2), 3);
                add_ssssssaaaaaaaaaaaa(r5, r4, r3, r2, r1, r0,
                                       r5, r4, r3, r2, r1, r0,
                                       t[5], t[4], t[3], t[2], t[1], t[0]);
            }
        }
        else
        {
            for (j = jlo; 2 * j < idx; j++)
            {
                flint_mpn_mul_n(t, G + 3 * j, G + 3 * (idx - j), 3);
                add_sssssssaaaaaaaaaaaaaa(r6, r5, r4, r3, r2, r1, r0,
                                          r6, r5, r4, r3, r2, r1, r0,
                            0, t[5], t[4], t[3], t[2], t[1], t[0]);
            }

            add_sssssssaaaaaaaaaaaaaa(r6, r5, r4, r3, r2, r1, r0,
                                      r6, r5, r4, r3, r2, r1, r0,
                                      r6, r5, r4, r3, r2, r1, r0);

            if (idx % 2 == 0)
            {
                flint_mpn_sqr(t, G + 3 * (idx / 2), 3);
                add_sssssssaaaaaaaaaaaaaa(r6, r5, r4, r3, r2, r1, r0,
                                          r6, r5, r4, r3, r2, r1, r0,
                            0, t[5], t[4], t[3], t[2], t[1], t[0]);
            }
        }

        c[s * idx + 0] = r0;
        c[s * idx + 1] = r1;
        c[s * idx + 2] = r2;
        c[s * idx + 3] = r3;
        c[s * idx + 4] = r4;
        c[s * idx + 5] = r5;
        if (s == 7)
            c[s * idx + 6] = r6;
    }
}

/* Generic kernels reusing the mullow dot products. */
static void
_mpn_mod_unred_mul_acc(nn_ptr c, slong s, nn_srcptr G, slong len1,
                                    nn_srcptr H, slong len2, slong n)
{
    slong idx, jlo, jhi, dlen;

    for (idx = 0; idx <= len1 + len2 - 2; idx++)
    {
        nn_srcptr a, b;

        jlo = FLINT_MAX(0, idx - len2 + 1);
        jhi = FLINT_MIN(idx, len1 - 1);
        dlen = jhi - jlo + 1;
        a = G + jlo * n;
        b = H + (idx - jhi) * n;

        if (n == 3 && s == 5)
            _mpn_dot_rev_3x3_5(c + idx * s, a, b, dlen);
        else if (s == 2 * n - 1)
            _mpn_dot_rev_nxn_2nm1(c + idx * s, a, b, dlen, n);
        else if (s == 2 * n)
            _mpn_dot_rev_nxn_2n(c + idx * s, a, b, dlen, n);
        else
            _mpn_dot_rev_nxn_2np1(c + idx * s, a, b, dlen, n);
    }
}

static void
_mpn_mod_unred_sqr_acc(nn_ptr c, slong s, nn_srcptr G, slong len, slong n)
{
    slong idx, start, stop, dlen;
    ulong t[2 * MPN_MOD_MAX_LIMBS];
    ulong cy;

    for (idx = 0; idx <= 2 * len - 2; idx++)
    {
        nn_ptr ci = c + idx * s;

        /* cross terms G_j G_{idx - j} with j < idx - j, counted twice */
        start = FLINT_MAX(0, idx - len + 1);
        stop = FLINT_MIN(len - 1, (idx + 1) / 2 - 1);
        dlen = stop - start + 1;

        if (dlen <= 0)
        {
            flint_mpn_zero(ci, s);
        }
        else
        {
            if (n == 3 && s == 5)
                _mpn_dot_rev_3x3_5(ci, G + start * n, G + (idx - stop) * n, dlen);
            else if (s == 2 * n - 1)
                _mpn_dot_rev_nxn_2nm1(ci, G + start * n, G + (idx - stop) * n, dlen, n);
            else if (s == 2 * n)
                _mpn_dot_rev_nxn_2n(ci, G + start * n, G + (idx - stop) * n, dlen, n);
            else
                _mpn_dot_rev_nxn_2np1(ci, G + start * n, G + (idx - stop) * n, dlen, n);

            cy = mpn_lshift(ci, ci, s, 1);
            FLINT_ASSERT(cy == 0);
            (void) cy;
        }

        /* square term */
        if (idx % 2 == 0)
        {
            flint_mpn_sqr(t, G + (idx / 2) * n, n);

            if (s > 2 * n)
                ci[2 * n] += mpn_add_n(ci, ci, t, 2 * n);
            else
                mpn_add_n(ci, ci, t, s);
        }
    }
}

void
_mpn_mod_poly_mul_unreduced(nn_ptr res, slong slimbs, nn_srcptr poly1,
        slong len1, nn_srcptr poly2, slong len2, gr_ctx_t ctx)
{
    slong n = MPN_MOD_CTX_NLIMBS(ctx);

    FLINT_ASSERT(len1 >= 1 && len2 >= 1);
    FLINT_ASSERT(slimbs >= _mpn_mod_poly_mul_unreduced_slimbs(FLINT_MIN(len1, len2), ctx));
    FLINT_ASSERT(slimbs >= 2 * n - 1 && slimbs <= 2 * n + 1);

    if (n == 2 && len1 == len2)
    {
        if (slimbs == 3)
            _mpn_mod_unred_mul_acc_2_3(res, poly1, poly2, len1);
        else if (slimbs == 4)
            _mpn_mod_unred_mul_acc_2_4(res, poly1, poly2, len1);
        else
            _mpn_mod_unred_mul_acc_2_5(res, poly1, poly2, len1);
    }
    else if (n == 3 && slimbs >= 6 && len1 == len2)
        _mpn_mod_unred_mul_acc_3(res, poly1, poly2, len1, slimbs);
    else
        _mpn_mod_unred_mul_acc(res, slimbs, poly1, len1, poly2, len2, n);
}

void
_mpn_mod_poly_sqr_unreduced(nn_ptr res, slong slimbs, nn_srcptr poly,
                                                slong len, gr_ctx_t ctx)
{
    slong n = MPN_MOD_CTX_NLIMBS(ctx);

    FLINT_ASSERT(len >= 1);
    FLINT_ASSERT(slimbs >= _mpn_mod_poly_mul_unreduced_slimbs(len, ctx));
    FLINT_ASSERT(slimbs >= 2 * n - 1 && slimbs <= 2 * n + 1);

    if (n == 2)
    {
        if (slimbs == 3)
            _mpn_mod_unred_sqr_acc_2_3(res, poly, len);
        else if (slimbs == 4)
            _mpn_mod_unred_sqr_acc_2_4(res, poly, len);
        else
            _mpn_mod_unred_sqr_acc_2_5(res, poly, len);
    }
    else if (n == 3 && slimbs >= 6)
        _mpn_mod_unred_sqr_acc_3(res, poly, len, slimbs);
    else
        _mpn_mod_unred_sqr_acc(res, slimbs, poly, len, n);
}
