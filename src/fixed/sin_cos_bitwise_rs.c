/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpn_extras.h"
#include "arb.h"
#include "fixed.h"

/* sin and cos on [0, 1) by the rotation analog of the bitwise exp
   algorithm, following the identity of [Joh2022]: with
   alpha = 2 atan(b/a),
   exp(i x) = exp(i (x - k alpha)) ((a + b i)/(a - b i))^k, whose
   rotation unit is a ratio of conjugates and hence EXACTLY
   unimodular -- no normalizing square root arises.  With a = 2^i,
   b = 1 the rotation angles are 2 atan(2^-i).

   The TABLE, however, stores the UNDOUBLED angles
   A_i = atan(2^-i) (thread-locally cached like the exp
   logarithms), and phase 1 reduces x/2 rather than x.  This is not
   merely an economy of one table shared with
   fixed_atan_bitwise_rs: the windowed single-limb decision model of
   _fixed_bitwise_reduce requires tab_i < 2^-i, which A_i satisfies
   (as log(1 + 2^-i) does) but 2 atan(2^-i) ~ 2^(1-i) does NOT --
   with doubled entries the residual can keep one bit above the
   window limb after a boundary step, desynchronizing the model.
   Concavity of atan gives A_{i-1} < 2 A_i, so each index is used at
   most once, and sum_{i>=1} A_i ~ 0.8980 > 1/2 covers x/2.

   Writing x/2 = sum_{used} A_i + t', the used factors account for
   the rotation by sum 2 A_i and the residual angle is 2 t'.  The
   reduction therefore runs to index r + 1, so that
   2 t' < 2 A_{r+1} < 2^-r meets the series contract.

   Phase 2 evaluates sin and cos of the residual t < 2^-r: for short
   series both at once by fixed_sin_cos_rs (shared powers of t^2);
   once the series gets long, only the odd sine series, with
   cos t = sqrt(1 - sin^2 t) by one squaring and one integer square
   root -- the same tradeoff and threshold shape as exp's sinh path
   (for small n the square root costs more than the extra terms).

   Phase 3 reconstructs: with W = prod (1 + i 2^-i) over the used
   factors (the real 2^i denominators cancel between W and its
   conjugate), each factor is two shifts and two add/subtracts on
   the old components, and

       (cos x, sin x) = components of (c + i s) W^2 / |W|^2.

   Two squarings give Re W^2 = wx^2 - wy^2 and |W|^2 = wx^2 + wy^2
   from the same pair, one product gives Im W^2 = 2 wx wy, a complex
   multiplication by (c + i s) and two divisions by the real |W|^2
   finish.  The self-normalization is structural: the ratio of the
   SAME truncated W against its own conjugate is exactly unimodular
   up to the final roundings, so the reconstruction truncations
   perturb only the angle of the correction (O(r) ulp), never its
   modulus.  |W| < prod sqrt(1 + 4^-i) = 1.16444, so wx in
   (0.72, 1.17) and wy in [0, 0.92) fit n + 1 limbs, and
   |W|^2 in [1, 1.3559).

   Error accounting (output ulp 2^(-FLINT_BITS n)): reduction table
   and creep < 2 (r + 2) ulp on the angle; series <= 15 ulp per
   component (plus squaring and square root floors on the sqrt
   path); reconstruction truncations < 2 ulp of angle per used
   factor; the final squarings, products and division floors are
   constants.  Total below
   FIXED_SIN_COS_BITWISE_RS_MAX_ERR(n, r) = 6 r + 128. */

FLINT_TLS_PREFIX nn_ptr _fixed_atans = NULL;
FLINT_TLS_PREFIX slong _fixed_atans_n = 0;
FLINT_TLS_PREFIX slong _fixed_atans_r = 0;

/* free the cached table; see _fixed_exp_logs_clear */
void
_fixed_atans_clear(void)
{
    flint_free(_fixed_atans);
    _fixed_atans = NULL;
    _fixed_atans_n = 0;
    _fixed_atans_r = 0;
}
static FLINT_TLS_PREFIX int _fixed_atans_cleanup_registered = 0;

static void
_fixed_atans_cleanup(void)
{
    flint_free(_fixed_atans);
    _fixed_atans = NULL;
    _fixed_atans_n = 0;
    _fixed_atans_r = 0;
    _fixed_atans_cleanup_registered = 0;
}

#define ATANS_A(i) (_fixed_atans + (i) * _fixed_atans_n)


/* acc (nc limbs) += or -= t >> s; the bits shifted below position 0
   fall below the guard limb and are dropped.  scratch: nc limbs. */
static void
_fixed_tab_addshift(nn_ptr acc, nn_srcptr t, slong nc, slong s,
    int sub, nn_ptr scratch)
{
    slong q = s / FLINT_BITS, m = nc - q;
    int b = (int) (s % FLINT_BITS);
    nn_srcptr v;

    if (m <= 0)
        return;

    if (b == 0)
        v = t + q;
    else
    {
        mpn_rshift(scratch, t + q, m, b);
        v = scratch;
    }

    if (sub)
    {
        /* subtracted terms round UP (one extra guard ulp stands in
           for the dropped bits): together with the blanket bias
           applied to the leading term this keeps every entry AT OR
           BELOW the exact value, preserving the strict tab_i < 2^-i
           invariant of the reductions */
        mpn_sub(acc, acc, nc, v, m);
        mpn_sub_1(acc, acc, nc, 1);
    }
    else
        mpn_add(acc, acc, nc, v, m);
}

/* The table entries A_i = atan(2^-i), each nc limbs (the value
   scaled by 2^(FLINT_BITS nc), whose bottom limb is a guard), built
   in the same two tiers as the logarithm table in exp_bitwise_rs.c.

   Small i: binary splitting.  atan(2^-i) for i <= 3 falls out of the
   first four Gaussian-prime angles theta = {atan 1, atan 2,
   atan(3/2), atan 4}:

       atan(1)   = theta_0
       atan(1/2) = 2 theta_0 - theta_1        (complement of atan 2)
       atan(1/4) = 2 theta_0 - theta_3        (complement of atan 4)
       atan(1/8) = theta_1 - theta_2          (8 + i = -i(1+2i)(3+2i))

   and beyond that fixed_atan_2mexp_ui_bs (tab_bsplit.c) splits the
   series natively in mpn arithmetic, writing straight into the
   entry.  TODO: the remaining arb dependency of this tier is the
   gauss-prime vector combination for i <= 3.

   Large i: fixed-point multi-summation of

       atan(2^-i) = sum_{j >= 0} (-1)^j 2^(-i(2j+1)) / (2j+1),

   simpler than the logarithm's: every term has an odd denominator,
   so ONE reciprocal t = floor(2^wp / k) per odd k serves every i,
   added when k = 1 mod 4 and subtracted when k = 3 mod 4, and the
   number of divisions is about wp / (2 iter_start) for the whole
   tier.  The floored terms err below one guard ulp each, leaving
   the value limbs equal to the exact truncation of A_i except when
   A_i lies within about 2^-50 of a limb boundary, and then by one
   ulp at most, which the reduction's accounting absorbs on either
   side. */
void
_fixed_atans_ensure(slong nv, slong rc)
{
    slong nc, i, k, prec, iter_start;
    nn_ptr t, scratch, num;

    if (nv + 1 <= _fixed_atans_n && rc <= _fixed_atans_r)
        return;

    nc = FLINT_MAX(nv + 1, _fixed_atans_n);
    rc = FLINT_MAX(rc, _fixed_atans_r);

    flint_free(_fixed_atans);
    _fixed_atans = flint_malloc((rc + 1) * nc * sizeof(ulong));
    _fixed_atans_n = nc;
    _fixed_atans_r = rc;

    prec = FLINT_BITS * (nc - 1);

    if (prec <= 65536)
        iter_start = FLINT_MIN(n_sqrt(prec), rc + 1);
    else
        iter_start = rc + 1;

    /* tier 1: binary splitting; entry 0 = atan(1) = pi/4 is unused
       by the reductions (they start at i = 1) but tabulated anyway */
    {
        arb_ptr th;
        arb_t x;
        fmpz_t p, q;
        slong wp = FLINT_BITS * nc + 30;

        arb_init(x);
        fmpz_init(p);
        fmpz_init(q);
        th = _arb_vec_init(4);

        arb_atan_gauss_primes_vec_bsplit(th, 4, wp);

        for (i = 0; i < iter_start && i <= rc; i++)
        {
            switch (i)
            {
                case 0:
                    arb_set(x, th + 0);
                    break;
                case 1:
                    arb_mul_2exp_si(x, th + 0, 1);
                    arb_sub(x, x, th + 1, wp);
                    break;
                case 2:
                    arb_mul_2exp_si(x, th + 0, 1);
                    arb_sub(x, x, th + 3, wp);
                    break;
                case 3:
                    arb_sub(x, th + 1, th + 2, wp);
                    break;
                default:
                    /* native mpn binary splitting, straight into the
                       entry (see tab_bsplit.c) */
                    fixed_atan_2mexp_ui_bs(ATANS_A(i), (ulong) i, nc);
                    continue;
            }

            if (!_fixed_tab_store_floor(ATANS_A(i), x, nc, wp))
                _fixed_tab_entry_exact(ATANS_A(i), 1, (ulong) i, nc);
        }

        _arb_vec_clear(th, 4);
        arb_clear(x);
        fmpz_clear(p);
        fmpz_clear(q);
    }

    /* tier 2: fixed-point multi-summation, directly in the target
       storage */
    if (iter_start <= rc)
    {
        slong wp = FLINT_BITS * nc;

        t = flint_malloc((nc + 2) * sizeof(ulong));
        scratch = flint_malloc(nc * sizeof(ulong));
        num = flint_malloc((nc + 1) * sizeof(ulong));

        /* the k = 1 term is the single bit 2^-i */
        for (i = iter_start; i <= rc; i++)
        {
            nn_ptr acc = ATANS_A(i);

            flint_mpn_zero(acc, nc);
            acc[(wp - i) / FLINT_BITS] =
                UWORD(1) << ((wp - i) % FLINT_BITS);
            /* blanket bias: the series terms too small to represent
               (below the guard limb) sum to less than two guard
               ulps; over-subtracting them keeps the entry one-sided */
            mpn_sub_1(acc, acc, nc, 2);
        }

        for (k = 3; k * iter_start < wp; k += 2)
        {
            flint_mpn_zero(num, nc);
            num[nc] = 1;                    /* 2^wp */
            mpn_divrem_1(t, 0, num, nc + 1, (ulong) k);
            FLINT_ASSERT(t[nc] == 0);

            for (i = iter_start; i <= rc && i * k < wp; i++)
                _fixed_tab_addshift(ATANS_A(i), t, nc, i * k,
                    (k & 3) == 3, scratch);
        }

        flint_free(t);
        flint_free(scratch);
        flint_free(num);
    }

    if (!_fixed_atans_cleanup_registered)
    {
        flint_register_cleanup_function(_fixed_atans_cleanup);
        _fixed_atans_cleanup_registered = 1;
    }
    /* Direct tail band: for 3i >= 64 nc every term of
       atan(2^-i) = 2^-i - 2^(-3i)/3 + ... beyond the first falls
       below the entry, and the omitted mass r satisfies
       0 < r 2^(64 nc) <= 1/3, so the exact floor is the pure bit
       pattern 2^(64 nc - i) - 1.  Writing these directly (they are
       exactly the entries whose guard limb legitimately sits at
       0xfff..., which would otherwise trip the rescan below) keeps
       the arb fallback for genuinely borderline entries only. */
    {
        slong wp = FLINT_BITS * nc, band = (wp + 2) / 3;

        for (i = FLINT_MAX(band, 4); i <= rc; i++)
        {
            nn_ptr acc = ATANS_A(i);
            slong b = wp - i;   /* rc < 64 nv by the r convention,
                                   so b >= 1 */
            flint_mpn_zero(acc, nc);
            flint_mpn_store(acc, b / FLINT_BITS, ~UWORD(0));
            if (b % FLINT_BITS)
                acc[b / FLINT_BITS] =
                    (UWORD(1) << (b % FLINT_BITS)) - 1;
        }

        /* exactness rescan of the fast tiers below the band; see
           the log table's twin loop for the reasoning */
        for (i = 4; i < FLINT_MIN(band, rc + 1); i++)
        {
            ulong g = _fixed_atans[i * nc];
            if (g + FIXED_TAB_GUARD_SLACK < g)
                _fixed_tab_entry_exact(ATANS_A(i), 1, (ulong) i, nc);
        }
    }

}

/* Read-only view of the table for code outside the library (the
   thread-local storage itself is not exported; see fixed.h):
   the top n limbs of entry i, valid until the next ensure call on
   this thread. */
nn_srcptr
_fixed_atans_entry(slong i, slong n)
{
    FLINT_ASSERT(i >= 0 && i <= _fixed_atans_r);
    FLINT_ASSERT(n >= 1 && n <= _fixed_atans_n - 1);
    return _fixed_atans + i * _fixed_atans_n + (_fixed_atans_n - n);
}

/* the largest index the cached table currently covers */
slong
_fixed_atans_max_index(void)
{
    return _fixed_atans_r;
}

void
fixed_sin_cos_bitwise_rs(nn_ptr ysin, nn_ptr ycos, nn_srcptr x,
    slong n, int r)
{
    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(r == 0 || r >= 32);
    /* n = 1 is unsupported on 32-bit limbs: the clamp in the
       half-angle path would drop r under the r >= 32 contract of the
       residual series */
    FLINT_ASSERT(FLINT_BITS == 64 || n >= 2);

    /* Everything is done by the tangent half-angle reconstruction (see
       tan_bitwise_rs.c), which measured faster at every size tried,
       from n = 1 to 192.  The conjugate-ratio reconstruction that used
       to live here -- with its |W|^2, its two squarings and its two
       wide divisions -- is gone, and with it the doubled-angle
       reduction it needed. */
    _fixed_tan_halfangle(ysin, ycos, NULL, x, n, r);
}
