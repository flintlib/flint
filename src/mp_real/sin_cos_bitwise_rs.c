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
#include "mp_real.h"
#include "thread_support.h"
#include "impl.h"

/* The cached table of the angles A_i = atan(2^-i) (thread-local like the
   logarithms of exp_bitwise_rs.c; the reductions of atan, sin, cos and
   tan read it), and sin and cos on [0, 1) by bitwise reduction.

   The latter is the tangent half-angle reconstruction of
   tan_bitwise_rs.c (_mp_real_tan_halfangle): x/2 is reduced with the
   angles A_i, which satisfy the windowed reductions' requirement
   A_i < 2^-i (the doubled angles 2 atan(2^-i) of the rotation identity
   would not), the used factors build the rotation W, and the outputs
   come from the half-angle formulas, one reciprocal for both.  On 64-bit
   machines r = 0 at n <= 10 selects the per-size kernels instead (small.h
   at one to three limbs, the dense tables from four), and the reduction
   runs with signed digits on vectors (bitwise_vec.c) from
   MP_REAL_TRIG_BITWISE_VEC_MIN limbs.

   Error accounting (output ulp 2^(-FLINT_BITS n)): reduction table and
   creep < 2 (r + 2) ulp on the angle; the residual's series <= 15 ulp;
   rotation truncations < 2 ulp of angle per used factor; the products
   and the reciprocal of the reconstruction are constants.  Total below
   the returned bound *err = 6 r + 128. */

/* the shared cache of the table (thread_support.h; see the logarithm
   table of exp_bitwise_rs.c) */
static flint_shared_cache_struct _mp_real_atans_cache =
    FLINT_SHARED_CACHE_INIT(_mp_real_tab_free);

#define ATANS_A(i) (tab + (i) * nc)


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

   and beyond that _mp_real_atan_2mexp_ui_bs (tab_bsplit.c) splits the
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
   side.

   The entries i_from..rc of (tab, entries of nc limbs), rc < 64 nc. */
static void
_mp_real_atans_fill(nn_ptr tab, slong nc, slong i_from, slong rc)
{
    slong i, k, prec, iter_start, i_lo;
    nn_ptr t, scratch, num;

    prec = FLINT_BITS * (nc - 1);

    if (prec <= 65536)
        iter_start = FLINT_MIN(n_sqrt(prec), rc + 1);
    else
        iter_start = rc + 1;

    /* tier 1: binary splitting; entry 0 = atan(1) = pi/4 is unused
       by the reductions (they start at i = 1) but tabulated anyway */
    if (i_from < iter_start)
    {
        arb_ptr th;
        arb_t x;
        fmpz_t p, q;
        slong wp = FLINT_BITS * nc + 30;

        arb_init(x);
        fmpz_init(p);
        fmpz_init(q);
        th = _arb_vec_init(4);

        if (i_from <= 3)
            arb_atan_gauss_primes_vec_bsplit(th, 4, wp);

        for (i = i_from; i < iter_start && i <= rc; i++)
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
                    _mp_real_atan_2mexp_ui_bs(ATANS_A(i), (ulong) i, nc);
                    continue;
            }

            if (!_mp_real_tab_store_floor(ATANS_A(i), x, nc, wp))
                _mp_real_tab_entry_exact(ATANS_A(i), 1, (ulong) i, nc);
        }

        _arb_vec_clear(th, 4);
        arb_clear(x);
        fmpz_clear(p);
        fmpz_clear(q);
    }

    /* tier 2: fixed-point multi-summation, directly in the target
       storage */
    i_lo = FLINT_MAX(iter_start, i_from);
    if (i_lo <= rc)
    {
        slong wp = FLINT_BITS * nc;

        t = flint_malloc((nc + 2) * sizeof(ulong));
        scratch = flint_malloc(nc * sizeof(ulong));
        num = flint_malloc((nc + 1) * sizeof(ulong));

        /* the k = 1 term is the single bit 2^-i */
        for (i = i_lo; i <= rc; i++)
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

        for (k = 3; k * i_lo < wp; k += 2)
        {
            flint_mpn_zero(num, nc);
            num[nc] = 1;                    /* 2^wp */
            mpn_divrem_1(t, 0, num, nc + 1, (ulong) k);
            FLINT_ASSERT(t[nc] == 0);

            for (i = i_lo; i <= rc && i * k < wp; i++)
                _mp_real_tab_addshift(ATANS_A(i), t, nc, i * k,
                    (k & 3) == 3, scratch);
        }

        flint_free(t);
        flint_free(scratch);
        flint_free(num);
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

        for (i = FLINT_MAX(FLINT_MAX(band, 4), i_from); i <= rc; i++)
        {
            nn_ptr acc = ATANS_A(i);
            slong b = wp - i;   /* rc < 64 nc, so b >= 1 */
            flint_mpn_zero(acc, nc);
            flint_mpn_store(acc, b / FLINT_BITS, ~UWORD(0));
            if (b % FLINT_BITS)
                acc[b / FLINT_BITS] =
                    (UWORD(1) << (b % FLINT_BITS)) - 1;
        }

        /* exactness rescan of the fast tiers below the band; see
           the log table's twin loop for the reasoning */
        for (i = FLINT_MAX(4, i_from); i < FLINT_MIN(band, rc + 1); i++)
        {
            ulong g = tab[i * nc];
            if (g + MP_REAL_TAB_GUARD_SLACK < g)
                _mp_real_tab_entry_exact(ATANS_A(i), 1, (ulong) i, nc);
        }
    }
}

/* the shared cache protocol (see exp_bitwise_rs.c) */
static void *
_mp_real_atans_build(const void * old, void * args)
{
    _mp_real_tab_struct * T;
    slong i_from;

    T = _mp_real_tab_new(old, args, &i_from);
    _mp_real_atans_fill(T->tab, T->nc, i_from, T->num - 1);
    return T;
}

/* the version covering nv value limbs and the indices 0..rc */
const _mp_real_tab_struct *
_mp_real_atans_ensure(slong nv, slong rc)
{
    _mp_real_tab_req_struct R;

    R.nc = FLINT_MAX(nv + 1, rc / FLINT_BITS + 2);
    R.num = rc + 1;
    return flint_shared_cache_ensure(&_mp_real_atans_cache,
        _mp_real_tab_sufficient, _mp_real_atans_build, &R);
}

/* Read-only view of the table for code outside the library: the top
   n limbs of entry i of the current version (valid until
   flint_cleanup_master). */
nn_srcptr
_mp_real_atans_entry(slong i, slong n)
{
    const _mp_real_tab_struct * T =
        flint_shared_cache_get(&_mp_real_atans_cache);

    FLINT_ASSERT(T != NULL && i >= 0 && i < T->num);
    FLINT_ASSERT(n >= 1 && n <= T->nc - 1);
    return T->tab + i * T->nc + (T->nc - n);
}

/* free the cached table (tuning programs): not while another thread
   may be using it */
void
_mp_real_atans_clear(void)
{
    flint_shared_cache_clear(&_mp_real_atans_cache);
}

/* the largest index the current version covers (0 if none) */
slong
_mp_real_atans_max_index(void)
{
    const _mp_real_tab_struct * T =
        flint_shared_cache_get(&_mp_real_atans_cache);

    return (T == NULL) ? 0 : T->num - 1;
}

FLINT_FORCE_INLINE void
_mp_real_sin_cos_bitwise_rs_noerr(nn_ptr ysin, nn_ptr ycos, nn_srcptr x,
    slong n, int r)
{
    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(r == 0 || r >= 32);
    /* n = 1 is unsupported on 32-bit limbs: the clamp in the
       half-angle path would drop r under the r >= 32 contract of the
       residual series */
    FLINT_ASSERT(FLINT_BITS == 64 || n >= 2);

    _mp_real_tan_halfangle(ysin, ycos, NULL, x, n, r);
}

void
_mp_real_sin_cos_bitwise_rs(nn_ptr ysin, nn_ptr ycos, ulong * err, nn_srcptr x,
    slong n, int r)
{
#if FLINT_BITS == 64
    /* r = 0 at the sizes of the per-size kernels (small.h, the dense
       tables): those serve the call (any size they decline takes the
       bitwise path below) */
    if (r == 0 && n <= MP_REAL_TRIG_OPT_R_MAX)
    {
        ulong e, ts[MP_REAL_TRIG_OPT_R_MAX + 1], tc[MP_REAL_TRIG_OPT_R_MAX + 1];

        if (_mp_real_sin_cos_opt((ysin != NULL) ? ysin : ts,
                (ycos != NULL) ? ycos : tc, &e, x, n))
        {
            if (err != NULL)
                *err = e;
            return;
        }
    }
#endif

    _mp_real_sin_cos_bitwise_rs_noerr(ysin, ycos, x, n, r);
    if (err != NULL)
        *err = 6 * (r ? (ulong) r : (ulong) _mp_real_trig_bitwise_rs_default_r(n)) + 128;
}
