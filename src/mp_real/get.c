/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <string.h>
#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "arb.h"
#include "mp_real.h"
#include "impl.h"

/* Conversions: fixed-point export, verified floors, arb, printing. */

static double
_mp_real_get_fixed_d(nn_ptr y, slong wn, const mp_real_t x)
{
    double e;
    slong sh;

    /* radius in output ulps: err B^(exp - size + wn), saturating.
       Both saturation directions must OVERestimate: far above, the
       bound is useless anyway (HUGE_VAL); far below, err < B
       at k <= -k0 limbs is at most 2^(FLINT_BITS (1 - k0)), so the
       clamp constant must exceed that, while the ldexp in the live
       window (|FLINT_BITS k| <= 900) cannot underflow */
    if (x->err == 0)
        e = 0.0;
    else
    {
        slong k = x->exp - x->size + wn;
        if (FLINT_BITS * k > 900)
            e = HUGE_VAL;
        else if (FLINT_BITS * k < -900)
            e = 0x1p-700;   /* true bound < 2^(64 - 900) */
        else
            e = ldexp((double) x->err, (int) (FLINT_BITS * k)) * MP_REAL_EPS;
    }

    if (x->size == 0)
    {
        flint_mpn_zero(y, wn);
        return e + 1.0;
    }

    if (x->negative)
    {
        /* nonpositive point value in [0,1) ball: clamp; the true
           value is within (radius - value) <= radius of 0 */
        flint_mpn_zero(y, wn);
        return e;
    }

    FLINT_ASSERT(x->exp <= 0);  /* value < 1 */

    sh = wn + x->exp - x->size;
    flint_mpn_zero(y, wn);
    if (sh >= 0)
    {
        FLINT_ASSERT(sh + x->size <= wn);
        flint_mpn_copyi(y + sh, x->d, x->size);
        return e;
    }
    /* drop the low -sh limbs: one-sided < 1 output ulp */
    if (-sh < x->size)
        flint_mpn_copyi(y, x->d - sh, x->size + sh);
    return e + 1.0;
}

void
_mp_real_get_fixed(nn_ptr y, ulong * err, const mp_real_t x, slong n)
{
    double e = _mp_real_get_fixed_d(y, n, x);

    if (err != NULL)
    {
        /* round up to whole output ulps, saturating */
        if (!(e < 0x1p63))
            *err = UWORD_MAX;
        else
        {
            ulong u = (ulong) e;
            *err = u + ((double) u < e);
        }
    }
}

/* If the ball x, known to satisfy 0 <= x < 1, determines
   floor(x B^n) uniquely, write it to (y, n) and return 1; otherwise
   return 0 (caller retries at higher precision).  Everything is
   limb-aligned: the mantissa splits at the output grid into the kept
   top part M and a t-limb tail L, and the radius -- overestimated by
   the pure power of two 2^s, s = the bit length of err -- must fit
   strictly inside the tail on both
   sides: 2^s <= L (the true value cannot borrow below the grid line)
   and L < B^t - 2^s (nor carry above the next), both checked by bit
   scans without materializing L + 2^s. */
int
_mp_real_get_fixed_floor(nn_ptr y, slong n, const mp_real_t x)
{
    slong sh, t, s, i, hi;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(x->exp <= 0);
    FLINT_ASSERT(x->size == 0 || !x->negative);

    flint_mpn_zero(y, n);

    if (x->size == 0)
    {
        /* 0 +/- err: the floor is 0 iff the radius stays strictly
           below B^-n: err B^exp < B^-n */
        if (x->err == 0)
            return 1;
        return (slong) FLINT_BIT_COUNT(x->err) + FLINT_BITS * (x->exp + n) <= 0;
    }

    sh = x->exp - x->size + n;  /* x B^n = (d, size) B^sh */
    t = -sh;

    if (t <= 0)
    {
        /* mantissa entirely on or above the grid: exact values
           only (any radius straddles the grid line the value sits
           on) */
        if (x->err != 0)
            return 0;
        FLINT_ASSERT(x->size + sh <= n);
        flint_mpn_copyi(y + sh, x->d, x->size);
        return 1;
    }

    if (t < x->size)
        flint_mpn_copyi(y, x->d + t, x->size - t);
    /* else the whole value sits below the grid: floor 0, and the
       t-limb tail is the mantissa padded above with zero limbs */

    if (x->err == 0)
        return 1;

    s = (slong) FLINT_BIT_COUNT(x->err);

    if (s >= FLINT_BITS * t - 1)
        return 0;

    /* the tail L, read as t limbs with zeros above the mantissa */
#define MP_REAL_TAIL_LIMB(I_) ((I_) < x->size ? x->d[I_] : UWORD(0))

    /* low side: L >= 2^s, i.e. some bit at position >= s is set */
    hi = -1;
    for (i = t - 1; i >= 0; i--)
        if (MP_REAL_TAIL_LIMB(i) != 0)
        {
            hi = i;
            break;
        }
    if (hi < 0 || FLINT_BITS * hi
            + (slong) FLINT_BIT_COUNT(MP_REAL_TAIL_LIMB(hi)) <= s)
        return 0;

    /* high side: L < B^t - 2^s, i.e. some bit at position >= s is
       clear */
    {
        slong ls = s / FLINT_BITS;
        int sb = (int) (s % FLINT_BITS);

        for (i = t - 1; i > ls; i--)
            if (MP_REAL_TAIL_LIMB(i) != UWORD_MAX)
                return 1;
        return (MP_REAL_TAIL_LIMB(ls) >> sb) != (UWORD_MAX >> sb);
    }
#undef MP_REAL_TAIL_LIMB
}

void
mp_real_get_arb(arb_t res, const mp_real_t x)
{
    if (x->size == 0)
        arf_zero(arb_midref(res));
    else
    {
        /* the mantissa is exact at size limbs; value d B^(exp - size)
           has arf exponent FLINT_BITS exp + (bits shifted out) */
        slong fix;
        _arf_set_round_mpn(arb_midref(res), &fix, x->d, x->size,
            x->negative, FLINT_BITS * x->size, ARF_RND_DOWN);
        fmpz_set_si(ARF_EXPREF(arb_midref(res)),
            FLINT_BITS * x->exp + fix);
    }

    if (x->err == 0)
        mag_zero(arb_radref(res));
    else
    {
        mag_set_ui(arb_radref(res), x->err);
        mag_mul_2exp_si(arb_radref(res), arb_radref(res),
            FLINT_BITS * (x->exp - x->size));
    }
}

void
mp_real_print(const mp_real_t x)
{
    arb_t t;
    arb_init(t);
    mp_real_get_arb(t, x);
    flint_printf("[size %wd, exp %wd, err %wu] = ",
        x->size, x->exp, x->err);
    arb_printd(t, 20);
    flint_printf("\n");
    arb_clear(t);
}
