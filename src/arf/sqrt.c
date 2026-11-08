/*
    Copyright (C) 2014, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "mpn_extras/fp_small.h"
#include "arf.h"

/*
    arf_sqrt through the root limb-array kernels of arf (arf/mpn_div.c) (see arf_div), with the
    one-limb roots in registers.
*/
int
arf_sqrt(arf_ptr z, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
#if FLINT_BITS == 64
    slong xn;
    nn_srcptr xptr;
#endif
    int t, inexact;

    if (arf_is_special(x))
    {
        if (arf_is_neg_inf(x))
            arf_nan(z);
        else
            arf_set(z, x);
        return 0;
    }

    if (ARF_SGNBIT(x))
    {
        arf_nan(z);
        return 0;
    }

    if (ARF_IS_POW2(x) && fmpz_is_odd(ARF_EXPREF(x)))
    {
        arf_set(z, x);
        fmpz_cdiv_q_2exp(ARF_EXPREF(z), ARF_EXPREF(z), 1);
        return 0;
    }

    if (prec == ARF_PREC_EXACT)
        flint_throw(FLINT_ERROR, "arf_sqrt: exact precision not supported\n");

    t = fmpz_is_odd(ARF_EXPREF(x));

#if FLINT_BITS == 64
    ARF_GET_MPN_READONLY(xptr, xn, x);
    if (xn == 1 && prec <= FLINT_BITS && !COEFF_IS_MPZ(ARF_EXP(x)))
    {
        /* the exact floor s of sqrt(a B / 2^t) (one limb, top bit set) with
           its remainder r = a B / 2^t - s^2 <= 2 s, in registers: s is
           rounded to prec bits with [r != 0] as the sticky bit, and for
           prec = FLINT_BITS with the half bit [r > s] (sqrt >= s + 1/2
           iff r >= s + 1/4; a tie is impossible) */
        ulong s, r1, r0, low, half, mask;
        slong ee = (ARF_EXP(x) + t) >> 1;
        int up, tie_even;

        _flint_mpn_isqrt_2_rem(&s, &r1, &r0, xptr[0] >> t, (xptr[0] << (FLINT_BITS - 1)) & (-(ulong) t));
        inexact = ((r1 | r0) != 0);

        /* the bits of s below prec, and whether the value is above the
           midpoint of the rounding interval */
        if (prec == FLINT_BITS)
        {
            mask = 0;
            low = 0;
            up = (r1 != 0 || r0 > s);
            tie_even = 0;
        }
        else
        {
            mask = (UWORD(1) << (FLINT_BITS - prec)) - 1;
            half = UWORD(1) << (FLINT_BITS - prec - 1);
            low = s & mask;
            inexact |= (low != 0);
            up = (low > half) || (low == half && (r1 | r0) != 0);
            tie_even = (low == half && (r1 | r0) == 0 && ((s >> (FLINT_BITS - prec)) & 1));
        }

        switch (_arf_rnd_resolve(rnd, 0))
        {
            case ARF_RND_NEAR: up = up || tie_even; break;
            case ARF_RND_UP: up = inexact; break;
            default: up = 0;
        }

        s = (s & ~mask) + (up ? (mask + 1) : 0);
        if (s == 0)
        {
            s = UWORD(1) << (FLINT_BITS - 1);
            ee++;
        }

        /* z = s B^-1 2^ee */
        ARF_DEMOTE(z);
        ARF_NOPTR_D(z)[0] = s;
        ARF_XSIZE(z) = ARF_MAKE_XSIZE(1, 0);
        fmpz_set_si(ARF_EXPREF(z), ee);
        return inexact;
    }
#endif

    {
        nn_srcptr xd;
        nn_ptr zd;
        slong xn, rn, fix;
        fmpz_t e;
        ARF_KERNEL_DECL

        ARF_GET_MPN_READONLY(xd, xn, x);
        rn = (prec + FLINT_BITS - 1) / FLINT_BITS;
        ARF_KERNEL_ALLOC(zd, rn, z, x, x);

        /* the exponent (e + t) / 2 of the root (e + t is even) from the
           fmpz exponent, the kernel's correction relative to it */
        fmpz_init(e);
        if (!COEFF_IS_MPZ(ARF_EXP(x)))
            *e = (ARF_EXP(x) + t) >> 1;
        else
        {
            fmpz_add_ui(e, ARF_EXPREF(x), t);
            fmpz_tdiv_q_2exp(e, e, 1);
        }

        inexact = _arf_mpn_sqrt(zd, &fix, xd, xn, t, prec, _arf_rnd_resolve(rnd, 0));

        ARF_KERNEL_FINISH(z, zd, rn, 0);
        _fmpz_add_fast(ARF_EXPREF(z), e, fix);
        fmpz_clear(e);
        return inexact;
    }
}

int
arf_sqrt_ui(arf_t z, ulong x, slong prec, arf_rnd_t rnd)
{
    arf_t t;
    arf_init_set_ui(t, x); /* no need to free */
    return arf_sqrt(z, t, prec, rnd);
}

int
arf_sqrt_fmpz(arf_t z, const fmpz_t x, slong prec, arf_rnd_t rnd)
{
    arf_t t;
    int inexact;
    arf_init(t);
    arf_set_fmpz(t, x);
    inexact = arf_sqrt(z, t, prec, rnd);
    arf_clear(t);
    return inexact;
}
