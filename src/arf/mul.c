/*
    Copyright (C) 2014, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "arf.h"

/*
    arf_mul (arf_mul_rnd_down for the round-down mode of ball arithmetic,
    arf_mul_rnd_any for the others) in all rounding modes:

    - products of one and two limbs in registers (the original round-down
      code, and a register rounding for the other modes at prec <= 128);
    - the full product, rounded by _arf_set_round_mpn, when it is not
      much longer than the output or the operands are small
      (_arf_mul_prefer_full: up to about 7 x 7 limbs the full product is
      as fast as a high product, and exact);
    - otherwise the normalised high product (flint_mpn_mulhigh and its
      balanced kernels: rn limbs with a guard limb within a known number
      of units of the true product), rounded from the guard limb when it
      certifies the rounding (_arf_mpn_round_guarded: the bits of the
      guard limb above the deficit bound neither all zero nor all ones
      prove that the exact product rounds the same way and is inexact),
      and the full product otherwise (probability 2^-55 per product for
      random operands, always when the product is exact at prec bits).

    The exact product never needs MPFR.
*/

/*
    Rounds the normalised 2-limb value (hi, lo) (top bit set), with stop
    the limb below lo and sbot a sticky flag for anything further below,
    to prec <= 2 * FLINT_BITS bits in the given strict rounding mode, in
    registers. Sets the mantissa of z and its exponent to ex + ey +
    expfix (+1 on overflow to a power of two). Returns the inexact flag.
*/
FLINT_FORCE_INLINE int
_arf_mul_round_2(arf_ptr z, ulong hi, ulong lo, ulong stop, ulong sbot, int sgnbit,
    slong prec, arf_rnd_t rnd, slong expfix, arf_srcptr x, arf_srcptr y)
{
    slong k = 2 * FLINT_BITS - prec;   /* discarded bits */
    ulong half, rest, kept_lsb, mask, sticky = stop | sbot;
    int inexact, increment;
    nn_ptr zptr;

    if (k == 0)
    {
        /* the discarded part is the limb below lo (and further bits) */
        half = stop >> (FLINT_BITS - 1);
        rest = (stop << 1) | sbot;
        kept_lsb = lo & 1;
    }
    else if (k <= FLINT_BITS)
    {
        mask = (k == FLINT_BITS) ? ~UWORD(0) : ((UWORD(1) << k) - 1);
        half = (lo >> (k - 1)) & 1;
        rest = (lo & (mask >> 1)) | sticky;
        lo &= ~mask;
        kept_lsb = (k == FLINT_BITS) ? (hi & 1) : ((lo >> k) & 1);
    }
    else
    {
        slong kk = k - FLINT_BITS;   /* 1 .. FLINT_BITS - 1 */
        mask = (UWORD(1) << kk) - 1;
        half = (hi >> (kk - 1)) & 1;
        rest = (hi & (mask >> 1)) | lo | sticky;
        lo = 0;
        hi &= ~mask;
        kept_lsb = (hi >> kk) & 1;
    }

    inexact = (half | rest) != 0;

    if (rnd == ARF_RND_DOWN)
        increment = 0;
    else if (rnd == ARF_RND_NEAR)
        increment = half & ((rest != 0) | kept_lsb);
    else
        increment = inexact & (rnd == ARF_RND_UP || (rnd == ARF_RND_FLOOR ? sgnbit : !sgnbit));

    if (increment)
    {
        if (k < FLINT_BITS)
        {
            lo += UWORD(1) << k;
            hi += (lo < (UWORD(1) << k));
        }
        else
        {
            hi += UWORD(1) << (k - FLINT_BITS);
        }

        if (hi == 0)
        {
            /* overflow to a power of two */
            hi = UWORD(1) << (FLINT_BITS - 1);
            lo = 0;
            expfix++;
        }
    }

    ARF_DEMOTE(z);
    _fmpz_add2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), expfix);
    zptr = ARF_NOPTR_D(z);

    if (lo == 0)
    {
        ARF_XSIZE(z) = ARF_MAKE_XSIZE(1, sgnbit);
        zptr[0] = hi;
    }
    else
    {
        ARF_XSIZE(z) = ARF_MAKE_XSIZE(2, sgnbit);
        zptr[0] = lo;
        zptr[1] = hi;
    }

    return inexact;
}

/* the full product rounded by _arf_set_round_mpn (exact when prec is
   ARF_PREC_EXACT or the product fits); inline: the call and the stack
   frame cost 10% at three and four limbs */
FLINT_FORCE_INLINE int
_arf_mul_full(arf_ptr z, arf_srcptr x, arf_srcptr y, slong xn, slong yn,
    int sgnbit, slong prec, arf_rnd_t rnd)
{
    slong zn = xn + yn, expfix;
    nn_srcptr xd, yd;
    nn_ptr tmp;
    int inexact;
    ARF_MUL_TMP_DECL

    ARF_GET_MPN_READONLY(xd, xn, x);
    ARF_GET_MPN_READONLY(yd, yn, y);
    ARF_MUL_TMP_ALLOC(tmp, zn)
    FLINT_MPN_MUL_WITH_SPECIAL_CASES(tmp, xd, xn, yd, yn);
    inexact = _arf_set_round_mpn(z, &expfix, tmp, zn, sgnbit, prec, rnd);
    _fmpz_add2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), expfix);
    ARF_MUL_TMP_FREE(tmp, zn)
    return inexact;
}

/* the high product into the buffer of z, rounded from its guard limb;
   the full product when the guard limb cannot certify the rounding */
static int
_arf_mul_high(arf_ptr z, arf_srcptr x, arf_srcptr y, slong xn, slong yn,
    int sgnbit, slong prec, arf_rnd_t rnd)
{
    slong rn = (prec + FLINT_BITS - 1) / FLINT_BITS, fix;
    nn_srcptr xd, yd;
    nn_ptr zd;
    mp_limb_pair_t hi;
    int e, inexact;
    ARF_KERNEL_DECL

    ARF_GET_MPN_READONLY(xd, xn, x);
    ARF_GET_MPN_READONLY(yd, yn, y);
    ARF_KERNEL_ALLOC(zd, rn, z, x, y);

    hi = _arf_mpn_mulhigh_normalised(zd, xd, xn, yd, yn, rn);
    e = FLINT_BIT_COUNT(FLINT_MPN_MULHIGH_NORMALISED_GUARD_ERR_BOUND(rn));

    if (_arf_mpn_round_guarded(zd, &fix, zd, rn, hi.m1, e, prec, _arf_rnd_resolve(rnd, sgnbit), 0) >= 0)
    {
        inexact = 1;
        fix -= (slong) hi.m2;
    }
    else
    {
        /* inconclusive: the exact product decides */
        slong zn = xn + yn;
        nn_ptr tmp;
        ARF_MUL_TMP_DECL

        ARF_MUL_TMP_ALLOC(tmp, zn)
        FLINT_MPN_MUL_WITH_SPECIAL_CASES(tmp, xd, xn, yd, yn);
        inexact = _arf_mpn_set_round(zd, &fix, tmp, zn, prec, _arf_rnd_resolve(rnd, sgnbit));
        ARF_MUL_TMP_FREE(tmp, zn)
    }

    ARF_KERNEL_FINISH(z, zd, rn, sgnbit);
    _fmpz_add2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), fix);
    return inexact;
}

/* products of one- and two-limb operands in the strict modes other than
   round-down at prec <= 2 FLINT_BITS: the product of 2, 3 or 4 limbs
   rounded in registers */
static int
_arf_mul_small_rnd(arf_ptr z, arf_srcptr x, arf_srcptr y, slong xn, slong yn,
    int sgnbit, slong prec, arf_rnd_t rnd)
{
    /* stop and sbot are the limbs below lo */
    ulong hi, lo, stop, sbot;
    slong fix;

    if (xn == 1)
    {
        umul_ppmm(hi, lo, ARF_NOPTR_D(x)[0], ARF_NOPTR_D(y)[0]);
        stop = sbot = 0;
    }
    else
    {
        ulong zz[4];
        ulong x1 = ARF_NOPTR_D(x)[1], x0 = ARF_NOPTR_D(x)[0];

        if (yn == 2)
        {
            FLINT_MPN_MUL_2X2(zz[3], zz[2], zz[1], zz[0], x1, x0, ARF_NOPTR_D(y)[1], ARF_NOPTR_D(y)[0]);
            hi = zz[3];
            lo = zz[2];
            stop = zz[1];
            sbot = zz[0];
        }
        else
        {
            FLINT_MPN_MUL_2X1(zz[2], zz[1], zz[0], x1, x0, ARF_NOPTR_D(y)[0]);
            hi = zz[2];
            lo = zz[1];
            stop = zz[0];
            sbot = 0;
        }
    }

    /* normalise (branch free); the top bit of the limbs below moves
       into lo, the rest is only needed as a sticky flag */
    fix = !(hi >> (FLINT_BITS - 1));
    hi = (hi << fix) | ((lo >> (FLINT_BITS - 1)) & fix);
    lo = (lo << fix) | ((stop >> (FLINT_BITS - 1)) & fix);
    stop = (stop << fix) | ((sbot >> (FLINT_BITS - 1)) & fix);
    sbot = sbot << fix;

    return _arf_mul_round_2(z, hi, lo, stop, sbot, sgnbit, prec, rnd, -fix, x, y);
}

int
arf_mul_rnd_down(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec)
{
    slong xn, yn, rn, fix, expfix, zn;
    ulong hi, lo;
    int sgnbit, ret;
    nn_ptr zptr;

    xn = ARF_XSIZE(x);
    yn = ARF_XSIZE(y);
    sgnbit = (xn ^ yn) & 1;
    xn >>= 1;
    yn >>= 1;

    if (yn > xn)
    {
        FLINT_SWAP(arf_srcptr, x, y);
        FLINT_SWAP(slong, xn, yn);
    }

    /* Either operand is a special value. */
    if (yn == 0)
    {
        arf_mul_special(z, x, y);
        return 0;
    }

    /* products of 1- and 2-limb operands: the original register code */
    if (xn == 1)
    {
        lo = ARF_NOPTR_D(x)[0];
        hi = ARF_NOPTR_D(y)[0];

        umul_ppmm(hi, lo, hi, lo);

        /* Shift so that the top bit is set (branch free version). */
        fix = !(hi >> (FLINT_BITS - 1));
        hi = (hi << fix) | ((lo >> (FLINT_BITS - 1)) & fix);
        lo = (lo << fix);

        ARF_DEMOTE(z);

        if (lo == 0)
        {
            zn = 1;

            if (prec >= FLINT_BITS)
            {
                lo = hi;
                ret = 0;
            }
            else
            {
                lo = MASK_LIMB(hi, FLINT_BITS - prec);
                ret = (lo != hi);
            }
        }
        else
        {
            zn = 2;

            if (prec <= FLINT_BITS)  /* Must be inexact. */
            {
                lo = MASK_LIMB(hi, FLINT_BITS - prec);
                zn = ret = 1;
            }
            else if (prec >= 2 * FLINT_BITS) /* Must be exact. */
            {
                ret = 0;
            }
            else  /* prec < FLINT_BITS < 2 * FLINT_BITS */
            {
                ret = MASK_LIMB(lo, 2 * FLINT_BITS - prec) != lo;
                lo = MASK_LIMB(lo, 2 * FLINT_BITS - prec);

                if (lo == 0)
                {
                    zn = 1;
                    lo = hi;
                }
            }
        }

        _fmpz_add2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), -fix);
        ARF_XSIZE(z) = ARF_MAKE_XSIZE(zn, sgnbit);

        zptr = ARF_NOPTR_D(z);
        zptr[0] = lo;
        zptr[1] = hi;

        return ret;
    }
    else if (xn == 2)
    {
        ulong zz[4];
        ulong x1, x0, y1, y0;

        x0 = ARF_NOPTR_D(x)[0];
        x1 = ARF_NOPTR_D(x)[1];

        if (yn == 2)
        {
            y0 = ARF_NOPTR_D(y)[0];
            y1 = ARF_NOPTR_D(y)[1];

            FLINT_MPN_MUL_2X2(zz[3], zz[2], zz[1], zz[0], x1, x0, y1, y0);

            /* Likely case, must be inexact */
            if (prec <= 2 * FLINT_BITS)
            {
                ARF_DEMOTE(z);

                fix = !(zz[3] >> (FLINT_BITS - 1));
                zz[3] = (zz[3] << fix) | ((zz[2] >> (FLINT_BITS - 1)) & fix);
                zz[2] = (zz[2] << fix) | ((zz[1] >> (FLINT_BITS - 1)) & fix);

                _fmpz_add2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), -fix);

                /* Rounding */
                if (prec != 2 * FLINT_BITS)
                {
                    if (prec > FLINT_BITS)
                    {
                        zz[2] &= (LIMB_ONES << (2 * FLINT_BITS - prec));
                    }
                    else if (prec == FLINT_BITS)
                    {
                        zz[2] = 0;
                    }
                    else
                    {
                        zz[3] &= (LIMB_ONES << (FLINT_BITS - prec));
                        zz[2] = 0;
                    }
                }

                if (zz[2] == 0)
                {
                    ARF_XSIZE(z) = ARF_MAKE_XSIZE(1, sgnbit);
                    ARF_NOPTR_D(z)[0] = zz[3];
                }
                else
                {
                    ARF_XSIZE(z) = ARF_MAKE_XSIZE(2, sgnbit);
                    ARF_NOPTR_D(z)[0] = zz[2];
                    ARF_NOPTR_D(z)[1] = zz[3];
                }

                return 1;
            }
        }
        else
        {
            y0 = ARF_NOPTR_D(y)[0];
            FLINT_MPN_MUL_2X1(zz[2], zz[1], zz[0], x1, x0, y0);
        }

        zn = xn + yn;
        ret = _arf_set_round_mpn(z, &expfix, zz, zn, sgnbit, prec, ARF_RND_DOWN);
        _fmpz_add2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), expfix);
        return ret;
    }

    /* small products in full (the common case of ball arithmetic at
       moderate precision: the checks of the high product are not worth
       their instructions here; up to 16 limbs in total every shape has
       at most ARF_MUL_FULL_SMALL limb products) */
    if (xn + yn <= 16 || prec == ARF_PREC_EXACT)
        return _arf_mul_full(z, x, y, xn, yn, sgnbit, prec, ARF_RND_DOWN);

    rn = (prec + FLINT_BITS - 1) / FLINT_BITS;

    if (xn + yn <= rn + 1 || _arf_mul_prefer_full(xn, yn, rn))
        return _arf_mul_full(z, x, y, xn, yn, sgnbit, prec, ARF_RND_DOWN);

    return _arf_mul_high(z, x, y, xn, yn, sgnbit, prec, ARF_RND_DOWN);
}

int
arf_mul_rnd_any(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd)
{
    slong xn, yn, rn;
    int sgnbit;

    if (rnd == ARF_RND_DOWN || rnd == ARF_RND_FAST)
        return arf_mul_rnd_down(z, x, y, prec);

    xn = ARF_XSIZE(x);
    yn = ARF_XSIZE(y);
    sgnbit = (xn ^ yn) & 1;
    xn >>= 1;
    yn >>= 1;

    if (yn > xn)
    {
        FLINT_SWAP(arf_srcptr, x, y);
        FLINT_SWAP(slong, xn, yn);
    }

    /* Either operand is a special value. */
    if (yn == 0)
    {
        arf_mul_special(z, x, y);
        return 0;
    }

    if (rnd == ARF_RND_ACCURATE)
        rnd = ARF_RND_NEAR;

    /* products of one- and two-limb operands rounded in registers, or
       at a higher precision from the full product (exact or nearly so) */
    if (xn <= 2 && prec <= 2 * FLINT_BITS)
        return _arf_mul_small_rnd(z, x, y, xn, yn, sgnbit, prec, rnd);

    if (xn + yn <= 16 || prec == ARF_PREC_EXACT)
        return _arf_mul_full(z, x, y, xn, yn, sgnbit, prec, rnd);

    rn = (prec + FLINT_BITS - 1) / FLINT_BITS;

    if (xn + yn <= rn + 1 || _arf_mul_prefer_full(xn, yn, rn))
        return _arf_mul_full(z, x, y, xn, yn, sgnbit, prec, rnd);

    return _arf_mul_high(z, x, y, xn, yn, sgnbit, prec, rnd);
}

/* x y truncated to rn limbs (see arf.h): the exact product when it fits
   or the shape prefers it, else the normalised high product with its
   guard limb dropped (the deficit below 2^-(FLINT_BITS - 8) ulp) */
void
_arf_mpn_mul_approx(nn_ptr z, slong * zexp, nn_srcptr x, slong xn, nn_srcptr y, slong yn, slong prec)
{
    slong rn = (prec + FLINT_BITS - 1) / FLINT_BITS, zn = xn + yn;

    FLINT_ASSERT(prec >= 1);

    if (zn <= rn)
    {
        /* the exact product, zero padded */
        unsigned int lead;
        FLINT_MPN_MUL_WITH_SPECIAL_CASES(z + rn - zn, x, xn, y, yn);
        _flint_mpn_zero_small(z, rn - zn);
        lead = flint_clz(z[rn - 1]);
        if (lead != 0)
            mpn_lshift(z, z, rn, lead);
        *zexp = -(slong) lead;
    }
    else if (_arf_mul_prefer_full(xn, yn, rn))
    {
        unsigned int lead;
        nn_ptr t;
        ARF_MPN_TMP_DECL

        ARF_MPN_TMP_ALLOC(t, zn)
        FLINT_MPN_MUL_WITH_SPECIAL_CASES(t, x, xn, y, yn);
        lead = flint_clz(t[zn - 1]);
        if (lead != 0)
            mpn_lshift(t + zn - rn, t + zn - rn, rn, lead);
        flint_mpn_copyi(z, t + zn - rn, rn);
        if (lead != 0)
            z[0] |= t[zn - rn - 1] >> (FLINT_BITS - lead);
        *zexp = -(slong) lead;
    }
    else
    {
        mp_limb_pair_t hi = _arf_mpn_mulhigh_normalised(z, x, xn, y, yn, rn);
        *zexp = -(slong) hi.m2;
    }
}

int
arf_mul_mpz(arf_ptr z, arf_srcptr x, const mpz_t y, slong prec, arf_rnd_t rnd)
{
    slong xn, yn;
    slong fix, shift;
    int sgnbit, inexact;

    yn = FLINT_ABS(y->_mp_size);
    xn = ARF_XSIZE(x);
    xn >>= 1;
    sgnbit = ARF_SGNBIT(x) ^ (y->_mp_size < 0);

    /* Either operand is a special value. */
    if (xn == 0 || yn == 0)
    {
        if (arf_is_finite(x))
        {
            arf_zero(z);
        }
        else
        {
            arf_t t;
            arf_init_set_si(t, mpz_sgn(y));
            arf_mul(z, x, t, prec, rnd);
            arf_clear(t);
        }
        return 0;
    }
    else
    {
        slong zn, alloc;
        nn_srcptr xptr, yptr;
        nn_ptr tmp;
        ARF_MUL_TMP_DECL

        ARF_GET_MPN_READONLY(xptr, xn, x);
        yptr = y->_mp_d;

        alloc = zn = xn + yn;
        ARF_MUL_TMP_ALLOC(tmp, alloc)

        FLINT_MPN_MUL_WITH_SPECIAL_CASES(tmp, xptr, xn, yptr, yn);

        shift = yn * FLINT_BITS - (tmp[zn - 1] == 0) * FLINT_BITS;
        zn -= (tmp[zn - 1] == 0);

        inexact = _arf_set_round_mpn(z, &fix, tmp, zn, sgnbit, prec, rnd);

        _fmpz_add_fast(ARF_EXPREF(z), ARF_EXPREF(x), fix + shift);
        ARF_MUL_TMP_FREE(tmp, alloc)

        return inexact;
    }
}
