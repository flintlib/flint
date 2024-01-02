/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "arf.h"

void __gmpn_div_q(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t, mp_ptr);

void
arf_div_special(arf_t z, const arf_t x, const arf_t y)
{
    if ((arf_is_zero(x) && !arf_is_zero(y) && !arf_is_nan(y)) ||
        (arf_is_inf(y) && !arf_is_special(x)))
    {
        arf_zero(z);
    }
    else if (arf_is_zero(y) || (arf_is_special(x) && arf_is_special(y)) ||
        arf_is_nan(x) || arf_is_nan(y))
    {
        arf_nan(z);
    }
    else if (arf_sgn(x) == arf_sgn(y))
        arf_pos_inf(z);
    else
        arf_neg_inf(z);
}

int
arf_div(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd)
{
    mp_size_t xn, yn, zn, sn, tn, alloc;
    mp_srcptr xptr, yptr;
    mp_ptr tmp;
    mp_ptr tptr, zptr;
    int inexact;
    slong fix, fix2;
    ARF_MUL_TMP_DECL

    if (arf_is_special(x) || arf_is_special(y))
    {
        arf_div_special(z, x, y);
        return 0;
    }

    ARF_GET_MPN_READONLY(xptr, xn, x);
    ARF_GET_MPN_READONLY(yptr, yn, y);

    /* Division by a power of two */
    if (yn == 1 && yptr[0] == LIMB_TOP)
    {
        fmpz_t t;
        fmpz_init_set(t, ARF_EXPREF(y));

        if (ARF_SGNBIT(y))
            inexact = arf_neg_round(z, x, prec, rnd);
        else
            inexact = arf_set_round(z, x, prec, rnd);

        _fmpz_sub2_fast(ARF_EXPREF(z), ARF_EXPREF(z), t, 1);
        fmpz_clear(t);
        return inexact;
    }

    sn = FLINT_MAX(prec - xn * FLINT_BITS + yn * FLINT_BITS, 0) + 32;
    sn = (sn + FLINT_BITS - 1) / FLINT_BITS;

    /* Need space for numerator (tn), quotient (zn), and one extra limb at the
       top of tptr for the multiplication to check the remainder. */
    tn = xn + sn;
    zn = tn - yn + 1;
    alloc = zn + (tn + 1);
/* need tn + 1 extra temporary limbs, which we store at the end of tptr */
    alloc += tn + 1;

    ARF_MUL_TMP_ALLOC(tmp, alloc)

    zptr = tmp;
    tptr = tmp + zn;

    flint_mpn_zero(tptr, sn);
    flint_mpn_copyi(tptr + sn, xptr, xn);
/* uses tn + 1 extra temporary limbs, tn limbs after tptr */
    __gmpn_div_q(zptr, tptr, tn, yptr, yn, tptr + tn);

    if (zptr[zn - 1] == 0)
    {
        zn--;
        fix2 = 0;
    }
    else
    {
        fix2 = FLINT_BITS;
    }

    /* The remainder can only be zero if the last several bits of the
       extended quotient are zero. */
    if ((zptr[0] & ((LIMB_ONE << 24) - 1)) == 0)
    {
        /* The quotient is exact iff multiplying by y gives back the
           input. Note: the multiplication may write sn + xn + 1 limbs, but
           tptr[sn + xn] is guaranteed to be zero in that case since the
           approximate quotient cannot be larger than the true quotient. */
        if (zn >= yn)
            flint_mpn_mul(tptr, zptr, zn, yptr, yn);
        else
            flint_mpn_mul(tptr, yptr, yn, zptr, zn);

        /* The quotient is not exact. Perturbing the approximate quotient
           and rounding gives the correct the result. */
        if (!flint_mpn_zero_p(tptr, sn) ||
            mpn_cmp(tptr + sn, xptr, xn) != 0)
        {
            zptr[0]++;
        }
    }

    inexact = _arf_set_round_mpn(z, &fix, zptr, zn,
            ARF_SGNBIT(x) ^ ARF_SGNBIT(y), prec, rnd);
    _fmpz_sub2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), fix + fix2);
    ARF_MUL_TMP_FREE(tmp, alloc)

    return inexact;
}

