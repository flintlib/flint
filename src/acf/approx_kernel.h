/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACF_APPROX_KERNEL_H
#define ACF_APPROX_KERNEL_H

/* The complex kernels of arf (arf/mpn_complex.c) for acf_approx_inv / div / sqrt: the
   mantissa/exponent pairs of the parts, and the parts back from the rn
   limbs the kernels wrote. */

#include "acf.h"
#include "arf.h"

/* exponents for which the kernels' slong arithmetic (sums of two
   exponents, doubled) cannot overflow */
#define _ACF_KERNEL_EXP_MODERATE(e) (FLINT_ABS(e) < (WORD(1) << (FLINT_BITS - 8)))

/* whether the part is nonzero, finite and of moderate exponent */
FLINT_FORCE_INLINE int
_acf_part_kernel_ok(arf_srcptr x)
{
    return !ARF_IS_SPECIAL(x) && !COEFF_IS_MPZ(ARF_EXP(x)) && _ACF_KERNEL_EXP_MODERATE(ARF_EXP(x));
}

/* the part res from the rn limbs d a kernel wrote (clobbered), rounded to
   prec bits in rnd, with its exponent and sign bit; exp = WORD_MIN: zero */
FLINT_FORCE_INLINE void
_acf_set_part(arf_ptr res, nn_ptr d, slong rn, slong exp, int sgnbit, slong prec, arf_rnd_t rnd)
{
    slong fix, zn;
    arf_rnd_t mrnd;

    if (exp == WORD_MIN)
    {
        arf_zero(res);
        return;
    }

    mrnd = _arf_rnd_resolve(rnd, sgnbit);
    _arf_mpn_set_round(d, &fix, d, rn, prec, mrnd);
    zn = (prec + FLINT_BITS - 1) / FLINT_BITS;
    _arf_set_mpn_canonical_copy(res, d, zn, sgnbit);
    fmpz_set_si(ARF_EXPREF(res), exp + fix);
}

#endif
