/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acf.h"
#include "approx_kernel.h"

#define _ACF_KERNEL_STACK 40

/* through the kernel of arf (arf/mpn_complex.c): both parts nonzero,
   finite, with moderate exponents */
static int
_acf_approx_rsqrt_kernel(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x);
    nn_srcptr ad, bd;
    slong an, bn, rn, er, ei;
    ulong stack[2 * _ACF_KERNEL_STACK];
    nn_ptr re, im;

    if (!_acf_part_kernel_ok(a) || !_acf_part_kernel_ok(b) || prec < 2 || prec > (WORD(1) << (FLINT_BITS - 8)))
        return 0;

    ARF_GET_MPN_READONLY(ad, an, a);
    ARF_GET_MPN_READONLY(bd, bn, b);
    rn = (prec + FLINT_BITS - 1) / FLINT_BITS;
    re = (rn <= _ACF_KERNEL_STACK) ? stack : flint_malloc(2 * rn * sizeof(ulong));
    im = re + rn;

    _arf_mpn_crsqrt_approx(re, &er, im, &ei, ad, an, ARF_EXP(a), ARF_SGNBIT(a), bd, bn, ARF_EXP(b), prec);
    _acf_set_part(acf_realref(res), re, rn, er, 0, prec, rnd);
    _acf_set_part(acf_imagref(res), im, rn, ei, !ARF_SGNBIT(b), prec, rnd);

    if (re != stack)
        flint_free(re);
    return 1;
}

void
acf_approx_rsqrt(acf_t y, const acf_t x, slong prec, arf_rnd_t rnd)
{
    if (_acf_approx_rsqrt_kernel(y, x, prec, rnd))
        return;

    /* (the remaining cases, with a zero or extreme part, through the
       square root and one division; rsqrt(0) = +inf as arf_rsqrt) */
    if (arf_is_zero(acf_realref(x)) && arf_is_zero(acf_imagref(x)))
    {
        arf_pos_inf(acf_realref(y));
        arf_zero(acf_imagref(y));
        return;
    }

    {
        acf_t t;
        acf_init(t);
        acf_approx_sqrt(t, x, prec + 4, rnd);
        acf_approx_inv(y, t, prec, rnd);
        acf_clear(t);
    }
}
