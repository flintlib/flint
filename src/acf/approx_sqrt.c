/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acf.h"
#include "approx_kernel.h"

#define _ACF_KERNEL_STACK 40

/* through the kernel of arf (arf/mpn_complex.c): both parts nonzero, finite, with
   moderate exponents */
static int
_acf_approx_sqrt_kernel(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
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

    _arf_mpn_csqrt_approx(re, &er, im, &ei, ad, an, ARF_EXP(a), ARF_SGNBIT(a), bd, bn, ARF_EXP(b), prec);
    _acf_set_part(acf_realref(res), re, rn, er, 0, prec, rnd);
    _acf_set_part(acf_imagref(res), im, rn, ei, ARF_SGNBIT(b), prec, rnd);

    if (re != stack)
        flint_free(re);
    return 1;
}

void
acf_approx_sqrt(acf_t y, const acf_t x, slong prec, arf_rnd_t rnd)
{
    arf_t r, t, u;
    slong wp;
    int sgn;

#define a acf_realref(x)
#define b acf_imagref(x)
#define c acf_realref(y)
#define d acf_imagref(y)

    if (_acf_approx_sqrt_kernel(y, x, prec, rnd))
        return;

    if (arf_is_zero(b))
    {
        if (arf_sgn(a) >= 0)
        {
            arf_sqrt(c, a, prec, rnd);
            arf_zero(d);
            return;
        }
        else
        {
            arf_neg(d, a);
            arf_sqrt(d, d, prec, rnd);
            arf_zero(c);
            return;
        }
    }

    if (arf_is_zero(a))
    {
        if (arf_sgn(b) >= 0)
        {
            arf_mul_2exp_si(c, b, -1);
            arf_sqrt(c, c, prec, rnd);
            arf_set(d, c);
            return;
        }
        else
        {
            arf_mul_2exp_si(c, b, -1);
            arf_neg(c, c);
            arf_sqrt(c, c, prec, rnd);
            arf_neg(d, c);
            return;
        }
    }

    wp = prec + 4;

    arf_init(r);
    arf_init(t);
    arf_init(u);

    /* r = |a+bi| */
    arf_sosq(r, acf_realref(x), acf_imagref(x), wp, rnd);
    arf_sqrt(r, r, wp, rnd);

    if (arf_sgn(a) >= 0)
    {
        /* sqrt(a+bi) = sqrt((r+a)/2) + b/sqrt(2*(r+a))*i */

        arf_add(t, r, a, wp, rnd);

        arf_mul_2exp_si(u, t, 1);
        arf_sqrt(u, u, wp, rnd);
        arf_div(d, b, u, prec, rnd);

        arf_set_round(c, u, prec, rnd);
        arf_mul_2exp_si(c, c, -1);
    }
    else
    {
        /* sqrt(a+bi) = |b|/sqrt(2*(r-a)) + sgn(b)*sqrt((r-a)/2)*i */

        arf_sub(u, r, a, wp, rnd);

        sgn = arf_sgn(b);

        arf_mul_2exp_si(t, u, 1);
        arf_sqrt(t, t, wp, rnd);
        arf_div(c, b, t, prec, rnd);
        arf_abs(c, c);

        arf_set_round(d, t, prec, rnd);
        arf_mul_2exp_si(d, d, -1);
        if (sgn < 0)
            arf_neg(d, d);
    }

    arf_clear(r);
    arf_clear(t);
    arf_clear(u);

#undef a
#undef b
#undef c
#undef d
}
