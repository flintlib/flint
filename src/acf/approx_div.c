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

/* the parts through the kernels of arf (arf/mpn_complex.c): both parts nonzero, finite,
   with moderate exponents; the results at rn limbs rounded to prec */
#define _ACF_KERNEL_STACK 40

/* below three limbs the real operations (arf_sosq, acf_mul and two
   arf_div) are faster than the complex kernels: up to 1.5 times at one
   limb (measured) */
#define _ACF_DIV_KERNEL_MIN_PREC (2 * FLINT_BITS + 1)

static int
_acf_approx_inv_kernel(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x);
    nn_srcptr ad, bd;
    slong an, bn, rn, er, ei;
    ulong stack[2 * _ACF_KERNEL_STACK];
    nn_ptr re, im;

    if (prec < _ACF_DIV_KERNEL_MIN_PREC || prec > (WORD(1) << (FLINT_BITS - 8)) ||
        !_acf_part_kernel_ok(a) || !_acf_part_kernel_ok(b))
        return 0;

    ARF_GET_MPN_READONLY(ad, an, a);
    ARF_GET_MPN_READONLY(bd, bn, b);
    rn = (prec + FLINT_BITS - 1) / FLINT_BITS;
    re = (rn <= _ACF_KERNEL_STACK) ? stack : flint_malloc(2 * rn * sizeof(ulong));
    im = re + rn;

    _arf_mpn_cinv_approx(re, &er, im, &ei, ad, an, ARF_EXP(a), bd, bn, ARF_EXP(b), prec);
    _acf_set_part(acf_realref(res), re, rn, er, ARF_SGNBIT(a), prec, rnd);
    _acf_set_part(acf_imagref(res), im, rn, ei, !ARF_SGNBIT(b), prec, rnd);

    if (re != stack)
        flint_free(re);
    return 1;
}

static int
_acf_approx_div_kernel(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x), c = acf_realref(y), d = acf_imagref(y);
    nn_srcptr ad, bd, cd, dd;
    slong an, bn, cn, dn, rn, er, ei;
    ulong stack[2 * _ACF_KERNEL_STACK];
    nn_ptr re, im;
    int rneg, ineg;

    if (prec < _ACF_DIV_KERNEL_MIN_PREC || prec > (WORD(1) << (FLINT_BITS - 8)) ||
        !_acf_part_kernel_ok(a) || !_acf_part_kernel_ok(b) || !_acf_part_kernel_ok(c) || !_acf_part_kernel_ok(d))
        return 0;

    ARF_GET_MPN_READONLY(ad, an, a);
    ARF_GET_MPN_READONLY(bd, bn, b);
    ARF_GET_MPN_READONLY(cd, cn, c);
    ARF_GET_MPN_READONLY(dd, dn, d);
    rn = (prec + FLINT_BITS - 1) / FLINT_BITS;
    re = (rn <= _ACF_KERNEL_STACK) ? stack : flint_malloc(2 * rn * sizeof(ulong));
    im = re + rn;

    _arf_mpn_cdiv_approx(re, &er, &rneg, im, &ei, &ineg,
        ad, an, ARF_EXP(a), ARF_SGNBIT(a), bd, bn, ARF_EXP(b), ARF_SGNBIT(b),
        cd, cn, ARF_EXP(c), ARF_SGNBIT(c), dd, dn, ARF_EXP(d), ARF_SGNBIT(d), prec);
    _acf_set_part(acf_realref(res), re, rn, er, rneg, prec, rnd);
    _acf_set_part(acf_imagref(res), im, rn, ei, ineg, prec, rnd);

    if (re != stack)
        flint_free(re);
    return 1;
}

void
acf_approx_inv(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
{
    if (_acf_approx_inv_kernel(res, x, prec, rnd))
        return;

    if (arf_is_zero(acf_imagref(x)))
    {
        /* todo: arf_inv */
        arf_ui_div(acf_realref(res), 1, acf_realref(x), prec, rnd);
        arf_zero(acf_imagref(res));
    }
    else if (arf_is_zero(acf_realref(x)))
    {
        /* todo: arf_inv */
        arf_si_div(acf_imagref(res), -1, acf_imagref(x), prec, rnd);
        arf_zero(acf_realref(res));
    }
    else
    {
        arf_t t;
        arf_init(t);

        arf_sosq(t, acf_realref(x), acf_imagref(x), prec, rnd);
        arf_div(acf_realref(res), acf_realref(x), t, prec, rnd);
        arf_div(acf_imagref(res), acf_imagref(x), t, prec, rnd);
        arf_neg(acf_imagref(res), acf_imagref(res));

        arf_clear(t);
    }
}

#define a acf_realref(x)
#define b acf_imagref(x)
#define c acf_realref(y)
#define d acf_imagref(y)

void
acf_approx_div(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)
{
    if (_acf_approx_div_kernel(res, x, y, prec, rnd))
        return;

    if (arf_is_zero(d))
    {
        if (arf_is_zero(b))
        {
            arf_div(acf_realref(res), a, c, prec, rnd);
            arf_zero(acf_imagref(res));
        }
        else if (arf_is_zero(a))
        {
            arf_div(acf_imagref(res), b, c, prec, rnd);
            arf_zero(acf_realref(res));
        }
        else if (res != y)
        {
            arf_div(acf_realref(res), a, c, prec, rnd);
            arf_div(acf_imagref(res), b, c, prec, rnd);
        }
        else
        {
            arf_t t;
            arf_init(t);
            arf_set(t, c);
            arf_div(acf_realref(res), a, t, prec, rnd);
            arf_div(acf_imagref(res), b, t, prec, rnd);
            arf_clear(t);
        }
    }
    else if (arf_is_zero(c))
    {
        if (arf_is_zero(b))
        {
            arf_div(acf_imagref(res), a, d, prec, rnd);
            arf_neg(acf_imagref(res), acf_imagref(res));
            arf_zero(acf_realref(res));
        }
        else if (arf_is_zero(a))
        {
            arf_div(acf_realref(res), b, d, prec, rnd);
            arf_zero(acf_imagref(res));
        }
        else if (res != y)
        {
            arf_div(acf_realref(res), a, d, prec, rnd);
            arf_div(acf_imagref(res), b, d, prec, rnd);
            arf_swap(acf_realref(res), acf_imagref(res));
            arf_neg(acf_imagref(res), acf_imagref(res));
        }
        else
        {
            arf_t t;
            arf_init(t);
            arf_set(t, d);
            arf_div(acf_realref(res), a, t, prec, rnd);
            arf_div(acf_imagref(res), b, t, prec, rnd);
            arf_swap(acf_realref(res), acf_imagref(res));
            arf_neg(acf_imagref(res), acf_imagref(res));
            arf_clear(t);
        }
    }
    else
    {
        arf_t t;
        acf_t u;

        arf_init(t);

        arf_sosq(t, acf_realref(y), acf_imagref(y), prec, rnd);

        /* (res may be y: the conjugate as a copy then) */
        if (res != y)
        {
            arf_init_set_shallow(acf_realref(u), acf_realref(y));
            arf_init_neg_shallow(acf_imagref(u), acf_imagref(y));
            acf_mul(res, x, u, prec, rnd);
        }
        else
        {
            acf_init(u);
            arf_set(acf_realref(u), acf_realref(y));
            arf_neg(acf_imagref(u), acf_imagref(y));
            acf_mul(res, x, u, prec, rnd);
            acf_clear(u);
        }

        arf_div(acf_realref(res), acf_realref(res), t, prec, rnd);
        arf_div(acf_imagref(res), acf_imagref(res), t, prec, rnd);

        arf_clear(t);
    }
}

#undef a
#undef b
#undef c
#undef d
