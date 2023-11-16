/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_calc.h"

int arb_calc_refine_root_newton(arb_t r, arb_calc_func_t func,
    void * param, const arb_t start, const arb_t conv_region,
    const arf_t conv_factor, slong eval_extra_prec, slong prec)
{
    slong precs[FLINT_BITS];
    slong i, iters, wp, padding, start_prec;
    int result;

    start_prec = arb_rel_accuracy_bits(start);

    if (arb_calc_verbose)
        flint_printf("newton initial accuracy: %wd\n", start_prec);

    padding = arf_abs_bound_lt_2exp_si(conv_factor);
    padding = FLINT_MIN(padding, prec) + 5;
    padding = FLINT_MAX(0, padding);

    precs[0] = prec + padding;
    iters = 1;
    while ((iters < FLINT_BITS) && (precs[iters-1] + padding > 2*start_prec))
    {
        precs[iters] = (precs[iters-1] / 2) + padding;
        iters++;

        if (iters == FLINT_BITS)
        {
            return ARB_CALC_IMPRECISE_INPUT;
        }
    }

    arb_set(r, start);

    for (i = iters - 1; i >= 0; i--)
    {
        wp = precs[i] + eval_extra_prec;

        if (arb_calc_verbose)
            flint_printf("newton step: wp = %wd + %wd = %wd\n",
                precs[i], eval_extra_prec, wp);

        if ((result = arb_calc_newton_step(r, func, param,
            r, conv_region, conv_factor, wp)) != ARB_CALC_SUCCESS)
        {
            return result;
        }
    }

    return ARB_CALC_SUCCESS;
}

