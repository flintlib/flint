/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

static inline slong _arf_mag(const arf_t c)
{
    slong m = arf_abs_bound_lt_2exp_si(c);
    return FLINT_MAX(m, 0);
}

void
_arb_poly_newton_refine_root(arb_t r, arb_srcptr poly, slong len,
    const arb_t start,
    const arb_t convergence_interval,
    const arf_t convergence_factor,
    slong eval_extra_prec,
    slong prec)
{
    slong precs[FLINT_BITS];
    slong i, iters, wp, padding, start_prec;

    start_prec = arb_rel_accuracy_bits(start);

    padding = 5 + _arf_mag(convergence_factor);
    precs[0] = prec + padding;
    iters = 1;
    while ((iters < FLINT_BITS) && (precs[iters-1] + padding > 2*start_prec))
    {
        precs[iters] = (precs[iters-1] / 2) + padding;
        iters++;

        if (iters == FLINT_BITS)
        {
            flint_throw(FLINT_ERROR, "newton_refine_root: initial value too imprecise\n");
        }
    }

    arb_set(r, start);

    for (i = iters - 1; i >= 0; i--)
    {
        wp = precs[i] + eval_extra_prec;

        if (!_arb_poly_newton_step(r, poly, len, r, convergence_interval,
            convergence_factor, wp))
        {
            flint_printf("warning: newton_refine_root: improvement failed\n");
            break;
        }

    }
}

