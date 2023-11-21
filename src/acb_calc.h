/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_CALC_H
#define ACB_CALC_H

#include "acb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef int (*acb_calc_func_t)(acb_ptr out,
    const acb_t inp, void * param, slong order, slong prec);

/* Integration (old) */

void acb_calc_cauchy_bound(arb_t bound, acb_calc_func_t func,
    void * param, const acb_t x, const arb_t radius,
    slong maxdepth, slong prec);

int acb_calc_integrate_taylor(acb_t res,
    acb_calc_func_t func, void * param,
    const acb_t a, const acb_t b,
    const arf_t inner_radius,
    const arf_t outer_radius,
    slong accuracy_goal, slong prec);

/* Integration */

typedef struct
{
    slong deg_limit;
    slong eval_limit;
    slong depth_limit;
    int use_heap;
    int verbose;
}
acb_calc_integrate_opt_struct;

typedef acb_calc_integrate_opt_struct acb_calc_integrate_opt_t[1];

void acb_calc_integrate_opt_init(acb_calc_integrate_opt_t options);

int
acb_calc_integrate(acb_t res, acb_calc_func_t f, void * param,
    const acb_t a, const acb_t b,
    slong goal, const mag_t tol,
    const acb_calc_integrate_opt_t options,
    slong prec);

int
acb_calc_integrate_gl_auto_deg(acb_t res, slong * eval_count,
    acb_calc_func_t f, void * param,
    const acb_t a, const acb_t b, const mag_t tol,
    slong deg_limit, int verbose, slong prec);

#ifdef __cplusplus
}
#endif

#endif

