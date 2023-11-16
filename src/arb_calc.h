/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ARB_CALC_H
#define ARB_CALC_H

#include "arb.h"

#ifdef __cplusplus
extern "C" {
#endif

extern FLINT_TLS_PREFIX int arb_calc_verbose;

typedef int (*arb_calc_func_t)(arb_ptr out,
    const arb_t inp, void * param, slong order, slong prec);

#define ARB_CALC_SUCCESS 0
#define ARB_CALC_IMPRECISE_INPUT 1
#define ARB_CALC_NO_CONVERGENCE 2

/* Root-finding */

static inline void
arf_interval_init(arf_interval_t v)
{
    arf_init(&v->a);
    arf_init(&v->b);
}

static inline void
arf_interval_clear(arf_interval_t v)
{
    arf_clear(&v->a);
    arf_clear(&v->b);
}

static inline arf_interval_ptr
_arf_interval_vec_init(slong n)
{
    slong i;
    arf_interval_ptr v = (arf_interval_ptr) flint_malloc(sizeof(arf_interval_struct) * n);

    for (i = 0; i < n; i++)
        arf_interval_init(v + i);

    return v;
}

static inline void
_arf_interval_vec_clear(arf_interval_ptr v, slong n)
{
    slong i;
    for (i = 0; i < n; i++)
        arf_interval_clear(v + i);
    flint_free(v);
}

static inline void
arf_interval_set(arf_interval_t v, const arf_interval_t u)
{
    arf_set(&v->a, &u->a);
    arf_set(&v->b, &u->b);
}

static inline void
arf_interval_swap(arf_interval_t v, arf_interval_t u)
{
    arf_swap(&v->a, &u->a);
    arf_swap(&v->b, &u->b);
}

static inline void
arf_interval_get_arb(arb_t x, const arf_interval_t v, slong prec)
{
    arb_set_interval_arf(x, &v->a, &v->b, prec);
}

#ifdef FLINT_HAVE_FILE
void arf_interval_fprintd(FILE * file, const arf_interval_t v, slong n);
#endif

void arf_interval_printd(const arf_interval_t v, slong n);

/* bisection */

int arb_calc_partition(arf_interval_t L, arf_interval_t R,
    arb_calc_func_t func, void * param, const arf_interval_t block, slong prec);

slong arb_calc_isolate_roots(arf_interval_ptr * blocks, int ** flags,
    arb_calc_func_t func, void * param,
    const arf_interval_t block, slong maxdepth, slong maxeval, slong maxfound,
    slong prec);

int arb_calc_refine_root_bisect(arf_interval_t r, arb_calc_func_t func,
    void * param, const arf_interval_t start, slong iter, slong prec);

/* newton iteration */

void arb_calc_newton_conv_factor(arf_t conv_factor,
    arb_calc_func_t func, void * param, const arb_t conv_region, slong prec);

int arb_calc_newton_step(arb_t xnew, arb_calc_func_t func, void * param,
    const arb_t x, const arb_t conv_region, const arf_t conv_factor, slong prec);

int arb_calc_refine_root_newton(arb_t r, arb_calc_func_t func,
    void * param, const arb_t start, const arb_t conv_region,
    const arf_t conv_factor, slong eval_extra_prec, slong prec);



#ifdef __cplusplus
}
#endif

#endif

