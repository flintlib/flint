/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_ext.h"
#include "ca_vec.h"

static ulong hash_func(calcium_func_code func, ca_srcptr args, slong nargs, ca_ctx_t ctx)
{
    ulong s;
    slong i;

    s = func;

    for (i = 0; i < nargs; i++)
        s = ca_hash_repr(args + i, ctx) * 1000003 + s;

    return s;
}

void
ca_ext_init_qqbar(ca_ext_t res, const qqbar_t x, ca_ctx_t ctx)
{
    fmpq_poly_t t;

    res->head = CA_QQBar;

    qqbar_init(CA_EXT_QQBAR(res));
    qqbar_set(CA_EXT_QQBAR(res), x);

    /* qqbar_cache_enclosure(CA_EXT_QQBAR(res), 2 * QQBAR_DEFAULT_PREC); */

    /* nf_init wants an fmpq_poly_t, so mock up one */
    t->coeffs = QQBAR_POLY(x)->coeffs;
    t->den[0] = 1;
    t->length = QQBAR_POLY(x)->length;
    t->alloc = QQBAR_POLY(x)->alloc;

    CA_EXT_QQBAR_NF(res) = flint_malloc(sizeof(nf_struct));
    nf_init(CA_EXT_QQBAR_NF(res), t);

    res->hash = qqbar_hash(CA_EXT_QQBAR(res));
    res->depth = 0;
}

slong ca_depth(const ca_t x, ca_ctx_t ctx);

static void _ca_ext_init_func(ca_ext_t res, ca_ctx_t ctx)
{
    slong i, d;

    CA_EXT_FUNC_PREC(res) = 0;
    acb_init(CA_EXT_FUNC_ENCLOSURE(res));
    acb_indeterminate(CA_EXT_FUNC_ENCLOSURE(res));

    res->hash = hash_func(CA_EXT_HEAD(res), CA_EXT_FUNC_ARGS(res), CA_EXT_FUNC_NARGS(res), ctx);
    res->depth = 0;

    for (i = 0; i < CA_EXT_FUNC_NARGS(res); i++)
    {
        d = ca_depth(CA_EXT_FUNC_ARGS(res) + i, ctx);
        res->depth = FLINT_MAX(res->depth, d + 1);
    }

    res->data.func_data.qqbar = NULL;
}

void
ca_ext_init_const(ca_ext_t res, calcium_func_code func, ca_ctx_t ctx)
{
    res->head = func;

    CA_EXT_FUNC_NARGS(res) = 0;
    CA_EXT_FUNC_ARGS(res) = NULL;

    _ca_ext_init_func(res, ctx);
}

void
ca_ext_init_fx(ca_ext_t res, calcium_func_code func, const ca_t x, ca_ctx_t ctx)
{
    res->head = func;
    res->hash = 0;

    CA_EXT_FUNC_NARGS(res) = 1;
    CA_EXT_FUNC_ARGS(res) = _ca_vec_init(1, ctx);
    ca_set(CA_EXT_FUNC_ARGS(res), x, ctx);

    _ca_ext_init_func(res, ctx);
}

void
ca_ext_init_fxy(ca_ext_t res, calcium_func_code func, const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    res->head = func;
    res->hash = 0;

    CA_EXT_FUNC_NARGS(res) = 2;
    CA_EXT_FUNC_ARGS(res) = _ca_vec_init(2, ctx);
    ca_set(CA_EXT_FUNC_ARGS(res), x, ctx);
    ca_set(CA_EXT_FUNC_ARGS(res) + 1, y, ctx);

    _ca_ext_init_func(res, ctx);
}

void
ca_ext_init_fxn(ca_ext_t res, calcium_func_code func, ca_srcptr x, slong nargs, ca_ctx_t ctx)
{
    res->head = func;

    CA_EXT_FUNC_NARGS(res) = nargs;

    if (nargs == 0)
        CA_EXT_FUNC_ARGS(res) = NULL;
    else
        CA_EXT_FUNC_ARGS(res) = _ca_vec_init(nargs, ctx);

    _ca_vec_set(CA_EXT_FUNC_ARGS(res), x, nargs, ctx);

    _ca_ext_init_func(res, ctx);
}

