/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "acb_modular.h"
#include "gr.h"

int _gr_psl2z_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Modular group (psl2z)");
    return GR_SUCCESS;
}

int
_gr_psl2z_init(psl2z_t res, gr_ctx_t ctx)
{
    psl2z_init(res);
    return GR_SUCCESS;
}

int
_gr_psl2z_clear(psl2z_t res, gr_ctx_t ctx)
{
    psl2z_clear(res);
    return GR_SUCCESS;
}

void
_gr_psl2z_swap(psl2z_t x, psl2z_t y, gr_ctx_t ctx)
{
    psl2z_swap(x, y);
}

int
_gr_psl2z_write(gr_stream_t out, psl2z_t x, gr_ctx_t ctx)
{
    gr_stream_write(out, "[[");
    gr_stream_write_fmpz(out, &x->a);
    gr_stream_write(out, ", ");
    gr_stream_write_fmpz(out, &x->b);
    gr_stream_write(out, "], [");
    gr_stream_write_fmpz(out, &x->c);
    gr_stream_write(out, ", ");
    gr_stream_write_fmpz(out, &x->d);
    gr_stream_write(out, "]]");
    return GR_SUCCESS;
}

int
_gr_psl2z_randtest(psl2z_t res, flint_rand_t state, gr_ctx_t ctx)
{
    if (n_randint(state, 4))
        psl2z_randtest(res, state, 8);
    else
        psl2z_randtest(res, state, 1 + n_randint(state, 100));
    return GR_SUCCESS;
}

truth_t
_gr_psl2z_equal(const psl2z_t x, const psl2z_t y, gr_ctx_t ctx)
{
    return psl2z_equal(x, y) ? T_TRUE : T_FALSE;
}

int
_gr_psl2z_set(psl2z_t res, const psl2z_t x, gr_ctx_t ctx)
{
    psl2z_set(res, x);
    return GR_SUCCESS;
}

int
_gr_psl2z_set_other(psl2z_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    if (x_ctx->which_ring == GR_CTX_PSL2Z)
    {
        psl2z_set(res, x);
        return GR_SUCCESS;
    }
    else if (x_ctx->which_ring == GR_CTX_GR_MAT && MATRIX_CTX(x_ctx)->base_ring->which_ring == GR_CTX_FMPZ)
    {
        fmpz_t det;

        if (fmpz_mat_nrows(x) != 2 || fmpz_mat_ncols(x) != 2)
            return GR_DOMAIN;

        fmpz_init(det);
        fmpz_mat_det(det, x);
        if (fmpz_is_one(det))
        {
            fmpz_set(&res->a, fmpz_mat_entry(x, 0, 0));
            fmpz_set(&res->b, fmpz_mat_entry(x, 0, 1));
            fmpz_set(&res->c, fmpz_mat_entry(x, 1, 0));
            fmpz_set(&res->d, fmpz_mat_entry(x, 1, 1));

            if (fmpz_sgn(&res->c) < 0 || (fmpz_is_zero(&res->c) && fmpz_sgn(&res->d) < 0))
            {
                fmpz_neg(&res->a, &res->a);
                fmpz_neg(&res->b, &res->b);
                fmpz_neg(&res->c, &res->c);
                fmpz_neg(&res->d, &res->d);
            }
            fmpz_clear(det);
            return GR_SUCCESS;
        }
        else
        {
            fmpz_clear(det);
            return GR_DOMAIN;
        }
    }
    else
    {
        return GR_UNABLE;
    }
}

int
_gr_psl2z_one(psl2z_t res, gr_ctx_t ctx)
{
    psl2z_one(res);
    return GR_SUCCESS;
}

truth_t
_gr_psl2z_is_one(const psl2z_t x, gr_ctx_t ctx)
{
    return psl2z_is_one(x) ? T_TRUE : T_FALSE;
}

int
_gr_psl2z_mul(psl2z_t res, const psl2z_t x, const psl2z_t y, gr_ctx_t ctx)
{
    psl2z_mul(res, x, y);
    return GR_SUCCESS;
}

/* todo: should be generic. also want left division */
int
_gr_psl2z_div(psl2z_t res, const psl2z_t x, const psl2z_t y, gr_ctx_t ctx)
{
    psl2z_t t;
    psl2z_init(t);
    psl2z_inv(t, y);
    psl2z_mul(res, x, t);
    psl2z_clear(t);
    return GR_SUCCESS;
}

int
_gr_psl2z_inv(psl2z_t res, const psl2z_t x, gr_ctx_t ctx)
{
    psl2z_inv(res, x);
    return GR_SUCCESS;
}


int _psl2z_methods_initialized = 0;

gr_static_method_table _psl2z_methods;

gr_method_tab_input _psl2z_methods_input[] =
{
    {GR_METHOD_CTX_IS_FINITE,
                            (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_MULTIPLICATIVE_GROUP,
                            (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) _gr_psl2z_ctx_write},
    {GR_METHOD_INIT,        (gr_funcptr) _gr_psl2z_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) _gr_psl2z_clear},
    {GR_METHOD_SWAP,        (gr_funcptr) _gr_psl2z_swap},
    {GR_METHOD_RANDTEST,    (gr_funcptr) _gr_psl2z_randtest},
    {GR_METHOD_WRITE,       (gr_funcptr) _gr_psl2z_write},
    {GR_METHOD_ONE,         (gr_funcptr) _gr_psl2z_one},
    {GR_METHOD_EQUAL,       (gr_funcptr) _gr_psl2z_equal},
    {GR_METHOD_SET,         (gr_funcptr) _gr_psl2z_set},
    {GR_METHOD_SET_OTHER,   (gr_funcptr) _gr_psl2z_set_other},
    {GR_METHOD_MUL,         (gr_funcptr) _gr_psl2z_mul},
    {GR_METHOD_INV,         (gr_funcptr) _gr_psl2z_inv},
    {GR_METHOD_DIV,         (gr_funcptr) _gr_psl2z_div},
    {0,                     (gr_funcptr) NULL},
};

void
gr_ctx_init_psl2z(gr_ctx_t ctx)
{
    ctx->which_ring = GR_CTX_PSL2Z;
    ctx->sizeof_elem = sizeof(psl2z_struct);
    ctx->size_limit = WORD_MAX;

    ctx->methods = _psl2z_methods;

    if (!_psl2z_methods_initialized)
    {
        gr_method_tab_init(_psl2z_methods, _psl2z_methods_input);
        _psl2z_methods_initialized = 1;
    }
}
