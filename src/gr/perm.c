/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "perm.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "gr.h"

#define PERM_N(ctx) (*((ulong *) (ctx)))

typedef struct
{
    slong * entries;
}
perm_struct;

typedef perm_struct perm_t[1];


int _gr_perm_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Symmetric group S_");
    gr_stream_write_ui(out, PERM_N(ctx));
    gr_stream_write(out, " (perm)");
    return GR_SUCCESS;
}

int
_gr_perm_init(perm_t res, gr_ctx_t ctx)
{
    res->entries = _perm_init(PERM_N(ctx));
    return GR_SUCCESS;
}

int
_gr_perm_clear(perm_t res, gr_ctx_t ctx)
{
    _perm_clear(res->entries);
    return GR_SUCCESS;
}

void
_gr_perm_swap(perm_t x, perm_t y, gr_ctx_t ctx)
{
    perm_struct t = *x;
    *x = *y;
    *y = t;
}

int
_gr_perm_write(gr_stream_t out, perm_t x, gr_ctx_t ctx)
{
    slong i;

    gr_stream_write(out, "[");

    for (i = 0; i < PERM_N(ctx); i++)
    {
        gr_stream_write_si(out, x->entries[i]);
        if (i + 1 < PERM_N(ctx))
            gr_stream_write(out, ", ");
    }

    gr_stream_write(out, "]");

    return GR_SUCCESS;
}

int
_gr_perm_randtest(perm_t res, flint_rand_t state, gr_ctx_t ctx)
{
    _perm_randtest(res->entries, PERM_N(ctx), state);
    return GR_SUCCESS;
}

truth_t
_gr_perm_equal(const perm_t x, const perm_t y, gr_ctx_t ctx)
{
    return _perm_equal(x->entries, y->entries, PERM_N(ctx)) ? T_TRUE : T_FALSE;
}

int
_gr_perm_set(perm_t res, const perm_t x, gr_ctx_t ctx)
{
    _perm_set(res->entries, x->entries, PERM_N(ctx));
    return GR_SUCCESS;
}

int
_gr_perm_set_other(perm_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    if (x_ctx->which_ring == GR_CTX_PERM)
    {
        if (PERM_N(x_ctx) == PERM_N(ctx))
            return _gr_perm_set(res, x, ctx);

        return GR_DOMAIN;
    }
    else if (x_ctx->which_ring == GR_CTX_GR_MAT && MATRIX_CTX(x_ctx)->base_ring->which_ring == GR_CTX_FMPZ)
    {
        slong i, j, n, c;
        const fmpz_mat_struct * mat = x;

        n = PERM_N(ctx);

        /* todo: factor out mat_is_permutation method */

        if (fmpz_mat_nrows(mat) != n || fmpz_mat_ncols(mat) != n)
            return GR_DOMAIN;

        for (i = 0; i < n; i++)
        {
            c = 0;
            for (j = 0; j < n; j++)
            {
                if (fmpz_is_zero(fmpz_mat_entry(mat, i, j)))
                    continue;

                c++;
                if (c != 1 || !fmpz_is_one(fmpz_mat_entry(mat, i, j)))
                    return GR_DOMAIN;
            }

            if (c == 0)
                return GR_DOMAIN;
        }

        for (i = 0; i < n; i++)
        {
            c = 0;
            for (j = 0; j < n; j++)
            {
                if (fmpz_is_zero(fmpz_mat_entry(mat, j, i)))
                    continue;

                c++;
                if (c != 1)
                    return GR_DOMAIN;
            }

            if (c == 0)
                return GR_DOMAIN;
        }

        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
            {
                if (fmpz_is_one(fmpz_mat_entry(mat, i, j)))
                {
                    res->entries[i] = j;
                    break;
                }
            }

        return GR_SUCCESS;
    }
    else
    {
        return GR_UNABLE;
    }
}

int
_gr_perm_one(perm_t res, gr_ctx_t ctx)
{
    _perm_one(res->entries, PERM_N(ctx));
    return GR_SUCCESS;
}

truth_t
_gr_perm_is_one(const perm_t x, gr_ctx_t ctx)
{
    slong i, n = PERM_N(ctx);

    for (i = 0; i < n; i++)
        if (x->entries[i] != i)
            return T_FALSE;

    return T_UNKNOWN;
}

int
_gr_perm_mul(perm_t res, const perm_t x, const perm_t y, gr_ctx_t ctx)
{
    _perm_compose(res->entries, x->entries, y->entries, PERM_N(ctx));
    return GR_SUCCESS;
}

/* todo: should be generic. also want left division */
int
_gr_perm_div(perm_t res, const perm_t x, const perm_t y, gr_ctx_t ctx)
{
    slong n = PERM_N(ctx);
    slong * t;
    t = _perm_init(n);
    _perm_inv(t, y->entries, n);
    _perm_compose(res->entries, x->entries, t, n);
    _perm_clear(t);
    return GR_SUCCESS;
}

int
_gr_perm_inv(perm_t res, const perm_t x, gr_ctx_t ctx)
{
    _perm_inv(res->entries, x->entries, PERM_N(ctx));
    return GR_SUCCESS;
}

/* todo: parity */

int _perm_methods_initialized = 0;

gr_static_method_table _perm_methods;

gr_method_tab_input _perm_methods_input[] =
{
    {GR_METHOD_CTX_IS_FINITE,
                            (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_MULTIPLICATIVE_GROUP,
                            (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) _gr_perm_ctx_write},
    {GR_METHOD_INIT,        (gr_funcptr) _gr_perm_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) _gr_perm_clear},
    {GR_METHOD_SWAP,        (gr_funcptr) _gr_perm_swap},
    {GR_METHOD_RANDTEST,    (gr_funcptr) _gr_perm_randtest},
    {GR_METHOD_WRITE,       (gr_funcptr) _gr_perm_write},
    {GR_METHOD_ONE,         (gr_funcptr) _gr_perm_one},
    {GR_METHOD_EQUAL,       (gr_funcptr) _gr_perm_equal},
    {GR_METHOD_SET,         (gr_funcptr) _gr_perm_set},
    {GR_METHOD_SET_OTHER,   (gr_funcptr) _gr_perm_set_other},
    {GR_METHOD_MUL,         (gr_funcptr) _gr_perm_mul},
    {GR_METHOD_INV,         (gr_funcptr) _gr_perm_inv},
    {GR_METHOD_DIV,         (gr_funcptr) _gr_perm_div},
    {0,                     (gr_funcptr) NULL},
};

void
gr_ctx_init_perm(gr_ctx_t ctx, ulong n)
{
    ctx->which_ring = GR_CTX_PERM;
    ctx->sizeof_elem = sizeof(perm_struct);
    ctx->size_limit = WORD_MAX;

    PERM_N(ctx) = n;

    ctx->methods = _perm_methods;

    if (!_perm_methods_initialized)
    {
        gr_method_tab_init(_perm_methods, _perm_methods_input);
        _perm_methods_initialized = 1;
    }
}
