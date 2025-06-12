/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"

/*
Todo:
* Generalize to rectangular matrices
* Parallel versions
* Make a version for integers optimized for small entries
* Add an option for sign reversals in the internal cofactor expansion
  and call from gr_mat_det_cofactor
* Exploit sparsity in the cofactor expansion

*/

static int
_gr_mat_permanent_cofactor_recursive(gr_ptr, const gr_mat_t, slong *, slong, gr_ctx_t);

static int
_gr_mat_permanent_cofactor_recursive(gr_ptr res, const gr_mat_t A, slong * rows, slong col, gr_ctx_t ctx)
{
    slong n = A->r - col;
    int status = GR_SUCCESS;
    slong r[16];
    slong i, j;

#define E(ii, jj) gr_mat_entry_srcptr(A, rows[ii], col + (jj), ctx)

    gr_ptr t;
    GR_TMP_INIT(t, ctx);

    if (n == 3)
    {
        status |= gr_mul(t, E(1, 0), E(2, 1), ctx);
        status |= gr_addmul(t, E(1, 1), E(2, 0), ctx);
        status |= gr_mul(res, t, E(0, 2), ctx);

        status |= gr_mul(t, E(1, 2), E(2, 0), ctx);
        status |= gr_addmul(t, E(1, 0), E(2, 2), ctx);
        status |= gr_addmul(res, t, E(0, 1), ctx);

        status |= gr_mul(t, E(1, 1), E(2, 2), ctx);
        status |= gr_addmul(t, E(1, 2), E(2, 1), ctx);
        status |= gr_addmul(res, t, E(0, 0), ctx);
    }
    else
    {
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < i; j++)
                r[j] = rows[j];
            for (j = i + 1; j < n; j++)
                r[j - 1] = rows[j];

            status |= _gr_mat_permanent_cofactor_recursive(t, A, r, col + 1, ctx);

            if (i == 0)
                status |= gr_mul(res, E(i, 0), t, ctx);
            else
                status |= gr_addmul(res, E(i, 0), t, ctx);
        }
    }

#undef E

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int
gr_mat_permanent_cofactor(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx)
{
    slong n = A->r;

    if (n > 16)
        return GR_UNABLE;

    if (n == 0)
    {
        return gr_one(res, ctx);
    }
    else if (n == 1)
    {
        return gr_set(res, A->entries, ctx);
    }
    else if (n == 2)
    {
        int status = GR_SUCCESS;
        status |= gr_mul(res, gr_mat_entry_srcptr(A, 0, 0, ctx),
                               gr_mat_entry_srcptr(A, 1, 1, ctx), ctx);
        status |= gr_addmul(res, gr_mat_entry_srcptr(A, 0, 1, ctx),
                                 gr_mat_entry_srcptr(A, 1, 0, ctx), ctx);
        return status;
    }
    else
    {
        slong rows[16];
        slong i;
        for (i = 0; i < n; i++)
            rows[i] = i;

        return _gr_mat_permanent_cofactor_recursive(res, A, rows, 0, ctx);
    }
}

int
gr_mat_permanent_ryser(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong n = A->r;
    slong i, j, jprev, d;
    gr_ptr s, p;

    if (n <= 0 || n > FLINT_BITS - 2)
        return GR_UNABLE;

    status |= gr_zero(res, ctx);

    jprev = 0;

    GR_TMP_INIT_VEC(s, n + 1, ctx);
    p = GR_ENTRY(s, n, ctx->sizeof_elem);

    for (i = 1; i < (WORD(1) << n); i++)
    {
        j = i ^ (i >> 1);

        if (j > jprev)
        {
            d = FLINT_BIT_COUNT(j - jprev) - 1;
            status |= _gr_vec_add(s, s, gr_mat_entry_srcptr(A, d, 0, ctx), n, ctx);
        }
        else
        {
            d = FLINT_BIT_COUNT(jprev - j) - 1;
            status |= _gr_vec_sub(s, s, gr_mat_entry_srcptr(A, d, 0, ctx), n, ctx);
        }

        status |= _gr_vec_product(p, s, n, ctx);

        if ((n + i) & 1)
            status |= gr_sub(res, res, p, ctx);
        else
            status |= gr_add(res, res, p, ctx);

        jprev = j;
    }

    GR_TMP_CLEAR_VEC(s, n + 1, ctx);

    return status;
}

int
gr_mat_permanent_glynn(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong n = A->r;
    slong i, j, jprev, d;
    gr_ptr s, p;
    gr_ptr A2;
    slong sz = ctx->sizeof_elem;

    if (n <= 0 || n > FLINT_BITS - 2)
        return GR_UNABLE;

    status |= gr_zero(res, ctx);

    jprev = 0;

    GR_TMP_INIT_VEC(s, n + 1, ctx);
    GR_TMP_INIT_VEC(A2, n * n, ctx);
    p = GR_ENTRY(s, n, ctx->sizeof_elem);

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            status |= gr_mul_2exp_si(GR_ENTRY(A2, j * n + i, sz), gr_mat_entry_srcptr(A, i, j, ctx), 1, ctx);

    for (i = 0; i < n; i++)
        status |= _gr_vec_sum(GR_ENTRY(s, i, sz), gr_mat_entry_srcptr(A, i, 0, ctx), n, ctx);

    status |= _gr_vec_product(res, s, n, ctx);

    for (i = 1; i < (WORD(1) << (n - 1)); i++)
    {
        j = i ^ (i >> 1);

        if (j > jprev)
        {
            d = FLINT_BIT_COUNT(j - jprev) - 1;
            status |= _gr_vec_sub(s, s, GR_ENTRY(A2, d * n, sz), n, ctx);
        }
        else
        {
            d = FLINT_BIT_COUNT(jprev - j) - 1;
            status |= _gr_vec_add(s, s, GR_ENTRY(A2, d * n, sz), n, ctx);
        }

        status |= _gr_vec_product(p, s, n, ctx);

        if (i & 1)
            status |= gr_sub(res, res, p, ctx);
        else
            status |= gr_add(res, res, p, ctx);

        jprev = j;
    }

    status |= gr_mul_2exp_si(res, res, -(n - 1), ctx);

    GR_TMP_CLEAR_VEC(s, n + 1, ctx);
    GR_TMP_CLEAR_VEC(A2, n * n, ctx);

    return status;
}

/* todo: builtin check for odd characteristic */
static int
can_maybe_divide_by_two(gr_ptr res, gr_ctx_t ctx)
{
    int finite, status;

    finite = gr_ctx_is_finite_characteristic(ctx);

    if (finite == T_TRUE)
    {
        status = gr_set_ui(res, 2, ctx);
        if (status == GR_SUCCESS)
        {
            if (gr_is_invertible(res, ctx) == T_TRUE)
                return 1;
        }

        return 0;
    }

    return 1;
}

int
gr_mat_permanent(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx)
{
    slong n = A->r;

    if (n != A->c)
        return GR_DOMAIN;

    if (n <= 4)
    {
        return gr_mat_permanent_cofactor(res, A, ctx);
    }
    else
    {
        if (can_maybe_divide_by_two(res, ctx))
        {
            if (gr_mat_permanent_glynn(res, A, ctx) == GR_SUCCESS)
                return GR_SUCCESS;
        }

        return gr_mat_permanent_ryser(res, A, ctx);
    }
}

