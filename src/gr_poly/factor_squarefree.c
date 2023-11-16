/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr_vec.h"
#include "gr_poly.h"

int
gr_poly_factor_squarefree(gr_ptr c, gr_vec_t fac, gr_vec_t exp, const gr_poly_t F, gr_ctx_t ctx)
{
    gr_poly_t f, d, t1;
    gr_poly_t v, w, s;
    slong i;
    int status = GR_SUCCESS;
    fmpz_t e;
    gr_ctx_t poly_ctx, fmpz_ctx;

    /* todo */
    if (gr_ctx_is_finite_characteristic(ctx) != T_FALSE)
        return GR_UNABLE;

    /* for vector of polynomials */
    gr_ctx_init_gr_poly(poly_ctx, ctx);
    /* for vector of exponents */
    gr_ctx_init_fmpz(fmpz_ctx);

    if (F->length == 0)
    {
        status |= gr_zero(c, ctx);
        gr_vec_set_length(fac, 0, poly_ctx);
        gr_vec_set_length(exp, 0, fmpz_ctx);
        gr_ctx_clear(poly_ctx);
        gr_ctx_clear(fmpz_ctx);
        return GR_SUCCESS;
    }

    status |= gr_poly_get_coeff_scalar(c, F, F->length - 1, ctx);

    if (gr_is_zero(c, ctx) != T_FALSE)
    {
        gr_ctx_clear(poly_ctx);
        gr_ctx_clear(fmpz_ctx);
        return GR_UNABLE;
    }

    if (F->length == 1)
    {
        gr_vec_set_length(fac, 0, poly_ctx);
        gr_vec_set_length(exp, 0, fmpz_ctx);
        gr_ctx_clear(poly_ctx);
        gr_ctx_clear(fmpz_ctx);
        return GR_SUCCESS;
    }

    fmpz_init(e);
    gr_poly_init(f, ctx);
    gr_poly_init(d, ctx);
    gr_poly_init(t1, ctx);
    gr_poly_init(v, ctx);
    gr_poly_init(w, ctx);
    gr_poly_init(s, ctx);

    status |= gr_poly_make_monic(f, F, ctx);
    status |= gr_poly_derivative(t1, f, ctx);
    status |= gr_poly_gcd(d, f, t1, ctx);

    if (status != GR_SUCCESS)
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    gr_vec_set_length(fac, 0, ctx);
    gr_vec_set_length(exp, 0, ctx);

    if (d->length == 1)
    {
        status |= gr_vec_append(fac, f, poly_ctx);
        fmpz_one(e);
        status |= gr_vec_append(exp, e, fmpz_ctx);
    }
    else
    {
        /* todo: should be divexact over field */
        status |= gr_poly_divrem(v, s, f, d, ctx);
        status |= gr_poly_divrem(w, s, t1, d, ctx);

        if (status != GR_SUCCESS)
            goto cleanup;

        /* invariant: v is monic, so we don't need to check if it is proper */
        for (i = 1; ; i++)
        {
            status |= gr_poly_derivative(t1, v, ctx);
            status |= gr_poly_sub(s, w, t1, ctx);

            /* check if polynomial is proper */
            if (s->length != 0 && gr_is_zero(gr_poly_entry_ptr(s, s->length - 1, ctx), ctx) != T_FALSE)
            {
                status = GR_UNABLE;
                goto cleanup;
            }

            if (s->length == 0)
            {
                if (v->length > 1)
                {
                    status |= gr_vec_append(fac, v, poly_ctx);
                    fmpz_set_ui(e, i); /* append_ui? */
                    status |= gr_vec_append(exp, e, fmpz_ctx);
                }
                break;
            }

            status |= gr_poly_gcd(d, v, s, ctx);

            if (status != GR_SUCCESS)
                goto cleanup;

            /* should be divexact over field */
            status |= gr_poly_divrem(v, t1, v, d, ctx);
            status |= gr_poly_divrem(w, t1, s, d, ctx);

            if (status != GR_SUCCESS)
                goto cleanup;

            if (d->length > 1)
            {
                status |= gr_vec_append(fac, d, poly_ctx);
                fmpz_set_ui(e, i); /* append_ui? */
                status |= gr_vec_append(exp, e, fmpz_ctx);
            }
        }
    }

cleanup:
    gr_poly_clear(f, ctx);
    gr_poly_clear(d, ctx);
    gr_poly_clear(t1, ctx);
    gr_poly_clear(v, ctx);
    gr_poly_clear(w, ctx);
    gr_poly_clear(s, ctx);
    fmpz_clear(e);

    gr_ctx_clear(poly_ctx);
    gr_ctx_clear(fmpz_ctx);

    return status;
}
