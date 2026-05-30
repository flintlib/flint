/*
    Copyright (C) 2023, 2025, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr_vec.h"
#include "gr_poly.h"

int
gr_poly_squarefree_part(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx)
{
    gr_poly_t t;
    int status = GR_SUCCESS;
    truth_t is_fin_char;

    if (poly->length <= 1)
        return gr_poly_one(res, ctx);

    is_fin_char = gr_ctx_is_finite_characteristic(ctx);

    if (is_fin_char == T_UNKNOWN)
        return GR_UNABLE;

    if (is_fin_char == T_TRUE)
    {
        /* The finite characteristic case is more complicated, so just call
           the factorization algorithm and multiply the distinct factors
           together. */
        gr_ptr c;
        gr_vec_t fac, exp;
        gr_ctx_t poly_ctx, fmpz_ctx;
        slong i;

        gr_ctx_init_gr_poly(poly_ctx, ctx);
        gr_ctx_init_fmpz(fmpz_ctx);
        gr_vec_init(fac, 0, poly_ctx);
        gr_vec_init(exp, 0, fmpz_ctx);
        c = gr_heap_init(ctx);

        status |= gr_poly_factor_squarefree(c, fac, exp, poly, ctx);

        if (status == GR_SUCCESS)
        {
            status |= _gr_vec_product(res, fac->entries, fac->length, poly_ctx);
            status |= gr_poly_canonical_associate(res, NULL, res, ctx);
        }

        gr_heap_clear(c, ctx);
        gr_vec_clear(fac, poly_ctx);
        gr_vec_clear(exp, fmpz_ctx);
        gr_ctx_clear(poly_ctx);
        gr_ctx_clear(fmpz_ctx);

        return status;
    }

    /* Characteristic 0: canonical_associate(poly / gcd(poly, poly')). */
    gr_poly_init(t, ctx);
    status |= gr_poly_derivative(t, poly, ctx);
    status |= gr_poly_gcd(t, poly, t, ctx);
    status |= gr_poly_divexact(res, poly, t, ctx);
    status |= gr_poly_canonical_associate(res, NULL, res, ctx);
    gr_poly_clear(t, ctx);

    return status;
}
