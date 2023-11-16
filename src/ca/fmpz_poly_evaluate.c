/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_fmpz_poly_evaluate(ca_t res, const fmpz_poly_t poly, const ca_t x, ca_ctx_t ctx)
{
    if (fmpz_poly_is_zero(poly))
    {
        ca_zero(res, ctx);
    }
    else if (fmpz_poly_length(poly) == 1)
    {
        ca_set_fmpz(res, poly->coeffs + 0, ctx);
    }
    else if (CA_IS_QQ(x, ctx))
    {
        fmpq_t t;
        fmpq_init(t);
        fmpz_poly_evaluate_fmpq(t, poly, CA_FMPQ(x));
        ca_set_fmpq(res, t, ctx);
        fmpq_clear(t);
    }
    /* todo: fast modular composition for number field elements */
    /* else if (CA_FIELD_IS_NF(CA_FIELD(x, ctx))) */
    else
    {
        ca_t t;  /* for aliasing */
        slong i, n;

        ca_init(t, ctx);

        /* todo: rectangular splitting? */
        n = fmpz_poly_degree(poly);
        ca_set_fmpz(t, poly->coeffs + n, ctx);

        for (i = n - 1; i >= 0; i--)
        {
            ca_mul(t, t, x, ctx);
            ca_add_fmpz(t, t, poly->coeffs + i, ctx);
        }

        ca_clear(t, ctx);
    }
}
