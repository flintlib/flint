/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "acb_poly.h"

static acb_poly_struct*
_acb_poly_vec_init(slong n)
{
    slong k;
    acb_poly_struct* ptr = flint_malloc(n * sizeof(acb_poly_struct));
    for (k = 0; k < n; k++)
    {
        acb_poly_init(&ptr[k]);
    }
    return ptr;
}

static void
_acb_poly_vec_clear(acb_poly_struct* ptr, slong n)
{
    slong k;
    for (k = 0; k < n; k++)
    {
        acb_poly_clear(&ptr[k]);
    }
    flint_free(ptr);
}

void
acb_poly_eval_fmpz_mpoly(acb_poly_t res, const fmpz_mpoly_t pol,
    const acb_poly_struct* val, const fmpz_mpoly_ctx_t ctx, slong prec)
{
    slong n = fmpz_mpoly_ctx_nvars(ctx);
    slong L = fmpz_mpoly_length(pol, ctx);
    slong* degrees = flint_malloc(n * sizeof(slong));
    slong j, k;
    acb_poly_struct** powers = flint_malloc(n * sizeof(acb_struct*));
    acb_poly_t ev, temp;
    fmpz_poly_t c;
    fmpz_t coeff;
    slong exp;

    fmpz_mpoly_degrees_si(degrees, pol, ctx);
    for (k = 0; k < n; k++)
    {
        powers[k] = _acb_poly_vec_init(degrees[k] + 2);
        acb_poly_one(&(powers[k][0]));
        for (j = 1; j <= degrees[k]; j++)
        {
            acb_poly_mul(&(powers[k][j]), &(powers[k][j - 1]), &val[k], prec);
        }
    }
    acb_poly_init(ev);
    acb_poly_init(temp);
    fmpz_init(coeff);
    fmpz_poly_init(c);

    acb_poly_zero(ev);
    for (j = 0; j < L; j++)
    {
        fmpz_mpoly_get_term_coeff_fmpz(coeff, pol, j, ctx);
        fmpz_poly_set_fmpz(c, coeff);
        acb_poly_set_fmpz_poly(temp, c, prec);
        for (k = 0; k < n; k++)
        {
            exp = fmpz_mpoly_get_term_var_exp_si(pol, j, k, ctx);
            acb_poly_mul(temp, temp, &(powers[k][exp]), prec);
        }
        acb_poly_add(ev, ev, temp, prec);
    }

    acb_poly_set(res, ev);

    acb_poly_clear(ev);
    acb_poly_clear(temp);
    fmpz_clear(coeff);
    fmpz_poly_clear(c);
    for (k = 0; k < n; k++)
    {
        _acb_poly_vec_clear(powers[k], degrees[k]+2);
    }
    flint_free(degrees);
    flint_free(powers);
}
