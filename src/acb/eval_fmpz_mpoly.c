/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_eval_fmpz_mpoly(acb_t res, const fmpz_mpoly_t pol, acb_srcptr val,
    const fmpz_mpoly_ctx_t ctx, slong prec)
{
    slong n = fmpz_mpoly_ctx_nvars(ctx);
    slong L = fmpz_mpoly_length(pol, ctx);
    slong* degrees = flint_malloc(n * sizeof(slong));
    slong j, k;
    acb_ptr* powers = flint_malloc(n * sizeof(acb_ptr));
    acb_t ev, temp;
    fmpz_t coeff;
    slong exp;

    fmpz_mpoly_degrees_si(degrees, pol, ctx);
    for (k = 0; k < n; k++)
    {
        powers[k] = _acb_vec_init(degrees[k] + 2);
        acb_one(&(powers[k][0]));
        for (j = 1; j <= degrees[k]; j++)
        {
            acb_mul(&(powers[k][j]), &(powers[k][j - 1]), &val[k], prec);
        }
    }
    acb_init(ev);
    acb_init(temp);
    fmpz_init(coeff);

    acb_zero(ev);
    for (j = 0; j < L; j++)
    {
        fmpz_mpoly_get_term_coeff_fmpz(coeff, pol, j, ctx);
        acb_set_fmpz(temp, coeff);
        for (k = 0; k < n; k++)
        {
            exp = fmpz_mpoly_get_term_var_exp_si(pol, j, k, ctx);
            acb_mul(temp, temp, &(powers[k][exp]), prec);
        }
        acb_add(ev, ev, temp, prec);
    }

    acb_set(res, ev);

    acb_clear(ev);
    acb_clear(temp);
    fmpz_clear(coeff);
    for (k = 0; k < n; k++)
    {
        _acb_vec_clear(powers[k], degrees[k]+2);
    }
    flint_free(degrees);
    flint_free(powers);
}
