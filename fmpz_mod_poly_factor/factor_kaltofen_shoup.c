/*
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "fmpz_mod_poly.h"
#include "thread_support.h"

void
fmpz_mod_poly_factor_kaltofen_shoup(fmpz_mod_poly_factor_t res,
                          const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
{
    slong i, j, k;
    slong num_threads = flint_get_num_threads();
    fmpz_mod_poly_t t, DDxp, EDxp;
    fmpz_mod_poly_factor_t SF, DD, ED;

    res->num = 0;

    fmpz_mod_poly_init(t, ctx);
    fmpz_mod_poly_make_monic(t, poly, ctx);
    
    if (poly->length <= 2)
    {
        fmpz_mod_poly_factor_insert(res, t, 1, ctx);
        fmpz_mod_poly_clear(t, ctx);
        return;        
    }

    fmpz_mod_poly_init(DDxp, ctx);
    fmpz_mod_poly_init(EDxp, ctx);
    fmpz_mod_poly_factor_init(SF, ctx);
    fmpz_mod_poly_factor_init(DD, ctx);
    fmpz_mod_poly_factor_init(ED, ctx);

    fmpz_mod_poly_factor_squarefree(SF, t, ctx);
    for (i = 0; i < SF->num; i++)
    {
        fmpz_mod_poly_struct * f = SF->poly + i;

        fmpz_mod_poly_reverse(t, f, f->length, ctx);
        fmpz_mod_poly_inv_series_newton(t, t, f->length, ctx);
        fmpz_mod_poly_powmod_x_fmpz_preinv(DDxp, fmpz_mod_ctx_modulus(ctx), f, t, ctx);

        /*
            TODO Since there is a fair amount of code duplicated in the
            threaded version, the two functions should be combined into one.
        */
        if (num_threads > 1 && f->length > 256*num_threads)
            fmpz_mod_poly_factor_distinct_deg_threaded_with_frob(DD, f, t, DDxp, ctx);
        else
            fmpz_mod_poly_factor_distinct_deg_with_frob(DD, f, t, DDxp, ctx);

        for (j = 0; j < DD->num; j++)
        {
            fmpz_mod_poly_divrem(t, EDxp, DDxp, DD->poly + j, ctx);
            fmpz_mod_poly_factor_equal_deg_with_frob(ED, DD->poly + j,
                                                        DD->exp[j], EDxp, ctx);
            fmpz_mod_poly_factor_fit_length(res, res->num + ED->num, ctx);
            for (k = 0; k < ED->num; k++)
            {
                fmpz_mod_poly_swap(res->poly + res->num, ED->poly + k, ctx);
                res->exp[res->num] = SF->exp[i];
                res->num++;
            }
        }
    }

    fmpz_mod_poly_clear(t, ctx);
    fmpz_mod_poly_clear(DDxp, ctx);
    fmpz_mod_poly_clear(EDxp, ctx);
    fmpz_mod_poly_factor_clear(SF, ctx);
    fmpz_mod_poly_factor_clear(DD, ctx);
    fmpz_mod_poly_factor_clear(ED, ctx);
}
