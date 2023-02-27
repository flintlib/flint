/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

void
fmpz_mod_poly_factor_squarefree(fmpz_mod_poly_factor_t res,
                             const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_t f_d, g, g_1, r;
    fmpz_t x;
    slong deg, i, p_ui;

    res->num = 0;

    if (f->length <= 1)
        return;

    if (f->length == 2)
    {
        fmpz_mod_poly_factor_insert(res, f, 1, ctx);
        return;
    }

    p_ui = fmpz_get_ui(fmpz_mod_ctx_modulus(ctx)); /* TODO what? */
    deg = fmpz_mod_poly_degree(f, ctx);

    /* Step 1, look at f', if it is zero then we are done since f = h(x)^p
       for some particular h(x), clearly f(x) = sum a_k x^kp, k <= deg(f) */

    fmpz_init(x);
    fmpz_mod_poly_init(g_1, ctx);
    fmpz_mod_poly_init(f_d, ctx);
    fmpz_mod_poly_init(g, ctx);
    fmpz_mod_poly_derivative(f_d, f, ctx);

    /* Case 1 */
    if (fmpz_mod_poly_is_zero(f_d, ctx))
    {
        fmpz_mod_poly_factor_t new_res;
        fmpz_mod_poly_t h;

        fmpz_mod_poly_init(h, ctx);

        for (i = 0; i <= deg / p_ui; i++)   /* this will be an integer since f'=0 */
        {
            fmpz_mod_poly_get_coeff_fmpz(x, f, i * p_ui, ctx);
            fmpz_mod_poly_set_coeff_fmpz(h, i, x, ctx);
        }

        /* Now run squarefree on h, and return it to the pth power */
        fmpz_mod_poly_factor_init(new_res, ctx);

        fmpz_mod_poly_factor_squarefree(new_res, h, ctx);
        fmpz_mod_poly_factor_pow(new_res, p_ui, ctx);

        fmpz_mod_poly_factor_concat(res, new_res, ctx);
        fmpz_mod_poly_clear(h, ctx);
        fmpz_mod_poly_factor_clear(new_res, ctx);
    }
    else
    {
        fmpz_mod_poly_t h, z;

        fmpz_mod_poly_init(r, ctx);

        fmpz_mod_poly_gcd(g, f, f_d, ctx);
        fmpz_mod_poly_divrem(g_1, r, f, g, ctx);

        i = 1;

        fmpz_mod_poly_init(h, ctx);
        fmpz_mod_poly_init(z, ctx);

        /* Case 2 */
        while (g_1->length > 1)
        {
            fmpz_mod_poly_gcd(h, g_1, g, ctx);
            fmpz_mod_poly_divrem(z, r, g_1, h, ctx);

            /* out <- out.z */
            if (z->length > 1)
            {
                fmpz_mod_poly_factor_insert(res, z, 1, ctx);
                fmpz_mod_poly_make_monic(res->poly + (res->num - 1),
                                         res->poly + (res->num - 1), ctx);
                if (res->num)
                    res->exp[res->num - 1] *= i;
            }

            i++;
            fmpz_mod_poly_set(g_1, h, ctx);
            fmpz_mod_poly_divrem(g, r, g, h, ctx);
        }

        fmpz_mod_poly_clear(h, ctx);
        fmpz_mod_poly_clear(z, ctx);
        fmpz_mod_poly_clear(r, ctx);

        fmpz_mod_poly_make_monic(g, g, ctx);

        if (g->length > 1)
        {
            /* so now we multiply res with squarefree(g^1/p) ^ p  */
            fmpz_mod_poly_t g_p;    /* g^(1/p) */
            fmpz_mod_poly_factor_t new_res_2;

            fmpz_mod_poly_init(g_p, ctx);

            for (i = 0; i <= fmpz_mod_poly_degree(g, ctx) / p_ui; i++)
            {
                fmpz_mod_poly_get_coeff_fmpz(x, g, i * p_ui, ctx);
                fmpz_mod_poly_set_coeff_fmpz(g_p, i, x, ctx);
            }

            fmpz_mod_poly_factor_init(new_res_2, ctx);

            /* squarefree(g^(1/p)) */
            fmpz_mod_poly_factor_squarefree(new_res_2, g_p, ctx);
            fmpz_mod_poly_factor_pow(new_res_2, p_ui, ctx);

            fmpz_mod_poly_factor_concat(res, new_res_2, ctx);
            fmpz_mod_poly_clear(g_p, ctx);
            fmpz_mod_poly_factor_clear(new_res_2, ctx);
        }
    }

    fmpz_clear(x);
    fmpz_mod_poly_clear(g_1, ctx);
    fmpz_mod_poly_clear(f_d, ctx);
    fmpz_mod_poly_clear(g, ctx);
}
