/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include "ulong_extras.h"

void
TEMPLATE(T, poly_factor_squarefree) (TEMPLATE(T, poly_factor_t) res,
                                     const TEMPLATE(T, poly_t) f,
                                     const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_t) f_d, g, g_1, r;
    TEMPLATE(T, t) x;
    fmpz_t p;
    slong deg, i, p_ui;

    if (f->length <= 1)
    {
        res->num = 0;
        return;
    }

    if (f->length == 2)
    {
        TEMPLATE(T, poly_factor_insert) (res, f, 1, ctx);
        return;
    }

    fmpz_init(p);
    fmpz_set(p, TEMPLATE(T, ctx_prime) (ctx));

    deg = TEMPLATE(T, poly_degree) (f, ctx);

    /* Step 1, look at f', if it is zero then we are done since f = h(x)^p
       for some particular h(x), clearly f(x) = sum a_k x^kp, k <= deg(f) */

    TEMPLATE(T, init) (x, ctx);
    TEMPLATE(T, poly_init) (g_1, ctx);
    TEMPLATE(T, poly_init) (f_d, ctx);
    TEMPLATE(T, poly_init) (g, ctx);
    TEMPLATE(T, poly_derivative) (f_d, f, ctx);

    /* Case 1 */
    if (TEMPLATE(T, poly_is_zero) (f_d, ctx))
    {
        TEMPLATE(T, poly_factor_t) new_res;
        TEMPLATE(T, poly_t) h;

        /* We can do this since deg is a multiple of p in this case */
        p_ui = fmpz_get_ui(p);

        TEMPLATE(T, poly_init) (h, ctx);

        for (i = 0; i <= deg / p_ui; i++)   /* this will be an integer since f'=0 */
        {
            TEMPLATE(T, poly_get_coeff) (x, f, i * p_ui, ctx);
            TEMPLATE(T, pth_root) (x, x, ctx);
            TEMPLATE(T, poly_set_coeff) (h, i, x, ctx);
        }

        /* Now run squarefree on h, and return it to the pth power */
        TEMPLATE(T, poly_factor_init) (new_res, ctx);

        TEMPLATE(T, poly_factor_squarefree) (new_res, h, ctx);
        TEMPLATE(T, poly_factor_pow) (new_res, p_ui, ctx);

        TEMPLATE(T, poly_factor_concat) (res, new_res, ctx);
        TEMPLATE(T, poly_clear) (h, ctx);
        TEMPLATE(T, poly_factor_clear) (new_res, ctx);
    }
    else
    {
        TEMPLATE(T, poly_t) h, z;

        TEMPLATE(T, poly_init) (r, ctx);

        TEMPLATE(T, poly_gcd) (g, f, f_d, ctx);
        TEMPLATE(T, poly_divrem) (g_1, r, f, g, ctx);

        i = 1;

        TEMPLATE(T, poly_init) (h, ctx);
        TEMPLATE(T, poly_init) (z, ctx);

        /* Case 2 */
        while (g_1->length > 1)
        {
            TEMPLATE(T, poly_gcd) (h, g_1, g, ctx);
            TEMPLATE(T, poly_divrem) (z, r, g_1, h, ctx);

            /* out <- out.z */
            if (z->length > 1)
            {
                TEMPLATE(T, poly_factor_insert) (res, z, 1, ctx);
                TEMPLATE(T, poly_make_monic) (res->poly + (res->num - 1),
                                              res->poly + (res->num - 1), ctx);
                if (res->num)
                    res->exp[res->num - 1] *= i;
            }

            i++;
            TEMPLATE(T, poly_set) (g_1, h, ctx);
            TEMPLATE(T, poly_divrem) (g, r, g, h, ctx);
        }

        TEMPLATE(T, poly_clear) (h, ctx);
        TEMPLATE(T, poly_clear) (z, ctx);
        TEMPLATE(T, poly_clear) (r, ctx);

        TEMPLATE(T, poly_make_monic) (g, g, ctx);

        if (g->length > 1)
        {
            /* so now we multiply res with squarefree(g^1/p) ^ p  */
            TEMPLATE(T, poly_t) g_p;    /* g^(1/p) */
            TEMPLATE(T, poly_factor_t) new_res_2;

            TEMPLATE(T, poly_init) (g_p, ctx);

            p_ui = fmpz_get_ui(p);

            for (i = 0; i <= TEMPLATE(T, poly_degree) (g, ctx) / p_ui; i++)
            {
                TEMPLATE(T, poly_get_coeff) (x, g, i * p_ui, ctx);
                TEMPLATE(T, pth_root) (x, x, ctx);
                TEMPLATE(T, poly_set_coeff) (g_p, i, x, ctx);
            }

            TEMPLATE(T, poly_factor_init) (new_res_2, ctx);

            /* squarefree(g^(1/p)) */
            TEMPLATE(T, poly_factor_squarefree) (new_res_2, g_p, ctx);
            TEMPLATE(T, poly_factor_pow) (new_res_2, p_ui, ctx);

            TEMPLATE(T, poly_factor_concat) (res, new_res_2, ctx);
            TEMPLATE(T, poly_clear) (g_p, ctx);
            TEMPLATE(T, poly_factor_clear) (new_res_2, ctx);
        }
    }

    fmpz_clear(p);
    TEMPLATE(T, clear) (x, ctx);
    TEMPLATE(T, poly_clear) (g_1, ctx);
    TEMPLATE(T, poly_clear) (f_d, ctx);
    TEMPLATE(T, poly_clear) (g, ctx);
}


#endif
