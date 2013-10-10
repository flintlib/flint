/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include "fq_poly.h"
#include "ulong_extras.h"

void
fq_poly_factor_squarefree(fq_poly_factor_t res, const fq_poly_t f,
                          const fq_ctx_t ctx)
{
    fq_poly_t f_d, g, g_1, r;
    fq_t x;
    fmpz_t p;
    slong deg, i, p_ui;

    if (f->length <= 1)
    {
        res->num = 0;
        return;
    }

    if (f->length == 2)
    {
        fq_poly_factor_insert(res, f, 1, ctx);
        return;
    }

    fmpz_init(p);
    fmpz_set(p, fq_ctx_prime(ctx));
    
    deg = fq_poly_degree(f);

    /* Step 1, look at f', if it is zero then we are done since f = h(x)^p
       for some particular h(x), clearly f(x) = sum a_k x^kp, k <= deg(f) */

    fq_init(x);
    fq_poly_init(g_1);
    fq_poly_init(f_d);
    fq_poly_init(g);
    fq_poly_derivative(f_d, f, ctx);

    /* Case 1 */
    if (fq_poly_is_zero(f_d))
    {
        fq_poly_factor_t new_res;
        fq_poly_t h;

        /* We can do this since deg is a multiple of p in this case */
        p_ui = fmpz_get_ui(p);

        fq_poly_init(h);

        for (i = 0; i <= deg / p_ui; i++)   /* this will be an integer since f'=0 */
        {
            fq_poly_get_coeff(x, f, i * p_ui, ctx);
            fq_pth_root(x, x, ctx);
            fq_poly_set_coeff(h, i, x, ctx);
        }

        /* Now run squarefree on h, and return it to the pth power */
        fq_poly_factor_init(new_res, ctx);

        fq_poly_factor_squarefree(new_res, h, ctx);
        fq_poly_factor_pow(new_res, p_ui);

        fq_poly_factor_concat(res, new_res, ctx);
        fq_poly_clear(h);
        fq_poly_factor_clear(new_res);
    }
    else
    {
        fq_poly_t h, z;

        fq_poly_init(r);

        fq_poly_gcd(g, f, f_d, ctx);
        fq_poly_divrem(g_1, r, f, g, ctx);

        i = 1;

        fq_poly_init(h);
        fq_poly_init(z);

        /* Case 2 */
        while (g_1->length > 1)
        {
            fq_poly_gcd(h, g_1, g, ctx);
            fq_poly_divrem(z, r, g_1, h, ctx);

            /* out <- out.z */
            if (z->length > 1)
            {
                fq_poly_factor_insert(res, z, 1, ctx);
                fq_poly_make_monic(res->poly + (res->num - 1),
                                   res->poly + (res->num - 1),
                                   ctx);
                if (res->num)
                    res->exp[res->num - 1] *= i;
            }

            i++;
            fq_poly_set(g_1, h, ctx);
            fq_poly_divrem(g, r, g, h, ctx);
        }

        fq_poly_clear(h);
        fq_poly_clear(z);
        fq_poly_clear(r);

        fq_poly_make_monic(g, g, ctx);

        if (g->length > 1)
        {
            /* so now we multiply res with squarefree(g^1/p) ^ p  */
            fq_poly_t g_p;    /* g^(1/p) */
            fq_poly_factor_t new_res_2;

            fq_poly_init(g_p);
            
            p_ui = fmpz_get_ui(p);

            for (i = 0; i <= fq_poly_degree(g) / p_ui; i++)
            {
                fq_poly_get_coeff(x, g, i * p_ui, ctx);
                fq_pth_root(x, x, ctx);
                fq_poly_set_coeff(g_p, i, x, ctx);
            }

            fq_poly_factor_init(new_res_2, ctx);

            /* squarefree(g^(1/p)) */
            fq_poly_factor_squarefree(new_res_2, g_p, ctx);
            fq_poly_factor_pow(new_res_2, p_ui);

            fq_poly_factor_concat(res, new_res_2, ctx);
            fq_poly_clear(g_p);
            fq_poly_factor_clear(new_res_2);
        }
    }

    fmpz_clear(p);
    fq_clear(x);
    fq_poly_clear(g_1);
    fq_poly_clear(f_d);
    fq_poly_clear(g);
}
