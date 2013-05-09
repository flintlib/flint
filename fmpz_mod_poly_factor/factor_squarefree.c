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

    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

******************************************************************************/

#include "fmpz_mod_poly_factor.h"
#include "ulong_extras.h"

void
fmpz_mod_poly_factor_squarefree(fmpz_mod_poly_factor_t res,
                                const fmpz_mod_poly_t f)
{
    fmpz_mod_poly_t f_d, g, g_1, r;
    fmpz_t p, x;
    len_t deg, i, p_ui;

    if (f->length <= 1)
    {
        res->num = 0;
        return;
    }

    if (f->length == 2)
    {
        fmpz_mod_poly_factor_insert(res, f, 1);
        return;
    }

    fmpz_init(p);
    fmpz_set(p, &f->p);
    p_ui = fmpz_get_ui(p);
    deg = fmpz_mod_poly_degree(f);

    /* Step 1, look at f', if it is zero then we are done since f = h(x)^p
       for some particular h(x), clearly f(x) = sum a_k x^kp, k <= deg(f) */

    fmpz_init(x);
    fmpz_mod_poly_init(g_1, p);
    fmpz_mod_poly_init(f_d, p);
    fmpz_mod_poly_init(g, p);
    fmpz_mod_poly_derivative(f_d, f);

    /* Case 1 */
    if (fmpz_mod_poly_is_zero(f_d))
    {
        fmpz_mod_poly_factor_t new_res;
        fmpz_mod_poly_t h;

        fmpz_mod_poly_init(h, p);

        for (i = 0; i <= deg / p_ui; i++)   /* this will be an integer since f'=0 */
        {
            fmpz_mod_poly_get_coeff_fmpz(x, f, i * p_ui);
            fmpz_mod_poly_set_coeff_fmpz(h, i, x);
        }

        /* Now run squarefree on h, and return it to the pth power */
        fmpz_mod_poly_factor_init(new_res);

        fmpz_mod_poly_factor_squarefree(new_res, h);
        fmpz_mod_poly_factor_pow(new_res, p_ui);

        fmpz_mod_poly_factor_concat(res, new_res);
        fmpz_mod_poly_clear(h);
        fmpz_mod_poly_factor_clear(new_res);
    }
    else
    {
        fmpz_mod_poly_t h, z;

        fmpz_mod_poly_init(r, p);

        fmpz_mod_poly_gcd(g, f, f_d);
        fmpz_mod_poly_divrem(g_1, r, f, g);

        i = 1;

        fmpz_mod_poly_init(h, p);
        fmpz_mod_poly_init(z, p);

        /* Case 2 */
        while (g_1->length > 1)
        {
            fmpz_mod_poly_gcd(h, g_1, g);
            fmpz_mod_poly_divrem(z, r, g_1, h);

            /* out <- out.z */
            if (z->length > 1)
            {
                fmpz_mod_poly_factor_insert(res, z, 1);
                fmpz_mod_poly_make_monic(res->poly + (res->num - 1),
                                         res->poly + (res->num - 1));
                if (res->num)
                    res->exp[res->num - 1] *= i;
            }

            i++;
            fmpz_mod_poly_set(g_1, h);
            fmpz_mod_poly_divrem(g, r, g, h);
        }

        fmpz_mod_poly_clear(h);
        fmpz_mod_poly_clear(z);
        fmpz_mod_poly_clear(r);

        fmpz_mod_poly_make_monic(g, g);

        if (g->length > 1)
        {
            /* so now we multiply res with squarefree(g^1/p) ^ p  */
            fmpz_mod_poly_t g_p;    /* g^(1/p) */
            fmpz_mod_poly_factor_t new_res_2;

            fmpz_mod_poly_init(g_p, p);

            for (i = 0; i <= fmpz_mod_poly_degree(g) / p_ui; i++)
            {
                fmpz_mod_poly_get_coeff_fmpz(x, g, i * p_ui);
                fmpz_mod_poly_set_coeff_fmpz(g_p, i, x);
            }

            fmpz_mod_poly_factor_init(new_res_2);

            /* squarefree(g^(1/p)) */
            fmpz_mod_poly_factor_squarefree(new_res_2, g_p);
            fmpz_mod_poly_factor_pow(new_res_2, p_ui);

            fmpz_mod_poly_factor_concat(res, new_res_2);
            fmpz_mod_poly_clear(g_p);
            fmpz_mod_poly_factor_clear(new_res_2);
        }
    }

    fmpz_clear(p);
    fmpz_clear(x);
    fmpz_mod_poly_clear(g_1);
    fmpz_mod_poly_clear(f_d);
    fmpz_mod_poly_clear(g);
}
