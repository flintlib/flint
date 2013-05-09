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

******************************************************************************/

#include "nmod_poly.h"
#include "ulong_extras.h"

void
nmod_poly_factor_squarefree(nmod_poly_factor_t res, const nmod_poly_t f)
{
    nmod_poly_t f_d, g, g_1;
    mp_limb_t p;
    len_t deg, i;

    if (f->length <= 1) 
    {
        res->num = 0;
        return;
    }

    if (f->length == 2)
    {
        nmod_poly_factor_insert(res, f, 1);
        return;
    }

    p = nmod_poly_modulus(f);
    deg = nmod_poly_degree(f);

    
    /* Step 1, look at f', if it is zero then we are done since f = h(x)^p
       for some particular h(x), clearly f(x) = sum a_k x^kp, k <= deg(f) */

    nmod_poly_init(g_1, p);
    nmod_poly_init(f_d, p);
    nmod_poly_init(g, p);
    nmod_poly_derivative(f_d, f);

    /* Case 1 */
    if (nmod_poly_is_zero(f_d))
    {
        nmod_poly_factor_t new_res;
        nmod_poly_t h;

        nmod_poly_init(h, p);

        for (i = 0; i <= deg / p; i++) /* this will be an integer since f'=0 */
        {
            nmod_poly_set_coeff_ui(h, i, nmod_poly_get_coeff_ui(f, i * p));
        }
        
        /* Now run square-free on h, and return it to the pth power */
        nmod_poly_factor_init(new_res);

        nmod_poly_factor_squarefree(new_res, h);
        nmod_poly_factor_pow(new_res, p);

        nmod_poly_factor_concat(res, new_res);
        nmod_poly_clear(h);
        nmod_poly_factor_clear(new_res);
   }
   else 
   { 
        nmod_poly_t h, z;

        nmod_poly_gcd(g, f, f_d);
        nmod_poly_div(g_1, f, g);

        i = 1;

        nmod_poly_init(h, p);
        nmod_poly_init(z, p);

        /* Case 2 */
        while (!nmod_poly_is_one(g_1)) 
        {
            nmod_poly_gcd(h, g_1, g);
            nmod_poly_div(z, g_1, h);

            /* out <- out.z */
            if (z->length > 1)
            {
                nmod_poly_factor_insert(res, z, 1);
                nmod_poly_make_monic(res->p + (res->num - 1),
                                     res->p + (res->num - 1));
                if (res->num)
                    res->exp[res->num - 1] *= i;
            }

            i++;
            nmod_poly_set(g_1, h);
            nmod_poly_div(g, g, h);
        }

        nmod_poly_clear(h);
        nmod_poly_clear(z);
        
        nmod_poly_make_monic(g, g);

        if (!nmod_poly_is_one(g))
        {
            /* so now we multiply res with square-free(g^1/p) ^ p  */
            nmod_poly_t g_p; /* g^(1/p) */
            nmod_poly_factor_t new_res_2;

            nmod_poly_init(g_p, p);

            for (i = 0; i <= nmod_poly_degree(g) / p; i++)
                nmod_poly_set_coeff_ui(g_p, i, nmod_poly_get_coeff_ui(g, i*p));

            nmod_poly_factor_init(new_res_2);

            /* square-free(g^(1/p)) */
            nmod_poly_factor_squarefree(new_res_2, g_p);
            nmod_poly_factor_pow(new_res_2, p);

            nmod_poly_factor_concat(res, new_res_2);
            nmod_poly_clear(g_p);
            nmod_poly_factor_clear(new_res_2);
        }
   }

    nmod_poly_clear(g_1);
    nmod_poly_clear(f_d);
    nmod_poly_clear(g);
}
