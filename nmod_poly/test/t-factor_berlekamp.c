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

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int iter;
    flint_rand_t state;
    flint_randinit(state);

    printf("factor_berlekamp....");
    fflush(stdout);

    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        int result = 1;
        nmod_poly_t pol1, poly, quot, rem;
        nmod_poly_factor_t res;
        mp_limb_t modulus;
        len_t i, length, num;

        modulus = n_randtest_prime(state, 0);

        nmod_poly_init(pol1, modulus);
        nmod_poly_init(poly, modulus);
        nmod_poly_init(quot, modulus);
        nmod_poly_init(rem, modulus);
     
        length = n_randint(state, 10) + 2;
        do 
        {
            nmod_poly_randtest(pol1, state, length);
            if (pol1->length)
                nmod_poly_make_monic(pol1, pol1);
        }
        while ((!nmod_poly_is_irreducible(pol1)) || (pol1->length < 2));

        num = n_randint(state, 5) + 1;
        for (i = 1; i < num; i++)
        {
            do 
            {
                length = n_randint(state, 10) + 2;
                nmod_poly_randtest(poly, state, length); 
                if (poly->length)
                {
                    nmod_poly_make_monic(poly, poly);
                    nmod_poly_divrem(quot, rem, pol1, poly);
                }
            }
            while ((!nmod_poly_is_irreducible(poly)) || (poly->length < 2)
                || (rem->length == 0));
            nmod_poly_mul(pol1, pol1, poly);
        }
     
        nmod_poly_factor_init(res);
        nmod_poly_factor_berlekamp(res, pol1);

        result = (res->num == num);
        if (!result)
        {
            printf("FAIL: %lu, %ld, %ld\n", modulus, num, res->num);
            abort();
        }
      
        nmod_poly_clear(quot);
        nmod_poly_clear(rem);
        nmod_poly_clear(pol1);
        nmod_poly_clear(poly);
        nmod_poly_factor_clear(res);
    }

    flint_randclear(state);
    printf("PASS\n");
    return 0;
}
