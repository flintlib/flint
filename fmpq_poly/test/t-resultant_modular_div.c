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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"
#include "ulong_extras.h"

#pragma GCC diagnostic ignored "-Woverlength-strings"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("resultant_div....");
    fflush(stdout);

    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g, h;
        fmpq_t x, y, z, zz;
        fmpz_t den;
        slong nbits;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);
        fmpq_init(zz);
        
        fmpz_init(den);

        fmpq_poly_randtest(f, state, n_randint(state, 60), 60);
        fmpq_poly_randtest(g, state, n_randint(state, 60), 60);
        fmpq_poly_randtest(h, state, n_randint(state, 60), 60);

        fmpz_set(den, fmpq_poly_denref(f));
        fmpq_poly_scalar_mul_fmpz(f, f, den);

        fmpz_set(den, fmpq_poly_denref(g));
        fmpq_poly_scalar_mul_fmpz(g, g, den);

        fmpz_set(den, fmpq_poly_denref(h));
        fmpq_poly_scalar_mul_fmpz(h, h, den);

        flint_printf("\nf:\n");
        fmpq_poly_print_pretty(f, "x");
        flint_printf("\ng:\n");
        fmpq_poly_print_pretty(g, "x");
        flint_printf("\nh:\n");
        fmpq_poly_print_pretty(h, "x");
        flint_printf("\n");


        fmpq_poly_resultant(y, f, g);
        
        flint_printf("y: "), fmpq_print(y), flint_printf("\n");

        if (!fmpz_is_one(fmpq_denref(y)))
        {
            flint_printf("FAIL resultant not integral\n");
            flint_printf("f = "), fmpq_poly_print(f), flint_printf("\n\n");
            flint_printf("g = "), fmpq_poly_print(g), flint_printf("\n\n");
            flint_printf("y = "), fmpq_print(y), flint_printf("\n\n");
            abort();
        }

        fmpq_poly_resultant(z, h, g);
        
        if (!fmpz_is_one(fmpq_denref(z)))
        {
            flint_printf("FAIL resultant not integral\n");
            flint_printf("h = "), fmpq_poly_print(f), flint_printf("\n\n");
            flint_printf("g = "), fmpq_poly_print(g), flint_printf("\n\n");
            flint_printf("z = "), fmpq_print(y), flint_printf("\n\n");
            abort();
        }

        flint_printf("y: "), fmpq_print(y), flint_printf("\n");

        fmpq_mul(y, y, z);

        flint_printf("y: "), fmpq_print(y), flint_printf("\n");
        
        /* y is the correct value we want to compute.
         * Get the number of bits and add some random value. */

        if (!fmpz_is_one(fmpq_denref(y)))
        {
            flint_printf("FAIL resultant not integral\n");
            flint_printf("f = "), fmpq_poly_print(f), flint_printf("\n\n");
            flint_printf("g = "), fmpq_poly_print(g), flint_printf("\n\n");
            flint_printf("y = "), fmpq_print(y), flint_printf("\n\n");
            abort();
        }
     

        nbits = (slong)fmpz_bits(fmpq_numref(y));

        /* Compute res(f h, g) using the divisor res(f, g) and the
         * correct number of bits of the output. */

        fmpq_poly_mul(f, f, h);
        fmpq_poly_resultant(zz, f, g);

        flint_printf("res(f h, g) = "), fmpq_print(zz), flint_printf("\n\n");
        fmpq_poly_resultant_div(x, f, g, fmpq_numref(z), nbits);

        result = fmpq_equal(x, y);
        
        if (!result)
        {
            flint_printf("FAIL (res_modular_div(f, g) == res(f, g):\n");
            flint_printf("f = "), fmpq_poly_print_pretty(f, "x"), flint_printf("\n\n");
            flint_printf("g = "), fmpq_poly_print_pretty(g, "x"), flint_printf("\n\n");
            flint_printf("x = "), fmpq_print(x), flint_printf("\n\n");
            flint_printf("nbits = %wu\n\n");
            flint_printf("div = "), fmpq_print(z), flint_printf("\n\n");
            flint_printf("y = "), fmpq_print(y), flint_printf("\n\n");
            abort();
        }
        
        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
        fmpq_clear(zz);
        
        fmpz_clear(den);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

