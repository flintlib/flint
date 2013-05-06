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

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "long_extras.h"
#include "fmpz.h"

int main(void)
{
    int i;
    flint_rand_t state;
    
    printf("init_set_readonly....");
    fflush(stdout);
    
    flint_randinit(state);

    /* Create some small fmpz integers, clear the mpz_t */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f;
        mpz_t z;

        *f = z_randint(state, COEFF_MAX + 1);

        mpz_init(z);
        fmpz_get_mpz(z, f);

        {
            fmpz_t g;

            fmpz_init_set_readonly(g, z);
            fmpz_clear_readonly(g);
        }

        mpz_clear(z);
    }

    /* Create some small fmpz integers, do *not* clear the mpz_t */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f;
        mpz_t z;

        *f = z_randint(state, COEFF_MAX + 1);

        mpz_init(z);
        fmpz_get_mpz(z, f);

        {
            fmpz_t g;

            fmpz_init_set_readonly(g, z);
        }

        mpz_clear(z);
    }

    /* Create some more fmpz integers */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f;
        mpz_t z;

        fmpz_init(f);
        fmpz_randtest(f, state, 2 * FLINT_BITS);
        mpz_init(z);
        fmpz_get_mpz(z, f);

        {
            fmpz_t g, h;

            fmpz_init_set_readonly(g, z);
            fmpz_init(h);
            fmpz_set_mpz(h, z);

            if (!fmpz_equal(g, h))
            {
                printf("FAIL:\n\n");
                printf("g = "), fmpz_print(g), printf("\n");
                printf("h = "), fmpz_print(h), printf("\n");
                gmp_printf("z = %Zd\n", z);
            }

            fmpz_clear_readonly(g);
            fmpz_clear(h);
        }

        fmpz_clear(f);
        mpz_clear(z);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
