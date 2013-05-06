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
    
    printf("mpz_init_set_readonly....");
    fflush(stdout);
    
    flint_randinit(state);

    /* Create some small fmpz integers, clear the mpz_t */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f;
        mpz_t z;

        *f = z_randint(state, COEFF_MAX + 1);

        flint_mpz_init_set_readonly(z, f);
        flint_mpz_clear_readonly(z);
    }

    /* Create some large fmpz integers, do not clear the mpz_t */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f;
        mpz_t z;

        fmpz_init(f);
        fmpz_randtest(f, state, 2 * FLINT_BITS);

        if (COEFF_IS_MPZ(*f))
        {
            flint_mpz_init_set_readonly(z, f);
        }

        fmpz_clear(f);
    }

    /* Create some more fmpz integers */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f, g;
        mpz_t z;

        fmpz_init(f);
        fmpz_init(g);
        fmpz_randtest(f, state, 2 * FLINT_BITS);

        flint_mpz_init_set_readonly(z, f);
        fmpz_set_mpz(g, z);

        if (!fmpz_equal(f, g))
        {
            printf("FAIL:\n\n");
            printf("f = "), fmpz_print(f), printf("\n");
            printf("g = "), fmpz_print(g), printf("\n");
            gmp_printf("z = %Zd\n", z);
        }

        fmpz_clear(f);
        fmpz_clear(g);
        flint_mpz_clear_readonly(z);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
