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

    Copyright (C) 2012 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("is_probabprime_BPSW....");
    fflush(stdout);

    /* test primes always pass */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p, F;

        fmpz_init(p);
        fmpz_init(F);

        do {
           fmpz_randbits(p, state, n_randint(state, 330) + 2);
           fmpz_abs(p, p);
        } while (!fmpz_is_probabprime(p) || fmpz_cmp_ui(p, 2) == 0);

        result = fmpz_is_probabprime_BPSW(p);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_print(p); printf("\n");
            abort();
        }

        fmpz_clear(p);
        fmpz_clear(F);
    }

    /* test composites never pass */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p, a, F;

        fmpz_init(p);
        fmpz_init(a);
        fmpz_init(F);

        do {
           fmpz_randbits(p, state, n_randint(state, 100) + 2);
        } while (fmpz_cmp_ui(p, 2) < 0);
        do {
           fmpz_randbits(a, state, n_randint(state, 100) + 2);
        } while (fmpz_cmp_ui(a, 2) < 0);
        
        fmpz_mul(p, p, a);

        result = !fmpz_is_probabprime_BPSW(p);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_print(p); printf("\n");
            abort();
        }

        fmpz_clear(p);
        fmpz_clear(a);
        fmpz_clear(F);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
