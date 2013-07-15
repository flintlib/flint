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

    Copyright (C) 2009 William Hart

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
    flint_rand_t state;

    printf("jacobi....");
    fflush(stdout);

    flint_randinit(state);
    _flint_rand_init_gmp(state);

    for (i = 0; i < 3000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, p;
        mpz_t b, q;
        int r1, r2;

        fmpz_init(a);
        fmpz_init(p);

        mpz_init(b);
        mpz_init(q);

        mpz_rrandomb(q, state->gmp_state, n_randint(state, 200) + 1);
#ifdef mpz_next_likely_prime
        mpz_next_likely_prime(q, q, state->gmp_state);
#else
        mpz_nextprime(q, q);
#endif
        fmpz_set_mpz(p, q);

        mpz_rrandomb(b, state->gmp_state, n_randint(state, 200) + 1);
        mpz_mod(b, b, q);
        if (n_randint(state, 2))
            mpz_neg(b, b);
        fmpz_set_mpz(a, b);

        r1 = fmpz_jacobi(a, p);
        r2 = mpz_jacobi(b, q);
        result = (r1 == r2);

        if (!result)
        {
            printf("FAIL:\n");
            gmp_printf("b = %Zd, q = %Zd\n", b, q);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(p);

        mpz_clear(b);
        mpz_clear(q);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
