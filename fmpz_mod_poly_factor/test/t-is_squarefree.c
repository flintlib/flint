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

    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly_factor.h"
#include "ulong_extras.h"

int
main(void)
{
    int iter;
    flint_rand_t state;
    flint_randinit(state);

    printf("is_squarefree....");
    fflush(stdout);

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        fmpz_mod_poly_t poly, Q, R, t;
        fmpz_t modulus;
        mp_limb_t mod;
        long i, num_factors, exp, max_exp;
        int v, result;

        mod = n_randtest_prime(state, 0);
        fmpz_init_set_ui(modulus, mod);

        fmpz_mod_poly_init(poly, modulus);
        fmpz_mod_poly_init(t, modulus);
        fmpz_mod_poly_init(Q, modulus);
        fmpz_mod_poly_init(R, modulus);

        fmpz_mod_poly_set_coeff_ui(poly, 0, n_randint(state, mod));
        num_factors = n_randint(state, 5);

        max_exp = 0;
        for (i = 0; i < num_factors; i++)
        {
            do {
                fmpz_mod_poly_randtest(t, state, n_randint(state, 10));
            } while (!fmpz_mod_poly_is_irreducible(t) ||
                    (fmpz_mod_poly_length(t) < 2));

            exp = n_randint(state, 4) + 1;
            if (n_randint(state, 2) == 0)
                exp = 1;

            fmpz_mod_poly_divrem(Q, R, poly, t);
            if (!fmpz_mod_poly_is_zero(R))
            {
                fmpz_mod_poly_pow(t, t, exp);
                fmpz_mod_poly_mul(poly, poly, t);
                max_exp = FLINT_MAX(exp, max_exp);
            }
        }

        v = fmpz_mod_poly_is_squarefree(poly);

        if (v == 1)
            result = (max_exp <= 1 && !fmpz_mod_poly_is_zero(poly));
        else
            result = (max_exp > 1 || fmpz_mod_poly_is_zero(poly));

        if (!result)
        {
            printf("FAIL: ");
            fmpz_print(modulus);
            printf(" %ld, %d\n", max_exp, v);
            fmpz_mod_poly_print(poly); printf("\n");
            abort();
        }

        fmpz_clear(modulus);
        fmpz_mod_poly_clear(poly);
        fmpz_mod_poly_clear(t);
        fmpz_mod_poly_clear(Q);
        fmpz_mod_poly_clear(R);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
