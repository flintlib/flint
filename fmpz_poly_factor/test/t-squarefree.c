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
#include "fmpz_poly_factor.h"

int main(void)
{
    int i, result;
    flint_rand_t state;

    printf("squarefree....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 1000; i++)
    {
        fmpz_poly_t f, g[10], h, t;
        fmpz_poly_factor_t fac;
        long k, l, n = n_randint(state, 10) + 1;

        fmpz_poly_init(f);
        fmpz_poly_init(h);
        fmpz_poly_init(t);
        fmpz_poly_factor_init(fac);

        fmpz_poly_one(f);
        for (k = 0; k < n; k++)
        {
            fmpz_poly_init(g[k]);
            fmpz_poly_randtest_not_zero(g[k], state, n_randint(state, 40)+1, 20);
            l = n_randint(state, 2) + 1;
            while (l--)
                fmpz_poly_mul(f, f, g[k]);
        }
        fmpz_poly_factor_squarefree(fac, f);

        /* Squarefree? */
        result = 1;
        for (k = 0; k < fac->num && result; k++)
        {
            fmpz_poly_derivative(h, fac->p + k);
            fmpz_poly_gcd(t, h, fac->p + k);
            result &= fmpz_poly_is_one(t);
        }

        /* Product? */
        fmpz_poly_set_fmpz(h, &(fac->c));
        for (k = 0; k < fac->num; k++)
        {
            if (fac->exp[k] == 1)
                fmpz_poly_mul(h, h, fac->p + k);
            else
            {
                fmpz_poly_pow(t, fac->p + k, fac->exp[k]);
                fmpz_poly_mul(h, h, t);
            }
        }

        result &= fmpz_poly_equal(f, h);
        if (!result)
        {
            printf("FAIL:\n");
            printf("f = "), fmpz_poly_print_pretty(f, "x"), printf("\n\n");
            for (k = 0; k < n; k++)
            {
                printf("g[%ld] = ", k), fmpz_poly_print_pretty(g[k], "x"), printf("\n\n");
            }
            printf("h = "), fmpz_poly_print_pretty(h, "x"), printf("\n\n");
            fmpz_poly_factor_print(fac);
            abort();
        }

        fmpz_poly_clear(f);
        for (k = 0; k < n; k++)
            fmpz_poly_clear(g[k]);
        fmpz_poly_clear(h);
        fmpz_poly_clear(t);
        fmpz_poly_factor_clear(fac);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}

