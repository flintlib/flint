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
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("evaluate_mod....");
    fflush(stdout);

    flint_randinit(state);

    /* Compare with evaluation over the integers */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t b, s;
        fmpz_poly_t f;
        mp_limb_t a, n, r;

        fmpz_poly_init(f);
        fmpz_poly_randtest(f, state, n_randint(state, 10), 20);

        n = n_randtest_not_zero(state);
        a = n_randint(state, n);

        fmpz_init(b);
        fmpz_init(s);
        fmpz_set_ui(b, a);

        r = fmpz_poly_evaluate_mod(f, a, n);
        fmpz_poly_evaluate_fmpz(s, f, b);

        result = (r == fmpz_mod_ui(s, s, n));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_print(f), printf("\n\n");
            gmp_printf("a = %Mu\n\n", a);
            gmp_printf("n = %Mu\n\n", n);
            gmp_printf("r = %Mu\n\n", r);
            printf("s = "), fmpz_print(s), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(f);
        fmpz_clear(b);
        fmpz_clear(s);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
