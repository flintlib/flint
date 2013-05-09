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

    Copyright (C) 2010, 2012 William Hart
    Copyright (C) 2011, 2012 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result = 1;
    flint_rand_t state;
    flint_randinit(state);
    
    printf("evaluate_fmpz_vec_fast....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t P;
        fmpz * x, * y, * z;
        fmpz_t mod;
        len_t j, n, npoints;

        fmpz_init(mod);
        
        do 
        {
           fmpz_randtest_unsigned(mod, state, 5);
           fmpz_add_ui(mod, mod, 2);
        } while (!fmpz_is_probabprime(mod));
        
        npoints = n_randint(state, 10);
        n = n_randint(state, 10);

        fmpz_mod_poly_init(P, mod);
        x = _fmpz_vec_init(npoints);
        y = _fmpz_vec_init(npoints);
        z = _fmpz_vec_init(npoints);

        fmpz_mod_poly_randtest(P, state, n);

        for (j = 0; j < npoints; j++)
            fmpz_randtest_mod(x + j, state, mod);

        fmpz_mod_poly_evaluate_fmpz_vec_iter(y, P, x, npoints);
        fmpz_mod_poly_evaluate_fmpz_vec_fast(z, P, x, npoints);
        
        result = _fmpz_vec_equal(y, z, npoints);

        if (!result)
        {
            printf("FAIL:\n");
            printf("mod=");
            fmpz_print(mod);
            printf(", n=%ld, npoints=%ld\n\n", n, npoints);
            printf("P: "); fmpz_mod_poly_print(P); printf("\n\n");
            for (j = 0; j < npoints; j++)
               fmpz_print(y + j), printf(" ");
            printf("\n");
            for (j = 0; j < npoints; j++)
               fmpz_print(z + j), printf(" ");
            printf("\n");
            abort();
        }

        fmpz_clear(mod);
        fmpz_mod_poly_clear(P);
        _fmpz_vec_clear(x, npoints);
        _fmpz_vec_clear(y, npoints);
        _fmpz_vec_clear(z, npoints);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
