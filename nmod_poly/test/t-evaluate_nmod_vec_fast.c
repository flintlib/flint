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

    Copyright (C) 2010 William Hart
    Copyright (C) 2011, 2012 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result = 1;
    flint_rand_t state;
    flint_randinit(state);
    
    printf("evaluate_nmod_vec_fast....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t P, Q;
        mp_ptr x, y, z;
        mp_limb_t mod;
        long j, n, npoints;

        mod = n_randtest_prime(state, 0);
        npoints = n_randint(state, 100);
        n = n_randint(state, 100);

        nmod_poly_init(P, mod);
        nmod_poly_init(Q, mod);
        x = _nmod_vec_init(npoints);
        y = _nmod_vec_init(npoints);
        z = _nmod_vec_init(npoints);

        nmod_poly_randtest(P, state, n);

        for (j = 0; j < npoints; j++)
            x[j] = n_randint(state, mod);

        nmod_poly_evaluate_nmod_vec_iter(y, P, x, npoints);
        nmod_poly_evaluate_nmod_vec_fast(z, P, x, npoints);

        result = _nmod_vec_equal(y, z, npoints);

        if (!result)
        {
            printf("FAIL:\n");
            printf("mod=%lu, n=%ld, npoints=%ld\n\n", mod, n, npoints);
            printf("P: "); nmod_poly_print(P); printf("\n\n");
            abort();
        }

        nmod_poly_clear(P);
        nmod_poly_clear(Q);
        _nmod_vec_clear(x);
        _nmod_vec_clear(y);
        _nmod_vec_clear(z);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
