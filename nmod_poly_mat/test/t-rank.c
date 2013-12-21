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

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "nmod_mat.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"
#include "fmpz.h"

int
main(void)
{
    slong i;

    FLINT_TEST_INIT(state);

    flint_printf("rank....");
    fflush(stdout);    

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A;
        nmod_mat_t Ax;
        mp_limb_t mod, x;
        slong j, m, n, deg, rank, zrank;
        float density;

        /* Don't pick a too small modulus, to avoid failure in
           the probabilistic rank computation (todo: test
           for small moduli) */
        do {
            mod = n_randtest_prime(state, 0);
        } while (mod < 20);

        m = n_randint(state, 15);
        n = n_randint(state, 15);
        deg = 1 + n_randint(state, 5);
        density = n_randint(state, 100) * 0.01;

        nmod_poly_mat_init(A, m, n, mod);
        nmod_mat_init(Ax, m, n, mod);

        nmod_poly_mat_randtest_sparse(A, state, deg, density);

        /* Probabilistic rank computation */
        zrank = 0;
        for (j = 0; j < 5; j++)
        {
            slong r;
            x = n_randint(state, mod);
            nmod_poly_mat_evaluate_nmod(Ax, A, x);
            r = nmod_mat_rank(Ax);
            zrank = FLINT_MAX(zrank, r);
        }

        rank = nmod_poly_mat_rank(A);

        if (rank != zrank)
        {
            flint_printf("FAIL:\n");
            flint_printf("A:\n");
            nmod_poly_mat_print(A, "x");
            flint_printf("Computed rank: %wd (zrank = %wd)\n", rank, zrank);
            abort();
        }

        nmod_mat_clear(Ax);
        nmod_poly_mat_clear(A);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
