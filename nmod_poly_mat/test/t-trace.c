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
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    len_t i;
    flint_rand_t state;

    printf("trace....");
    fflush(stdout);

    flint_randinit(state);

    /* Test trace(AB) = trace(BA) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A, B, AB, BA;
        nmod_poly_t trab, trba;
        mp_limb_t mod;
        len_t m, n;

        mod = n_randtest_prime(state, 0);
        m = n_randint(state, 10);
        n = n_randint(state, 10);

        nmod_poly_mat_init(A, m, n, mod);
        nmod_poly_mat_init(B, n, m, mod);
        nmod_poly_mat_init(AB, m, m, mod);
        nmod_poly_mat_init(BA, n, n, mod);

        nmod_poly_init(trab, mod);
        nmod_poly_init(trba, mod);

        nmod_poly_mat_randtest(A, state, 1 + n_randint(state, 10));
        nmod_poly_mat_randtest(B, state, 1 + n_randint(state, 10));

        nmod_poly_mat_mul(AB, A, B);
        nmod_poly_mat_mul(BA, B, A);

        nmod_poly_mat_trace(trab, AB);
        nmod_poly_mat_trace(trba, BA);

        if (!nmod_poly_equal(trab, trba))
        {
            printf("FAIL:\n");
            nmod_poly_mat_print(A, "x"), printf("\n");
            nmod_poly_mat_print(B, "x"), printf("\n");
            nmod_poly_mat_print(AB, "x"), printf("\n");
            nmod_poly_mat_print(BA, "x"), printf("\n");
            printf("tr(AB): "),  nmod_poly_print(trab),    printf("\n");
            printf("tr(BA): "),  nmod_poly_print(trba),    printf("\n");
            abort();
        }

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(AB);
        nmod_poly_mat_clear(BA);
        nmod_poly_clear(trab);
        nmod_poly_clear(trba);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
