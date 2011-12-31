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
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    flint_randinit(state);

    printf("rem_q1....");
    fflush(stdout);

    /* Check result of divrem */
    for (i = 0; i < 5000; i++)
    {
        nmod_poly_t a, b, q0, r0, r;

        mp_limb_t n = n_randprime(state, n_randint(state,FLINT_BITS-1)+2, 0);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(q0, n);
        nmod_poly_init(r0, n);
        nmod_poly_init(r, n);
        do nmod_poly_randtest(a, state, n_randint(state, 1000));
        while (a->length < 2);
        nmod_poly_fit_length(b, a->length - 1);
        mpn_zero(b->coeffs, a->length - 1);
        nmod_poly_randtest_not_zero(b, state, n_randint(state, 1000) + 1);
        do b->coeffs[a->length - 2] = n_randint(state, n);
        while (b->coeffs[a->length - 2] == 0);
        b->length = a->length - 1;

        nmod_poly_divrem_q1(q0, r0, a, b);
        nmod_poly_rem_q1(r, a, b);

        result = (nmod_poly_equal(r0, r));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(q0), printf("\n\n");
            nmod_poly_print(r0), printf("\n\n");
            nmod_poly_print(r), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }
        
        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(q0);
        nmod_poly_clear(r0);
        nmod_poly_clear(r);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
