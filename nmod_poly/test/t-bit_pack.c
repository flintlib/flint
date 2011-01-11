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

    printf("bit_pack/bit_unpack....");
    fflush(stdout);

    /* Check aliasing of a and c */
    for (i = 0; i < 10000; i++)
    {
        nmod_poly_t a, b;
        mp_limb_t n;
        ulong bits;
        mp_ptr mpn;

        do
        {
            n = n_randtest_not_zero(state);
        } while (n == 1);
        bits = 2 * FLINT_BIT_COUNT(n) + n_randint(state, FLINT_BITS);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        do
        {
            nmod_poly_randtest(a, state, n_randint(state, 100));
        } while (a->length == 0);

        mpn =
            malloc(sizeof(mp_limb_t) *
                   ((bits * a->length - 1) / FLINT_BITS + 1));

        _nmod_poly_bit_pack(mpn, a->coeffs, a->length, bits);
        nmod_poly_fit_length(b, a->length);
        _nmod_poly_bit_unpack(b->coeffs, mpn, a->length, bits, a->mod);
        b->length = a->length;

        result = (nmod_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
