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

    printf("2norm_normalised_bits....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c;
        fmpz_poly_t f;
        mp_bitcnt_t b1, b2;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_poly_init(f);
        do {
            fmpz_poly_randtest(f, state, n_randint(state, 100) + 1, 200);
        } while (f->length == 0);

        fmpz_poly_2norm(a, f);
        fmpz_abs(b, fmpz_poly_lead(f));
        fmpz_fdiv_q(c, a, b);
        b1 = fmpz_bits(c);

        b2 = _fmpz_poly_2norm_normalised_bits(f->coeffs, f->length);

        result = (b1 == b2 || b1 + 1 == b2);
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_print(a), printf("\n\n");
            fmpz_print(b), printf("\n\n");
            fmpz_print(c), printf("\n\n");
            printf("b1 = %ld, b2 = %ld\n", b1, b2);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_poly_clear(f);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
