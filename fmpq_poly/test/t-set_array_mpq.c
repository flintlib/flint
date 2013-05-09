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
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    ulong cflags = 0UL;

    printf("set_array_mpq....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        len_t j, n = 100;
        fmpq_poly_t f, g;
        mpq_t * a;

        a = (mpq_t *) flint_malloc(n * sizeof(mpq_t));
        for (j = 0; j < n; j++)
            mpq_init(a[j]);

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_randtest(f, state, n_randint(state, n), 200);
        for (j = 0; j < f->length; j++)
            fmpq_poly_get_coeff_mpq(a[j], f, j);

        fmpq_poly_set_array_mpq(g, (const mpq_t *) a, n);

        cflags |= fmpq_poly_is_canonical(g) ? 0 : 1;
        result = (fmpq_poly_equal(f, g) && !cflags);
        if (!result)
        {
            printf("FAIL:\n");
            printf("f = "), fmpq_poly_debug(f), printf("\n\n");
            printf("g = "), fmpq_poly_debug(g), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        for (j = 0; j < n; j++)
            mpq_clear(a[j]);
        flint_free(a);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
