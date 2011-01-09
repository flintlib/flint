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
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    ulong cflags = 0UL;
    fmpz_randstate_t state;

    printf("rescale....");
    fflush(stdout);

    fmpq_poly_randinit(state);

    /* Check aliasing */
    for (i = 0; i < 1000; i++)
    {
        fmpq_poly_t f, g;
        fmpz_t anum, aden;
        mpq_t a;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_randtest(f, state, n_randint(100), 100);
        
        fmpz_init(anum);
        fmpz_init(aden);
        fmpz_randtest(anum, state, 100);
        fmpz_randtest_not_zero(aden, state, 100);
        
        mpq_init(a);
        fmpz_get_mpz(mpq_numref(a), anum);
        fmpz_get_mpz(mpq_denref(a), aden);
        mpq_canonicalize(a);
        
        fmpq_poly_rescale(g, f, a);
        fmpq_poly_rescale(f, f, a);

        cflags |= fmpq_poly_is_canonical(f) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(g) ? 0 : 2;
        result = (fmpq_poly_equal(f, g) && !cflags);
        if (!result)
        {
            printf("FAIL (aliasing):\n");
            fmpq_poly_print(f), printf("\n\n");
            fmpq_poly_print(g), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpz_clear(anum);
        fmpz_clear(aden);
        mpq_clear(a);
    }

    /* Check that rescaling by a and then by 1/a is the identity */
    for (i = 0; i < 1000; i++)
    {
        fmpq_poly_t f, g;
        fmpz_t anum, aden;
        mpq_t a;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_randtest(f, state, n_randint(100), 100);
        
        fmpz_init(anum);
        fmpz_init(aden);
        fmpz_randtest_not_zero(anum, state, 100);
        fmpz_randtest_not_zero(aden, state, 100);
        
        mpq_init(a);
        fmpz_get_mpz(mpq_numref(a), anum);
        fmpz_get_mpz(mpq_denref(a), aden);
        mpq_canonicalize(a);
        
        fmpq_poly_rescale(g, f, a);
        mpq_inv(a, a);
        fmpq_poly_rescale(g, g, a);

        cflags |= fmpq_poly_is_canonical(f) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(g) ? 0 : 2;
        result = (fmpq_poly_equal(f, g) && !cflags);
        if (!result)
        {
            printf("FAIL (composition of a and 1/a):\n");
            fmpq_poly_print(f), printf("\n\n");
            fmpq_poly_print(g), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpz_clear(anum);
        fmpz_clear(aden);
        mpq_clear(a);
    }

    fmpq_poly_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
