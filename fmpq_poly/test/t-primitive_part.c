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
    flint_rand_t state;
    ulong cflags = 0UL;

    printf("primitive_part....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t a, b;
        fmpq_poly_t f, g;
        mpq_t x;

        mpq_init(x);
        fmpz_init(a);
        fmpz_init(b);
        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_randtest(f, state, n_randint(state, 100), 200);
        fmpz_randtest(a, state, 100);
        fmpz_randtest_not_zero(b, state, 100);

        fmpz_get_mpz(mpq_numref(x), a);
        fmpz_get_mpz(mpq_denref(x), b);
        mpq_canonicalize(x);
        fmpq_poly_scalar_mul_mpq(f, f, x);

        fmpq_poly_primitive_part(g, f);
        fmpq_poly_primitive_part(f, f);

        cflags |= fmpq_poly_is_canonical(f) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(g) ? 0 : 2;
        result = (fmpq_poly_equal(f, g) && !cflags);
        if (!result)
        {
            printf("FAIL:\n");
            fmpq_poly_debug(f), printf("\n\n");
            fmpq_poly_debug(g), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        mpq_clear(x);
        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
    }

    /* Check that content(f) primitive_part(f) = +- f */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t a, b;
        fmpq_poly_t f, g;
        mpq_t x, y;

        mpq_init(x);
        mpq_init(y);
        fmpz_init(a);
        fmpz_init(b);
        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_randtest(f, state, n_randint(state, 100), 200);
        fmpz_randtest(a, state, 100);
        fmpz_randtest_not_zero(b, state, 100);

        fmpz_get_mpz(mpq_numref(x), a);
        fmpz_get_mpz(mpq_denref(x), b);
        mpq_canonicalize(x);
        fmpq_poly_scalar_mul_mpq(f, f, x);

        fmpq_poly_content(y, f);
        fmpq_poly_primitive_part(g, f);
        fmpq_poly_scalar_mul_mpq(g, g, y);

        if (!fmpq_poly_is_zero(f) && fmpz_sgn(f->coeffs + (f->length - 1)) < 0)
            fmpq_poly_neg(g, g);

        cflags |= fmpq_poly_is_canonical(f) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(g) ? 0 : 2;
        result = (fmpq_poly_equal(f, g) && !cflags);
        if (!result)
        {
            printf("FAIL:\n");
            fmpq_poly_debug(f), printf("\n\n");
            fmpq_poly_debug(g), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        mpq_clear(x);
        mpq_clear(y);
        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
