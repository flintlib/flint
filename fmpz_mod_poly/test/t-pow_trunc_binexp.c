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
    Copyright (C) 2012 Lina Kulakova

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    flint_randinit(state);

    printf("pow_trunc_binexp....");
    fflush(stdout);

    /* Check aliasing */
    for (i = 0; i < 10000; i++)
    {
        fmpz_mod_poly_t a, b, c;
        fmpz_t p;
        len_t e, trunc;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 30));
        e = n_randint(state, 20);
        trunc = n_randint(state, 30);

        fmpz_mod_poly_set(c, a);

        fmpz_mod_poly_pow_trunc_binexp(b, a, e, trunc);
        fmpz_mod_poly_pow_trunc_binexp(c, c, e, trunc);

        result = (fmpz_mod_poly_equal(b, c));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a->length = %ld, p = %lu, exp = %ld, trunc = %ld\n",
                a->length, a->p, e, trunc);
            printf("a:\n"); fmpz_mod_poly_print(a), printf("\n\n");
            printf("b:\n"); fmpz_mod_poly_print(b), printf("\n\n");
            printf("c:\n"); fmpz_mod_poly_print(c), printf("\n\n");
            abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(c);
    }

    /* Check powering against naive method */
    for (i = 0; i < 10000; i++)
    {
        fmpz_mod_poly_t a, b, c;
        fmpz_t p;
        len_t e, trunc;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 30));
        e = n_randint(state, 20);
        trunc = n_randint(state, 30);

        fmpz_mod_poly_pow_trunc_binexp(b, a, e, trunc);
        fmpz_mod_poly_pow(c, a, e);
        fmpz_mod_poly_truncate(c, trunc);
        
        result = (fmpz_mod_poly_equal(b, c)
            || (a->length == 0 && e == 0 && c->length == 1 && c->coeffs[0] == 1));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a->length = %ld, p = %lu, exp = %ld, trunc = %ld\n",
                a->length, a->p, e, trunc);
            printf("a:\n"); fmpz_mod_poly_print(a), printf("\n\n");
            printf("b:\n"); fmpz_mod_poly_print(b), printf("\n\n");
            printf("c:\n"); fmpz_mod_poly_print(c), printf("\n\n");
            abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(c);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
