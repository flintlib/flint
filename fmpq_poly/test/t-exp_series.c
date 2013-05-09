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
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    ulong cflags = 0UL;

    printf("exp_series....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and c */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        len_t n = n_randint(state, 50) + 1;

        fmpq_poly_init(a);
        fmpq_poly_init(b);

        fmpq_poly_randtest_not_zero(a, state, n_randint(state, 50) + 1, 50);
        fmpq_poly_set_coeff_ui(a, 0, 0UL);

        fmpq_poly_canonicalise(a);

        fmpq_poly_exp_series(b, a, n);
        fmpq_poly_exp_series(a, a, n);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        result = (fmpq_poly_equal(a, b) && !cflags);
        if (!result)
        {
            printf("FAIL:\n");
            fmpq_poly_debug(a), printf("\n\n");
            fmpq_poly_debug(b), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Check exp(a+b) = exp(a) * exp(b) */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, ab, expa, expb, expab, expa_expb;
        len_t n = n_randint(state, 80) + 1;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(ab);
        fmpq_poly_init(expa);
        fmpq_poly_init(expb);
        fmpq_poly_init(expab);
        fmpq_poly_init(expa_expb);

        fmpq_poly_randtest_not_zero(a, state, n_randint(state, 60) + 1, 80);
        fmpq_poly_set_coeff_ui(a, 0, 0UL);

        fmpq_poly_randtest_not_zero(b, state, n_randint(state, 60) + 1, 80);
        fmpq_poly_set_coeff_ui(b, 0, 0UL);

        fmpq_poly_add(ab, a, b);

        fmpq_poly_exp_series(expab, ab, n);
        fmpq_poly_exp_series(expa, a, n);
        fmpq_poly_exp_series(expb, b, n);
        fmpq_poly_mullow(expa_expb, expa, expb, n);

        cflags |= fmpq_poly_is_canonical(expa) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(expb) ? 0 : 2;
        cflags |= fmpq_poly_is_canonical(expab) ? 0 : 4;
        result = (fmpq_poly_equal(expab, expa_expb) && !cflags);
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), fmpq_poly_debug(a), printf("\n\n");
            printf("b = "), fmpq_poly_debug(b), printf("\n\n");
            printf("exp(a) = "), fmpq_poly_debug(expa), printf("\n\n");
            printf("exp(b) = "), fmpq_poly_debug(expb), printf("\n\n");
            printf("exp(ab) = "), fmpq_poly_debug(expab), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(ab);
        fmpq_poly_clear(expa);
        fmpq_poly_clear(expb);
        fmpq_poly_clear(expab);
        fmpq_poly_clear(expa_expb);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
