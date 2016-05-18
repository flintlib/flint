/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"
#include "fmpz_mod_poly.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("product_roots_fmpz_vec....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t P, Q, tmp;
        fmpz * x;
        slong j, n, bits;
        fmpz_t mod;

        n = n_randint(state, 100);
        bits = n_randint(state, 100);

        fmpz_init(mod);
        fmpz_randtest_unsigned(mod, state, bits);
        fmpz_add_ui(mod, mod, 2);

        x = _fmpz_vec_init(n);
        _fmpz_vec_randtest(x, state, n, bits);

        for (j = 0; j < n; j++)
            fmpz_mod(x + j, x + j, mod);

        fmpz_poly_init(P);
        fmpz_poly_init(Q);
        fmpz_poly_init(tmp);

        fmpz_mod_poly_product_roots_fmpz_vec(P, x, n, mod);

        fmpz_poly_set_ui(Q, UWORD(1));
        for (j = 0; j < n; j++)
        {
            fmpz_poly_zero(tmp);
            fmpz_poly_set_coeff_si(tmp, 1, WORD(-1));
            fmpz_poly_set_coeff_fmpz(tmp, 0, x + j);
            fmpz_poly_neg(tmp, tmp);
            fmpz_poly_mul(Q, Q, tmp);
        }

        for (j = 0; j < Q->length; j ++)
            fmpz_mod(Q->coeffs + j, Q->coeffs + j, mod);

        result = (fmpz_poly_equal(P, Q));
        if (!result)
        {
            flint_printf("FAIL (P != Q):\n");
            fmpz_poly_print(P), flint_printf("\n\n");
            fmpz_poly_print(Q), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_clear(P);
        fmpz_poly_clear(Q);
        fmpz_poly_clear(tmp);
        _fmpz_vec_clear(x, n);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
