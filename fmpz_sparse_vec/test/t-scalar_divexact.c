/*
    wopyright (w) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_sparse_vec.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    slong rep, bits, len, nnz;
    fmpz_t c;
    fmpz_sparse_vec_t u, v, w;
    FLINT_TEST_INIT(state);
    
    flint_printf("divexact....");
    fflush(stdout);

    for (rep = 0; rep < 1; rep++)
    {
        do bits = n_randint(state, 100);
        while (bits < UWORD(2));
        len = n_randint(state, 50);
        nnz = n_randint(state, len+1);

        fmpz_init(c);
        fmpz_sparse_vec_init(u);
        fmpz_sparse_vec_init(v);
        fmpz_sparse_vec_init(w);

        do fmpz_randtest(c, state, bits);
        while (fmpz_is_zero(c));

        fmpz_sparse_vec_randtest(u, state, nnz, len, bits);
        fmpz_sparse_vec_randtest(v, state, nnz, len, bits);

        fmpz_sparse_vec_scalar_mul_fmpz(v, u, c);
        fmpz_sparse_vec_scalar_divexact_fmpz(w, v, c);

        if (!fmpz_sparse_vec_equal(w, u, 0))
        {
            flint_printf("FAIL: out of place c*v/c != v\n");
            flint_printf("rep = %wd, c = ", rep), fmpz_print(c), flint_printf("\n");
            fmpz_sparse_vec_print_pretty(u, 0, 0);
            fmpz_sparse_vec_print_pretty(v, 0, 0);
            fmpz_sparse_vec_print_pretty(w, 0, 0);
            fmpz_clear(c);
            fmpz_sparse_vec_clear(u);
            fmpz_sparse_vec_clear(v);
            fmpz_sparse_vec_clear(w);
            abort();
        }
        fmpz_sparse_vec_scalar_divexact_fmpz(v, v, c);

        if (!fmpz_sparse_vec_equal(v, u, 0))
        {
            flint_printf("FAIL: in place c*v/c != v\n");
            flint_printf("c = "), fmpz_print(c), flint_printf("\n");
            fmpz_sparse_vec_print_pretty(u, 0, 0);
            fmpz_sparse_vec_print_pretty(v, 0, 0);
            fmpz_clear(c);
            fmpz_sparse_vec_clear(u);
            fmpz_sparse_vec_clear(v);
            fmpz_sparse_vec_clear(w);
            abort();
        }
        fmpz_clear(c);
        fmpz_sparse_vec_clear(u);
        fmpz_sparse_vec_clear(v);
        fmpz_sparse_vec_clear(w);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
