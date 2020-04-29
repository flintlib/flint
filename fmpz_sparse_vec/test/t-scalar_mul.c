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
    fmpz_sparse_vec_t u, v, w, x;
    FLINT_TEST_INIT(state);
    
    flint_printf("scalar_mul and muladd....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        do bits = n_randint(state, 200);
        while (bits < UWORD(2));
        len = n_randint(state, 50);
        nnz = n_randint(state, len+1);

        fmpz_init(c);
        fmpz_sparse_vec_init(u);
        fmpz_sparse_vec_init(v);
        fmpz_sparse_vec_init(w);
        fmpz_sparse_vec_init(x);

        fmpz_randtest(c, state, bits);
        fmpz_sparse_vec_randtest(u, state, nnz, len, bits);
        fmpz_sparse_vec_randtest(v, state, nnz, len, bits);

        fmpz_sparse_vec_scalar_addmul_fmpz(w, u, v, c);
        fmpz_sparse_vec_scalar_mul_fmpz(x, v, c);
        fmpz_sparse_vec_add(x, x, u);

        if (!fmpz_sparse_vec_equal(w, x, 0))
        {
            flint_printf("FAIL: u + c*v != u + (c*v)\n");
            abort();
        }

        fmpz_sparse_vec_scalar_addmul_fmpz(w, u, v, c);
        fmpz_sparse_vec_scalar_addmul_fmpz(u, u, v, c);

        if (!fmpz_sparse_vec_equal(u, w, 0))
        {
            flint_printf("FAIL: u + c*v != (u += c*v)\n");
            abort();
        }

        fmpz_sparse_vec_scalar_addmul_fmpz(w, u, v, c);
        fmpz_sparse_vec_scalar_addmul_fmpz(v, u, v, c);

        if (!fmpz_sparse_vec_equal(v, w, 0))
        {
            flint_printf("FAIL: u + c*v != (u += c*v)\n");
            abort();
        }

        fmpz_sparse_vec_scalar_submul_fmpz(w, u, v, c);
        fmpz_sparse_vec_scalar_submul_fmpz(u, u, v, c);

        if (!fmpz_sparse_vec_equal(u, w, 0))
        {
            flint_printf("FAIL: u + c*v != (u += c*v)\n");
            abort();
        }

        fmpz_sparse_vec_scalar_submul_fmpz(w, u, v, c);
        fmpz_sparse_vec_scalar_submul_fmpz(v, u, v, c);

        if (!fmpz_sparse_vec_equal(v, w, 0))
        {
            flint_printf("FAIL: u + c*v != (u += c*v)\n");
            abort();
        }

        fmpz_sparse_vec_scalar_mul_fmpz(x, v, c);
        fmpz_sparse_vec_scalar_mul_fmpz(v, v, c);

        if (!fmpz_sparse_vec_equal(v, x, 0))
        {
            flint_printf("FAIL: c*v != (c *= v)\n");
            abort();
        }
        fmpz_clear(c);
        fmpz_sparse_vec_clear(u);
        fmpz_sparse_vec_clear(v);
        fmpz_sparse_vec_clear(w);
        fmpz_sparse_vec_clear(x);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
