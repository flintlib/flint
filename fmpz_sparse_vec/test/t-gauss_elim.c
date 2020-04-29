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
    fmpz_t a, b, c, d, g;
    fmpz_sparse_vec_t u, v, w, x, y;
    FLINT_TEST_INIT(state);
    
    flint_printf("Gaussian elimination....");
    fflush(stdout);

    for (rep = 0; rep < 10; rep++)
    {
        do bits = n_randint(state, 10);
        while (bits < UWORD(2));
        len = n_randint(state, 20);
        nnz = n_randint(state, len+1);

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);
        fmpz_init(g);
        fmpz_sparse_vec_init(u);
        fmpz_sparse_vec_init(v);
        fmpz_sparse_vec_init(w);
        fmpz_sparse_vec_init(x);
        fmpz_sparse_vec_init(y);

        fmpz_sparse_vec_randtest(u, state, nnz, len, bits);
        fmpz_sparse_vec_randtest(v, state, nnz, len, bits);
        fmpz_sparse_vec_set(w, u, 0);
        fmpz_sparse_vec_set(x, v, 0);
        fmpz_sparse_vec_gauss_elim(x, w);
        if (!fmpz_sparse_vec_equal(w, u, 0))
        {
            flint_printf("FAIL: pivot row modified\n");
            abort();
        }
        if (nnz != 0 && fmpz_sparse_vec_at(v, u->entries[0].ind) != NULL)
        {
            /* Undo elimination and check if matches original */
            fmpz_fdiv_q(c, *fmpz_sparse_vec_at(v, u->entries[0].ind), u->entries[0].val);
            fmpz_sparse_vec_scalar_addmul_fmpz(w, x, u, c);
            if (!fmpz_sparse_vec_equal(w, v, 0))
            {
                flint_printf("FAIL: Incorrect Gaussian elimination\n");
                fmpz_sparse_vec_print_pretty(u, 0, len);
                fmpz_sparse_vec_print_pretty(v, 0, len);
                fmpz_sparse_vec_print_pretty(x, 0, len);
                fmpz_sparse_vec_print_pretty(w, 0, len);
                abort();
            }
        }
        else if (!fmpz_sparse_vec_equal(x, v, 0))
        {
                flint_printf("FAIL: Incorrect Gaussian elimination on trivial vectors\n");
                fmpz_sparse_vec_print_pretty(x, 0, len);
                abort();
        }
        
        if (nnz != 0 && u->entries[0].ind == v->entries[0].ind)
        {
            fmpz_sparse_vec_set(w, u, 0);
            fmpz_sparse_vec_set(x, v, 0);
            fmpz_sparse_vec_gauss_elim_ext(x, w);
            
            fmpz_xgcd(g, a, b, u->entries[0].val, v->entries[0].val);
            fmpz_divexact(d, u->entries[0].val, g);
            fmpz_divexact(c, v->entries[0].val, g);
            fmpz_sparse_vec_scalar_mul_fmpz(y, w, d);
            fmpz_sparse_vec_scalar_submul_fmpz(y, y, x, b);
            if (!fmpz_sparse_vec_equal(y, u, 0))
            {
                flint_printf("FAIL: did not recover u when inverting extended Gaussian elimination\n");
                fmpz_sparse_vec_print_pretty(u, 0, len);
                fmpz_sparse_vec_print_pretty(v, 0, len);
                fmpz_sparse_vec_print_pretty(w, 0, len);
                fmpz_sparse_vec_print_pretty(x, 0, len);
                fmpz_sparse_vec_print_pretty(y, 0, len);
                abort();
            }
            fmpz_sparse_vec_scalar_mul_fmpz(y, w, c);
            fmpz_sparse_vec_scalar_addmul_fmpz(y, y, x, a);
            if (!fmpz_sparse_vec_equal(y, v, 0))
            {
                flint_printf("FAIL: did not recover v when inverting extended Gaussian elimination\n");
                fmpz_sparse_vec_print_pretty(u, 0, len);
                fmpz_sparse_vec_print_pretty(v, 0, len);
                fmpz_sparse_vec_print_pretty(w, 0, len);
                fmpz_sparse_vec_print_pretty(x, 0, len);
                fmpz_sparse_vec_print_pretty(y, 0, len);
                abort();
            }
        }
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_clear(g);
        fmpz_sparse_vec_clear(u);
        fmpz_sparse_vec_clear(v);
        fmpz_sparse_vec_clear(w);
        fmpz_sparse_vec_clear(x);
        fmpz_sparse_vec_clear(y);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
