/*
    Copyright (C) 2010 Fredrik Johansson

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
#include "ulong_extras.h"

int
main(void)
{
    slong rep, bits, len, nnz, i, j;
    fmpz_t mod, c;
    fmpz_sparse_vec_t u, v;
    FLINT_TEST_INIT(state);
    

    flint_printf("mod....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        do bits = n_randint(state, 10);
        while (bits < UWORD(3));
        len = n_randint(state, 50);
        nnz = n_randint(state, len+1);

        fmpz_init(mod);
        fmpz_init(c);
        fmpz_sparse_vec_init(u);
        fmpz_sparse_vec_init(v);

        fmpz_randbits(mod, state, bits-1);
        fmpz_sparse_vec_randtest(u, state, nnz, len, bits);

        fmpz_sparse_vec_scalar_mod_fmpz(v, u, mod);

        for (i = j = 0; i < u->nnz; ++i)
        {
            fmpz_mod(c, u->entries[i].val, mod);
            if (fmpz_is_zero(c)) continue;
            if (u->entries[i].ind != v->entries[j].ind)
            {
                flint_printf("Fail: v %% c has extra/missing entry\n");
                fmpz_sparse_vec_print_pretty(u, 0, len);
                fmpz_print(mod); flint_printf("\n");
                fmpz_sparse_vec_print_pretty(v, 0, len);
                abort();
            }
            if (fmpz_cmp(c, v->entries[j].val) != 0)
            {
                flint_printf("Fail: v %% c has incorrect entry\n");
                abort();
            }
            ++j;
        }
        if (j != v->nnz)
        {
            flint_printf("Fail: v %% c has incorrect number of nonzeroes\n");
            abort();
        }

        fmpz_sparse_vec_scalar_mods_fmpz(v, u, mod);
        for (i = j = 0; i < u->nnz; ++i)
        {
            fmpz_mods(c, u->entries[i].val, mod);
            if (fmpz_is_zero(c)) continue;
            if (u->entries[i].ind != v->entries[j].ind)
            {
                flint_printf("Fail: v %%s c has extra/missing entry\n");
                fmpz_sparse_vec_print_pretty(u, 0, len);
                fmpz_print(mod); flint_printf("\n");
                fmpz_sparse_vec_print_pretty(v, 0, len);
                abort();
            }
            if (fmpz_cmp(c, v->entries[j].val) != 0)
            {
                flint_printf("Fail: v %%s c has incorrect entry\n");
                abort();
            }
            ++j;
        }
        if (j != v->nnz)
        {
            flint_printf("Fail: v %%s c has incorrect number of nonzeroes\n");
            abort();
        }

        fmpz_sparse_vec_set(v, u, 0);
        fmpz_sparse_vec_scalar_mod_fmpz(v, v, mod);

        for (i = j = 0; i < u->nnz; ++i)
        {
            fmpz_mod(c, u->entries[i].val, mod);
            if (fmpz_is_zero(c)) continue;
            if (u->entries[i].ind != v->entries[j].ind)
            {
                flint_printf("Fail: v %%= c has extra/missing entry\n");
                abort();
            }
            if (fmpz_cmp(c, v->entries[j].val) != 0)
            {
                flint_printf("Fail: v %%= c has incorrect entry\n");
                abort();
            }
            ++j;
        }
        if (j != v->nnz)
        {
            flint_printf("Fail: v %%= c has incorrect number of nonzeroes\n");
            abort();
        }

        fmpz_sparse_vec_set(v, u, 0);
        fmpz_sparse_vec_scalar_mods_fmpz(v, v, mod);

        for (i = j = 0; i < u->nnz; ++i)
        {
            fmpz_mods(c, u->entries[i].val, mod);
            if (fmpz_is_zero(c)) continue;
            if (u->entries[i].ind != v->entries[j].ind)
            {
                flint_printf("Fail: v %%s= c has extra/missing entry\n");
                abort();
            }
            if (fmpz_cmp(c, v->entries[j].val) != 0)
            {
                flint_printf("Fail: v %%s= c has incorrect entry\n");
                abort();
            }
            ++j;
        }
        if (j != v->nnz)
        {
            flint_printf("Fail: v %%s= c has incorrect number of nonzeroes\n");
            abort();
        }

        fmpz_clear(mod);
        fmpz_clear(c);
        fmpz_sparse_vec_clear(u);
        fmpz_sparse_vec_clear(v);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
