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
    slong n, c;
    nmod_t mod;
    nmod_sparse_vec_t umod;
    fmpz_sparse_vec_t u;
    FLINT_TEST_INIT(state);
    

    flint_printf("nmod....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        do bits = n_randint(state, 256);
        while (bits < UWORD(2));
        len = n_randint(state, 50);
        nnz = n_randint(state, len+1);

        do n = n_randlimb(state);
        while (n < UWORD(2));
        nmod_init(&mod, n);

        nmod_sparse_vec_init(umod);
        fmpz_sparse_vec_init(u);

        fmpz_sparse_vec_randtest(u, state, nnz, len, bits);
        fmpz_sparse_vec_get_nmod_sparse_vec(umod, u, mod);

        for (i = j = 0; i < u->nnz; ++i)
        {
            c = fmpz_fdiv_ui(u->entries[i].val, n);
            if (c == UWORD(0)) continue;
            if (u->entries[i].ind != umod->entries[j].ind)
            {
                flint_printf("Fail: v %% c has extra/missing entry\n");
                abort();
            }
            if (c != umod->entries[j].val)
            {
                flint_printf("Fail: v %% c has incorrect entry\n");
                abort();
            }
            ++j;
        }
        if (j != umod->nnz)
        {
            flint_printf("Fail: v %% c has incorrect number of nonzeroes\n");
            abort();
        }

        nmod_sparse_vec_randtest(umod, state, nnz, len, mod);
        fmpz_sparse_vec_set_nmod_sparse_vec_unsigned(u, umod);

        for (i = 0; i < u->nnz; ++i)
        {
            if (u->entries[i].ind != umod->entries[i].ind ||
                fmpz_cmp_ui(u->entries[i].val, umod->entries[i].val))
            {
                flint_printf("FAIL: u = umod has bad entry\n");
                abort();
            }
        }
        if (u->nnz != umod->nnz)
        {
            flint_printf("Fail: u = umod has bad length\n");
            abort();
        }

        fmpz_sparse_vec_set_nmod_sparse_vec(u, umod, mod);

        for (i = 0; i < u->nnz; ++i)
        {
            if (u->entries[i].ind != umod->entries[i].ind ||
                (fmpz_sgn(u->entries[i].val) > 0 && fmpz_cmp_ui(u->entries[i].val, umod->entries[i].val)) ||
                (fmpz_sgn(u->entries[i].val) < 0 && fmpz_cmp_si(u->entries[i].val, umod->entries[i].val - n)))
            {
                flint_printf("FAIL: su = umod has bad entry\n");
                abort();
            }
        }
        if (u->nnz != umod->nnz)
        {
            flint_printf("Fail: su = umod has bad length\n");
            abort();
        }

        fmpz_sparse_vec_clear(u);
        nmod_sparse_vec_clear(umod);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
