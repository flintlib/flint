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
    slong rep, bits, mbits, max_bits, len, nnz, i;
    fmpz_sparse_vec_t u;
    FLINT_TEST_INIT(state);
    
    flint_printf("max bits....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        do bits = n_randint(state, 200);
        while (bits < UWORD(2));
        len = n_randint(state, 200);
        nnz = n_randint(state, len+1);

        fmpz_sparse_vec_init(u);
        fmpz_sparse_vec_randtest(u, state, nnz, len, bits);
        mbits = fmpz_sparse_vec_max_bits(u);
        mbits = FLINT_ABS(mbits);

        max_bits = 0;
        for (i = 0; i < nnz; ++i)
        {
            bits = fmpz_bits(u->entries[i].val);
            if (bits > mbits)
            {
                flint_printf("FAIL: entry %wd has bitcnt %wd, max was %wd\n", i, bits, mbits);
                fmpz_sparse_vec_print_pretty(u, 0, len);
                abort();
            }
            if (bits > max_bits) max_bits = bits;
        }
        if (mbits != max_bits)
        {
            flint_printf("FAIL: no entry has bitcnt %wd\n", mbits);
            fmpz_sparse_vec_print_pretty(u, 0, len);
            abort();
        }
        fmpz_sparse_vec_clear(u);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
