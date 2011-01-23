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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "ulong_extras.h"


int
main(void)
{
    nmod_mat_t A1, A2, Ar;
    long i, m, n, rank1, rank2, rank;
    mp_limb_t mod;
    flint_rand_t state;
    flint_randinit(state);

    printf("rowreduce_1/2....");
    fflush(stdout);

    for (i = 0; i < 100; i++)
    {
        m = n_randint(state, 200);
        n = n_randint(state, 200);
        mod = _nmod_mat_fast_rowreduce_modulus_1(m, n, 0);

        nmod_mat_init(A1, m, n, mod);
        nmod_mat_randtest(A1, state);
        nmod_mat_init_set(Ar, A1);

        rank1 = _nmod_mat_rowreduce_1(A1, 0);
        rank = _nmod_mat_rowreduce_r(Ar, 0);

        if (rank1 != rank)
        {
            printf("FAIL:\n");
            printf("rank %ld != %ld!\n", rank1, rank);
            abort();
        }

        if (!nmod_mat_equal(A1, Ar))
        {
            printf("FAIL:\n");
            printf("A1 != Ar!\n");
            abort();
        }

        nmod_mat_clear(A1);
        nmod_mat_clear(Ar);
    }

    for (i = 0; i < 100; i++)
    {
        m = n_randint(state, 200);
        n = n_randint(state, 200);
        mod = _nmod_mat_fast_rowreduce_modulus_2(m, n, 0);

        nmod_mat_init(A2, m, n, mod);
        nmod_mat_randtest(A2, state);
        nmod_mat_init_set(Ar, A2);

        rank2 = _nmod_mat_rowreduce_2(A2, 0);
        rank = _nmod_mat_rowreduce_r(Ar, 0);

        if (rank2 != rank)
        {
            printf("FAIL:\n");
            printf("rank %ld != %ld!\n", rank2, rank);
            abort();
        }

        if (!nmod_mat_equal(A2, Ar))
        {
            printf("FAIL:\n");
            printf("A2 != Ar!\n");
            abort();
        }

        nmod_mat_clear(A2);
        nmod_mat_clear(Ar);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
