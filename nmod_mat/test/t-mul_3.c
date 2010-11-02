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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    long i, j;
    printf("mul_3....");
    fflush(stdout);

    /* Full-limb entries */
    for (i = 0; i < 2000; i++)
    {
        nmod_mat_t A, B, C, D, E, F;
        mp_limb_t mod = ((mp_limb_t) -1) - n_randbits(FLINT_BITS - 1);

        long m = n_randint(50);
        long n = n_randint(50);
        long k = n_randint(50);

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(B, n, k, mod);
        nmod_mat_init(C, m, k, mod);
        nmod_mat_init(D, m, k, mod);
        nmod_mat_init(E, m, k, mod);
        nmod_mat_init(F, m, k, mod);

        /* TODO: randtest doesn't generate large enough values */
        for (j = 0; j < m * n; j++)
            A->entries[j] = mod - 1 - n_randint(1+mod/2);
        for (j = 0; j < n * k; j++)
            B->entries[j] = mod - 1 - n_randint(1+mod/2);

        _nmod_mat_mul_classical(C, A, B);

        /* Test of the test: _nmod_mat_mul_classical_2(D, A, B); should fail */
        _nmod_mat_mul_classical_3(D, A, B);

        _nmod_mat_mul_blocked_3(E, A, B, 1+n_randint(50));
        _nmod_mat_mul_blocked_3b(F, A, B, 1+n_randint(50));

        if (!nmod_mat_equal(C, D))
        {
            printf("FAIL: _nmod_mat_mul_classical_3\n");
            abort();
        }

        if (!nmod_mat_equal(C, E))
        {
            printf("FAIL: _nmod_mat_mul_blocked_3\n");
            abort();
        }

        if (!nmod_mat_equal(C, F))
        {
            printf("FAIL: _nmod_mat_mul_blocked_3b\n");
            abort();
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(D);
        nmod_mat_clear(E);
        nmod_mat_clear(F);
    }

    printf("PASS\n");
    return 0;
}
