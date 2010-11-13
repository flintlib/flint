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
#include <limits.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    long i;

    printf("mul....");
    fflush(stdout);

    for (i = 0; i < 10000; i++)
    {
        nmod_mat_t A, B, C, D;

        long m, k, n;

        m = n_randint(50);
        k = n_randint(50);
        n = n_randint(50);

        /* We want to generate matrices with many entries close to half
           or full limbs with high probability, to stress overflow handling */
        switch (n_randint(3))
        {
            case 0:
                nmod_mat_init(A, m, n, n_randtest_not_zero());
                break;
            case 1:
                nmod_mat_init(A, m, n, ULONG_MAX/2 + 1 - n_randbits(4));
                break;
            case 2:
                nmod_mat_init(A, m, n, ULONG_MAX - n_randbits(4));
                break;
        }

        switch (n_randint(3))
        {
            case 0:
                nmod_mat_init(B, n, k, n_randtest_not_zero());
                break;
            case 1:
                nmod_mat_init(B, n, k, ULONG_MAX/2 + 1 - n_randbits(4));
                break;
            case 2:
                nmod_mat_init(B, n, k, ULONG_MAX - n_randbits(4));
                break;
        }

        /* Small modulus with high probability, again to trigger reductions */
        switch (n_randint(2))
        {
            case 0:
                nmod_mat_init(C, m, k, n_randtest_not_zero());
                break;
            case 1:
                nmod_mat_init(C, m, k, 1+n_randbits(4));
                break;
        }

        nmod_mat_init(D, m, k, C->mod.n);

        if (n_randint(2))
            nmod_mat_randtest(A);
        else
            nmod_mat_randfull(A);

        if (n_randint(2))
            nmod_mat_randtest(B);
        else
            nmod_mat_randfull(B);

        nmod_mat_mul(C, A, B);
        _nmod_mat_mul_3(D, A, B);

        if (!nmod_mat_equal(C, D))
        {
            printf("FAIL: results not equal\n");
            nmod_mat_print_pretty(A);
            nmod_mat_print_pretty(B);
            nmod_mat_print_pretty(C);
            nmod_mat_print_pretty(D);
            abort();
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(D);
    }

    printf("PASS\n");
    return 0;
}
