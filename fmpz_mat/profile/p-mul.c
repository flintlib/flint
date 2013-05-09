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

    Copyright 2009 William Hart
    Copyright 2010,2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "flint.h"
#include "fmpz_mat.h"
#include "fmpz.h"
#include "ulong_extras.h"

typedef struct
{
    len_t m;
    len_t n;
    len_t k;
    int algorithm;
    len_t bits;
} mat_mul_t;


void sample(void * arg, ulong count)
{
    mat_mul_t * params = (mat_mul_t *) arg;
    len_t i, m = params->m, n = params->n, k = params->k;
    len_t bits = params->bits;
    int algorithm = params->algorithm;

    flint_rand_t rnd;
    fmpz_mat_t A, B, C;
    flint_rand_t state;
    flint_randinit(state);

    fmpz_mat_init(A, m, n);
    fmpz_mat_init(B, n, k);
    fmpz_mat_init(C, m, k);

    fmpz_mat_randbits(A, state, bits);
    fmpz_mat_randbits(B, state, bits);

    prof_start();

    if (algorithm == 0)
        for (i = 0; i < count; i++)
            fmpz_mat_mul(C, A, B);
    else if (algorithm == 1)
        for (i = 0; i < count; i++)
            fmpz_mat_mul_classical(C, A, B);
    else if (algorithm == 2)
        for (i = 0; i < count; i++)
            fmpz_mat_mul_classical_inline(C, A, B);
    else if (algorithm == 3)
        for (i = 0; i < count; i++)
            fmpz_mat_mul_multi_mod(C, A, B);

    prof_stop();

    fmpz_mat_clear(A);
    fmpz_mat_clear(B);
    fmpz_mat_clear(C);
    
    flint_randclear(state);
}

int main(void)
{
    double min_default, min_classical, min_inline, min_multi_mod, max;
    mat_mul_t params;
    len_t bits, dim;

    for (bits = 1; bits <= 2000; bits = (len_t) ((double) bits * 1.3) + 1)
    {
        params.bits = bits;

        printf("fmpz_mat_mul (bits = %ld):\n", params.bits);

        for (dim = 1; dim <= 512; dim = (len_t) ((double) dim * 1.3) + 1)
        {
            params.m = dim;
            params.n = dim;
            params.k = dim;

            params.algorithm = 0;
            prof_repeat(&min_default, &max, sample, &params);

            params.algorithm = 1;
            prof_repeat(&min_classical, &max, sample, &params);

            params.algorithm = 2;
            prof_repeat(&min_inline, &max, sample, &params);

            params.algorithm = 3;
            prof_repeat(&min_multi_mod, &max, sample, &params);

            printf("dim = %ld default/classical/inline/multi_mod %.2f %.2f %.2f %.2f (us)\n", 
                dim, min_default, min_classical, min_inline, min_multi_mod);

            if (min_multi_mod < 0.6*min_default)
                printf("BAD!\n");

            if (min_inline < 0.6*min_default)
                printf("BAD!\n");

            if (min_multi_mod < 0.7*min_inline)
                break;
        }
    }

    return 0;
}
