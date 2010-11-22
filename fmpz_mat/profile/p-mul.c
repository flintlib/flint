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
    Copyright 2010 Fredrik Johansson

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
    ulong dim;
    int algorithm;
    long bits;
} mat_mul_t;


void sample(void * arg, ulong count)
{
    mat_mul_t * params = (mat_mul_t *) arg;
    ulong i, dim = params->dim;
    long bits = params->bits;
    int algorithm = params->algorithm;

    fmpz_randstate_t rnd;
    fmpz_mat_t A, B, C;
    fmpz_randinit(rnd);

    fmpz_mat_init(A, dim, dim);
    fmpz_mat_init(B, dim, dim);
    fmpz_mat_init(C, dim, dim);

    fmpz_mat_randtest(A, rnd, bits);
    fmpz_mat_randtest(B, rnd, bits);

    prof_start();

    if (algorithm == 0)
        for (i = 0; i < count; i++)
            fmpz_mat_mul(C, A, B);
    else if (algorithm == 1)
        for (i = 0; i < count; i++)
            fmpz_mat_mul_classical(C, A, B);
    else if (algorithm == 2)
        for (i = 0; i < count; i++)
            fmpz_mat_mul_multi_mod(C, A, B);

    prof_stop();

    fmpz_mat_clear(A);
    fmpz_mat_clear(B);
    fmpz_mat_clear(C);
    fmpz_randclear(rnd);
}

int main(void)
{
    double min_default, min_classical, min_multi_mod, max;
    mat_mul_t params;
    long dim;

    params.bits = 200;

    printf("fmpz_mat_mul (bits = %ld):\n", params.bits);

    for (dim = 2; dim <= 512; dim = (long) ((double) dim * 1.1) + 1)
    {
        params.dim = dim;

        params.algorithm = 0;
        prof_repeat(&min_default, &max, sample, &params);

        params.algorithm = 1;
        prof_repeat(&min_classical, &max, sample, &params);

        params.algorithm = 2;
        prof_repeat(&min_multi_mod, &max, sample, &params);

        printf("dim = %ld mul/classical/multi_mod %.2f %.2f %.2f (us)\n", 
            dim, min_default, min_classical, min_multi_mod);
    }

    return 0;
}
