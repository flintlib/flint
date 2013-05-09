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
    ulong dim;
    int algorithm;
    len_t bits;
} mat_mul_t;


void sample(void * arg, ulong count)
{
    mat_mul_t * params = (mat_mul_t *) arg;
    ulong i, dim = params->dim;
    len_t bits = params->bits;
    int algorithm = params->algorithm;

    flint_rand_t rnd;
    fmpz_mat_t A;
    fmpz_t d;
    flint_rand_t state;
    flint_randinit(state);
   
    fmpz_mat_init(A, dim, dim);
    fmpz_init(d);

    fmpz_mat_randtest(A, state, bits);

    prof_start();

    if (algorithm == 0)
        for (i = 0; i < count; i++)
            fmpz_mat_det_bareiss(d, A);
    else if (algorithm == 1)
        for (i = 0; i < count; i++)
            fmpz_mat_det_modular(d, A, 1);
    else if (algorithm == 2)
        for (i = 0; i < count; i++)
            fmpz_mat_det_modular_accelerated(d, A, 1);

    prof_stop();

    fmpz_mat_clear(A);
    fmpz_clear(d);
    
    flint_randclear(state);
}

int main(void)
{
    double min_default, min_classical, min_modular, min_modular_2, max;
    mat_mul_t params;
    len_t dim, bits;

    params.bits = 200;

    for (bits = 2; bits <= 4096; bits *= 2)
    {
        params.bits = bits;
        printf("fmpz_mat_det (bits = %ld):\n", params.bits);

        for (dim = 2; dim <= 512; dim = (len_t) ((double) dim * 1.1) + 1)
        {
            params.dim = dim;

            params.algorithm = 0;
            prof_repeat(&min_classical, &max, sample, &params);

            params.algorithm = 1;
            prof_repeat(&min_modular, &max, sample, &params);

            params.algorithm = 2;
            prof_repeat(&min_modular_2, &max, sample, &params);

            printf("dim = %ld classical/modular/acc. %.2f %.2f %.2f (us)\n", 
                dim, min_classical, min_modular, min_modular_2);

            if (min_modular > 1.1*min_modular_2)
                break;
        }
    }

    return 0;
}
