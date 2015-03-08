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
    slong m;
    slong n;
    slong k;
    int algorithm;
    slong bits;
} mat_mul_t;


void sample(void * arg, ulong count)
{
    mat_mul_t * params = (mat_mul_t *) arg;
    slong i, n = params->n, m = params->m ;
    slong bits = params->bits;
    int algorithm = params->algorithm;

    flint_rand_t rnd;
    fmpz_mat_t A, C;
    FLINT_TEST_INIT(state);
    

    fmpz_mat_init(A, m, n);
    fmpz_mat_init(C, m, n);

    fmpz_mat_randbits(A, state, bits);

    prof_start();

    if (algorithm == 0)
        for (i = 0; i < count; i++)
            fmpz_mat_sqr_classical(C, A);
    else if (algorithm == 1)
        for (i = 0; i < count; i++)
            fmpz_mat_sqr(C, A);
    
    prof_stop();

    fmpz_mat_clear(A);
    fmpz_mat_clear(C);
    
    flint_randclear(state);
}

int main(void)
{
    double min_default, min_classical, min_inline, min_multi_mod, max;
    mat_mul_t params;
    slong bits, dim;

    for (bits = 1; bits <= 2000; bits = (slong) ((double) bits) + 100)
    {
        params.bits = bits;

        flint_printf("bits = %wd :\n", params.bits);

        for (dim = 1; dim <= 450; dim = (slong) ((double) dim * 3) + 1)
        {
            params.n = dim;
            params.m = dim;

            params.algorithm = 0;
            prof_repeat(&min_default, &max, sample, &params);

            params.algorithm = 1;
            prof_repeat(&min_classical, &max, sample, &params);

            flint_printf("dim = %wd Classical : %.2f Bodrato : %.2f        ", dim, min_default, min_classical);
            if (min_default < 0.8*min_classical)
                flint_printf("BAD! \n");
            else if(min_default <= min_classical)
                flint_printf("OK! \n");
            else if (min_default > min_classical)
                flint_printf("GREAT! \n");

        }
    }

    return 0;
}
