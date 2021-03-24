/*
    Copyright 2009 William Hart
    Copyright 2010,2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "flint.h"
#include "fmpz_mat.h"
#include "fmpz.h"
#include "ulong_extras.h"
#include "thread_support.h"
#if FLINT_USES_BLAS
#include "cblas.h"
#endif

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
    slong i, m = params->m, n = params->n, k = params->k;
    slong bits = params->bits;
    int algorithm = params->algorithm;

    fmpz_mat_t A, B, C;
    FLINT_TEST_INIT(state);
    

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
    else if (algorithm == 4)
	for (i = 0; i < count; i++)
	    fmpz_mat_mul_strassen(C, A, B);
#if FLINT_USES_BLAS
    else if (algorithm == 5)
	for (i = 0; i < count; i++)
	    fmpz_mat_mul_blas(C, A, B);
#endif

    prof_stop();

    fmpz_mat_clear(A);
    fmpz_mat_clear(B);
    fmpz_mat_clear(C);
    
    flint_randclear(state);
}

int main(void)
{
    double min_default, min_classical, min_inline, min_multi_mod, min_strassen, min_blas = 0.0, max;
    mat_mul_t params;
    slong bits, dim;
    int threads;

    for (bits = 1; bits <= 1024; bits = (slong) ((double) bits * 1.3) + 1)
    {
        params.bits = bits;

        flint_printf("fmpz_mat_mul (bits = %wd):\n", params.bits);

        for (dim = 1; dim <= (bits < 100 ? 1000 : bits < 500 ? 200 : 100); dim = (slong) ((double) dim * 1.3) + 1)
        {
            flint_printf("dim = %wd default/classical/inline/modular/strassen/blas:\n", dim);
	    
	    for (threads = 1; threads <= 8; threads++)
	    {
#if FLINT_USES_BLAS
    		    openblas_set_num_threads(threads);
#endif
    		    flint_set_num_threads(threads);

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
            
                params.algorithm = 4;
                prof_repeat(&min_strassen, &max, sample, &params);

#if FLINT_USES_BLAS
	        params.algorithm = 5;
	        prof_repeat(&min_blas, &max, sample, &params);
#endif

                flint_printf(" (%d): %.2f %.2f %.2f %.2f %.2f %.2f(us)\n", 
                    threads, min_default, min_classical, min_inline, min_multi_mod, min_strassen, min_blas);
	    }
        }
    }

    return 0;
}
