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
    Copyright 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "profiler.h"
#include "flint.h"
#include "fmpz_mat.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "arith.h"

typedef struct
{
    ulong n;
    int algorithm;
} bernoulli_vec_t;


void sample(void * arg, ulong count)
{
    fmpz * num;
    fmpz * den;
    bernoulli_vec_t * params = (bernoulli_vec_t *) arg;
    ulong n = params->n;
    len_t i;
    int algorithm = params->algorithm;

    num = _fmpz_vec_init(n);
    den = _fmpz_vec_init(n);

    prof_start();

    for (i = 0; i < count; i++)
    {
        if (algorithm == 0)
        {
            _arith_bernoulli_number_vec_recursive(num, den, n);
        }
        else if (algorithm == 1)
        {
            _arith_bernoulli_number_vec_multi_mod(num, den, n);
        }
        else if (algorithm == 2)
        {
            _arith_bernoulli_number_vec_zeta(num, den, n);
            mpfr_free_cache();
        }
    }

    prof_stop();

    _fmpz_vec_clear(num, n);
    _fmpz_vec_clear(den, n);
}

int main(void)
{
    double min_recursive, min_multi_mod, min_zeta, max;
    bernoulli_vec_t params;
    len_t n;

    printf("n / recursive / multi_mod / zeta / best [times in us]\n");

    for (n = 2; n <= 10000; n = (len_t) ((double) n * 1.2) + 1)
    {
        params.n = n;

        if (n < 1500)
        {
            params.algorithm = 0;
            prof_repeat(&min_recursive, &max, sample, &params);
        }
        else
            min_recursive = 0.0;

        params.algorithm = 1;
        prof_repeat(&min_multi_mod, &max, sample, &params);

        params.algorithm = 2;
        prof_repeat(&min_zeta, &max, sample, &params);

        printf("%ld %.2f %.2f %.2f ", 
            n, min_recursive, min_multi_mod, min_zeta);

        if (min_recursive && min_recursive < min_multi_mod && \
            min_recursive < min_zeta)
            printf("(recursive)\n");
        else if (min_multi_mod < min_zeta)
            printf("(multi_mod)\n");
        else
            printf("(zeta)\n");
    }

    return 0;
}
