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

    Copyright (C) 2007 David Harvey, William Hart
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

/*
    Demo FLINT program for computing the q-expansion of the delta function.
*/

#include <stdlib.h>
#include <stdio.h>
#include <mpir.h>
#include <math.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

int main(int argc, char* argv[])
{
    fmpz_t c;
    fmpz_poly_t F2, F4, F8;
    long i, j, N, stop, *values;

    if (argc != 2)
    {
        printf("Syntax: delta_qexp <integer>\n");
        printf("where <integer> is the number of terms to compute\n");
        return EXIT_FAILURE;
    }

    /* Number of terms to compute */
    N = atoi(argv[1]);

    if (N < 1)
    {
        printf("Syntax: delta_qexp <integer>\n");
        printf("where <integer> is the number of terms to compute\n");
        return EXIT_FAILURE;
    }

    /* Compute coefficients of F(q)^2 */
    values = calloc(N, sizeof(long));

    stop = (long) ceil((-1.0 + sqrt(1.0 + 8.0 * N)) / 2.0);

    for (i = 0; i <= stop; i++)
    {
        long index1, index2, value1, value2;

        index1 = i * (i + 1) / 2;
        value1 = (i & 1) ? (-2 * i - 1) : (2 * i + 1);
        for (j = 0; j <= stop; j++)
        {
            index2 = j * (j + 1) / 2;
            if (index1 + index2 >= N)
                break;
            value2 = (j & 1) ? (-2 * j - 1) : (2 * j + 1);
            values[index1 + index2] += value1 * value2;
        }
    }

    /* Create some polynomial objects */
    fmpz_init(c);
    fmpz_poly_init(F2);
    fmpz_poly_init(F4);
    fmpz_poly_init(F8);

    fmpz_poly_fit_length(F2, N);

    for (i = 0; i < N; i++)
      fmpz_poly_set_coeff_si(F2, i, values[i]);

    free(values);

    /* Compute F^4, truncated to length N */
    fmpz_poly_mullow(F4, F2, F2, N);

    /* Compute F^8, truncated to length N */
    fmpz_poly_mullow(F8, F4, F4, N);

    /* Print out last coefficient */
    fmpz_poly_get_coeff_fmpz(c, F8, N - 1);
    printf("Coefficient of q^%ld is ", N);
    fmpz_print(c);
    printf("\n");

    /* Clean up */
    fmpz_clear(c);
    fmpz_poly_clear(F2);
    fmpz_poly_clear(F4);
    fmpz_poly_clear(F8);

    return EXIT_SUCCESS;
}

