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
#include "arith.h"
#include "profiler.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq_poly.h"


int main()
{
    fmpz * nums;
    fmpz * dens;
    fmpz_t num;
    fmpz_t den;
    long n, N;

    printf("bernoulli_number_zeta....");
    fflush(stdout);

    N = 5000;

    nums = _fmpz_vec_init(N);
    dens = _fmpz_vec_init(N);
    fmpz_bernoulli_vec_2(nums, dens, N);

    for (n = 0; n < N; n++)
    {
        fmpz_init(num);
        fmpz_init(den);

        bernoulli_number(num, den, n);
        if (!fmpz_equal(num, nums + n) || !fmpz_equal(den, dens + n))
        {
            printf("FAIL: n = %ld\n", n);
            printf("Numerator: "); fmpz_print(num); printf("\n");
            printf("Denominator: "); fmpz_print(den); printf("\n");
            abort();
        }

        fmpz_clear(num);
        fmpz_clear(den);
    }

    _fmpz_vec_clear(nums, N);
    _fmpz_vec_clear(dens, N);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
