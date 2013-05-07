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
#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "arith.h"
#include "profiler.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq_poly.h"


int main()
{
    fmpz * ress;
    fmpz_t res;
    long n, N;

    printf("euler_number_zeta....");
    fflush(stdout);

    N = 3000;

    ress = _fmpz_vec_init(N);
    arith_euler_number_vec(ress, N);

    for (n = 0; n < N; n++)
    {
        fmpz_init(res);

        arith_euler_number(res, n);
        if (!fmpz_equal(res, ress + n))
        {
            printf("FAIL: n = %ld\n", n);
            printf("Value: "); fmpz_print(res); printf("\n");
            abort();
        }

        fmpz_clear(res);
    }

    _fmpz_vec_clear(ress, N);

    mpfr_free_cache();
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
