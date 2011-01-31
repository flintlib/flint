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
#include <limits.h>
#include <mpir.h>
#include "flint.h"
#include "arith.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"
#include "profiler.h"

int main(void)
{
    fmpz * b1;
    fmpz * b2;
    long n, k;

    const long maxn = 400;

    printf("bell....");
    fflush(stdout);

    b1 = _fmpz_vec_init(maxn);

    /* Consistency test */
    for (n = 0; n < maxn; n++)
        fmpz_bell(b1 + n, n);

    for (n = 0; n < maxn; n++)
    {
        b2 = _fmpz_vec_init(n);
        fmpz_bell_vec(b2, n);

        if (!_fmpz_vec_equal(b1, b2, n))
        {
            printf("FAIL:\n");
            printf("n = %ld\n", n);
            abort();
        }

        _fmpz_vec_clear(b2, n);
    }

    /* Compare with B_n = sum of Stirling numbers of 2nd kind */
    for (n = 0; n < 2500; n += (n < 50) ? + 1 : n/4)
    {
        b2 = _fmpz_vec_init(n+1);

        fmpz_stirling2_vec(b2, n, n+1);

        for (k = 1; k <= n; k++)
            fmpz_add(b2, b2, b2 + k);

        fmpz_bell(b1, n);

        if (!fmpz_equal(b1, b2))
        {
            printf("FAIL:\n");
            printf("n = %ld\n", n);
            fmpz_print(b1);
            printf("\n");
            fmpz_print(b2);
            printf("\n");
            abort();
        }

        _fmpz_vec_clear(b2, n+1);
    }

    _fmpz_vec_clear(b1, maxn);

    printf("PASS\n");
    return 0;
}
