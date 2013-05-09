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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz_poly.h"
#include "arith.h"
#include "ulong_extras.h"

void arith_divisors_naive(fmpz_poly_t p, len_t n)
{
    len_t k;
    len_t i = 0;

    n = FLINT_ABS(n);

    fmpz_poly_zero(p);
    for (k = 1; k <= n; k++)
    {
        if (n % k == 0)
        {
            fmpz_poly_set_coeff_si(p, i, k);
            i++;
        }
    }
}

int main(void)
{
    fmpz_t t;
    fmpz_poly_t a, b;
    len_t n;

    printf("divisors....");
    fflush(stdout);

    fmpz_init(t);
    fmpz_poly_init(a);
    fmpz_poly_init(b);

    for (n = -1000; n < 1000; n++)
    {
        fmpz_set_si(t, n);
        arith_divisors(a, t);
        arith_divisors_naive(b, n);
        if (!fmpz_poly_equal(a, b))
        {
            printf("FAIL:\n");
            printf("wrong value for n=%ld\n", n);
            abort();
        }
    }

    fmpz_clear(t);
    fmpz_poly_clear(a);
    fmpz_poly_clear(b);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
