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
#include "fmpz.h"
#include "arith.h"
#include "ulong_extras.h"

void fmpz_sigma_naive(fmpz_t x, ulong n, ulong k)
{
    len_t i = 0;

    fmpz_t t;
    fmpz_poly_t p;
    fmpz_init(t);
    fmpz_poly_init(p);
    fmpz_set_ui(t, n);
    arith_divisors(p, t);

    fmpz_zero(x);
    for (i = 0; i < p->length; i++)
    {
        fmpz_poly_get_coeff_fmpz(t, p, i);
        fmpz_pow_ui(t, t, k);
        fmpz_add(x, x, t);
    }

    fmpz_clear(t);
    fmpz_poly_clear(p);
}

int main(void)
{
    fmpz_t m, a, b;
    len_t n, k;

    printf("divisor_sigma....");
    fflush(stdout);

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(m);

    for (n = 0; n < 5000; n++)
    {
        for (k = 0; k < 10; k++)
        {
            fmpz_set_ui(m, n);
            arith_divisor_sigma(a, m, k);
            fmpz_sigma_naive(b, n, k);
            if (!fmpz_equal(a, b))
            {
                printf("FAIL:\n");
                printf("wrong value for n=%ld, k=%ld\n", n, k);
                abort();
            }
        }
    }

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(m);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
