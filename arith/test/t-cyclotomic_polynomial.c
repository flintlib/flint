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
#include "flint.h"
#include "arith.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

void cyclotomic_naive(fmpz_poly_t poly, ulong n)
{
    fmpz_poly_t t;
    long d;

    fmpz_poly_init(t);

    fmpz_poly_set_ui(poly, 1UL);
    for (d = 1; d <= n; d++)
    {
        if (n % d == 0)
        {
            if (n_moebius_mu(n / d) == 1)
            {
                fmpz_poly_zero(t);
                fmpz_poly_set_coeff_si(t, d, 1);
                fmpz_poly_set_coeff_si(t, 0, -1);
                fmpz_poly_mul(poly, poly, t);
            }
        }
    }

    for (d = 1; d <= n; d++)
    {
        if (n % d == 0)
        {
            if (n_moebius_mu(n / d) == -1)
            {
                fmpz_poly_zero(t);
                fmpz_poly_set_coeff_si(t, d, 1);
                fmpz_poly_set_coeff_si(t, 0, -1);
                fmpz_poly_div(poly, poly, t);
            }
        }
    }

    fmpz_poly_clear(t);
}

int main()
{
    fmpz_poly_t A, B;
    long n;

    printf("cyclotomic_polynomial....");
    fflush(stdout);

    for (n = 0; n <= 1000; n++)
    {
        fmpz_poly_init(A);
        fmpz_poly_init(B);

        arith_cyclotomic_polynomial(A, n);
        cyclotomic_naive(B, n);

        if (!fmpz_poly_equal(A, B))
        {
            printf("FAIL: wrong value of Phi_%ld(x)\n", n);
            printf("Computed:\n");
            fmpz_poly_print_pretty(A, "x");
            printf("\n\nExpected:\n");
            fmpz_poly_print_pretty(B, "x");
            printf("\n\n");
            abort();
        }

        fmpz_poly_clear(A);
        fmpz_poly_clear(B);
    }

    /* We verify the first value that does not fit on 32 bits.
       This exercises the slow path at least on a 32 bit system.
       Testing the 64 bit value is a bit too much to do by default
        as it requires ~2 GB of memory and takes a few minutes. */
    {
        fmpz_t h, ref;

        const ulong nn = 10163195UL;
        /* const ulong nn = 169828113UL;  64-bit case */

        fmpz_init(h);
        fmpz_init(ref);
        fmpz_set_str(ref, "1376877780831", 10);
        /* fmpz_set_str(ref, "31484567640915734941", 10);  64-bit case */

        fmpz_poly_init(A);
        arith_cyclotomic_polynomial(A, 10163195UL);
        fmpz_poly_height(h, A);

        if (!fmpz_equal(h, ref))
        {
            printf("Bad computation of Phi_%ld(x)\n", nn);
            printf("Computed height:\n");
            fmpz_print(h);
            printf("\nExpected height:\n");
            fmpz_print(ref);
            printf("\n\n");
            abort();
        }

        fmpz_poly_clear(A);
        fmpz_clear(h);
        fmpz_clear(ref);
    }

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
