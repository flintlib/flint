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

    Copyright (C) 2016  Ralf Stephan

******************************************************************************/

#include "fmpq_poly.h"

int main()
{
    fmpq_poly_t T0, T1, T2, t, tt;
    slong n;

    FLINT_TEST_INIT(state);

    flint_printf("legendre_pt....");
    fflush(stdout);

    fmpq_poly_init(T0);
    fmpq_poly_init(T1);
    fmpq_poly_init(T2);
    fmpq_poly_init(t);
    fmpq_poly_init(tt);

    fmpq_poly_legendre_p(T0, 0);
    fmpq_poly_legendre_p(T1, 1);

    for (n = 1; n <= 500; n++)
    {
        fmpq_poly_legendre_p(T2, n+1);
        fmpq_poly_set(t, T1);

        /* Verify (n+1)P_{n+1} = (2n+1)x P_n - nP_{n-1} */
        fmpq_poly_shift_left(t, t, 1);
        fmpq_poly_scalar_mul_ui(t, t, 2*n+1);
        fmpq_poly_scalar_mul_ui(tt, T0, n);
        fmpq_poly_sub(t, t, tt);
        fmpq_poly_scalar_mul_ui(tt, T2, n+1);

        if (!fmpq_poly_equal(t, tt))
        {
            flint_printf("FAIL: n = %wd\n", n);
            flint_printf("t: "); fmpq_poly_print_pretty(t, "x"); flint_printf("\n");
            flint_printf("tt: "); fmpq_poly_print_pretty(tt, "x"); flint_printf("\n");
            flint_printf("T2: "); fmpq_poly_print_pretty(T2, "x"); flint_printf("\n");
            abort();
        }

        fmpq_poly_swap(T0, T1);
        fmpq_poly_swap(T1, T2);
    }

    fmpq_poly_clear(T0);
    fmpq_poly_clear(T1);
    fmpq_poly_clear(T2);
    fmpq_poly_clear(t);
    fmpq_poly_clear(tt);

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
