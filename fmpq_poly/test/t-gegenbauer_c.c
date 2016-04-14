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
    fmpq_t a, rat;
    slong n, d;

    FLINT_TEST_INIT(state);

    flint_printf("gegenbauer_c....");
    fflush(stdout);

    fmpq_poly_init(T0);
    fmpq_poly_init(T1);
    fmpq_poly_init(t);
    fmpq_poly_init(tt);
    fmpq_init(a);
    fmpq_init(rat);

    for (d = 1; d < 11; d++)
    {
        fmpq_set_si(a, 1, d);
        fmpq_poly_gegenbauer_c(T0, 0, a);
        fmpq_poly_gegenbauer_c(T1, 1, a);

        for (n = 1; n <= 500; n++)
        {
            fmpq_poly_init(T2);
            fmpq_poly_gegenbauer_c(T2, n+1, a);
            fmpq_poly_set(t, T1);

            /* Verify (n+1)C^a_{n+1} = 2x(n+a) C^a_n - (n+2a-1)C^a_{n-1} */
            fmpq_poly_shift_left(t, t, 1);
            fmpq_set(rat, a);
            fmpq_add_si(rat, rat, n);
            fmpq_mul_2exp(rat, rat, 1);
            fmpq_poly_scalar_mul_fmpq(t, t, rat);

            fmpq_set(rat, a);
            fmpq_mul_2exp(rat, rat, 1);
            fmpq_add_si(rat, rat, n-1);
            fmpq_poly_scalar_mul_fmpq(tt, T0, rat);
            fmpq_poly_sub(t, t, tt);
            fmpq_poly_scalar_mul_si(tt, T2, n+1);

            fmpq_poly_canonicalise(t);
            fmpq_poly_canonicalise(tt);
            if (!fmpq_poly_equal(t, tt))
            {
                flint_printf("\nFAIL: n = %wd, a = ", n);
                fmpq_print(a); flint_printf("\n");
                flint_printf("t: "); fmpq_poly_print_pretty(t, "x"); flint_printf("\n");
                flint_printf("tt: "); fmpq_poly_print_pretty(tt, "x"); flint_printf("\n");
                abort();
            }

            fmpq_poly_clear(T0);
            fmpq_poly_init(T0);
            fmpq_poly_swap(T0, T1);
            fmpq_poly_clear(T1);
            fmpq_poly_init(T1);
            fmpq_poly_swap(T1, T2);
            fmpq_poly_clear(T2);
        }
    }

    fmpq_poly_clear(T0);
    fmpq_poly_clear(T1);
    fmpq_poly_clear(T2);
    fmpq_poly_clear(t);
    fmpq_poly_clear(tt);
    fmpq_clear(a);
    fmpq_clear(rat);

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
