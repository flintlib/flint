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
    fmpq_poly_t P, Q;

    long i, n;

    printf("bernoulli_vec....");
    fflush(stdout);

    for (n = 0; n <= 3200; n += (n<100) ? 1 : n/3)
    {
        fmpq_poly_init2(P, n);
        fmpq_poly_init2(Q, n);

        _fmpz_bernoulli_vec_recursive(P->den, P->coeffs, n);
        _fmpz_bernoulli_vec_series(Q->den, Q->coeffs, n);
        _fmpq_poly_set_length(P, n);
        _fmpq_poly_set_length(Q, n);

        if (!fmpq_poly_equal(P, Q))
        {
            mpq_t x, y;
            printf("ERROR: different results for %ld = \n", n);
            mpq_init(x);
            mpq_init(y);
            for (i = 0; i < n; i++)
            {
                fmpq_poly_get_coeff_mpq(x, P, i);
                fmpq_poly_get_coeff_mpq(y, Q, i);
                if (!mpq_equal(x, y))
                {
                    gmp_printf("P: %d: %Qd\n", i, x);
                    gmp_printf("Q: %d: %Qd\n", i, y);
                }
            }
            mpq_clear(y);
            mpq_clear(x);
            abort();
        }

        fmpq_poly_clear(P);
        fmpq_poly_clear(Q);
    }

    for (n = 0; n < 100; n++)
    {
        fmpq_poly_init2(P, n);
        fmpq_poly_init2(Q, n);

        _fmpz_bernoulli_vec_recursive(P->den, P->coeffs, n);
        fmpz_bernoulli_vec(Q->den, Q->coeffs, n);
        _fmpq_poly_set_length(P, n);
        _fmpq_poly_set_length(Q, n);

        if (!fmpq_poly_equal(P, Q))
        {
            printf("ERROR: different results for %ld = \n", n);
            abort();
        }

        fmpq_poly_clear(P);
        fmpq_poly_clear(Q);
    }

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
