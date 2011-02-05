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
#include "fmpz_poly.h"


int main()
{
    fmpq_poly_t P;
    fmpz_poly_t R, S;

    long i, n;

    printf("bernoulli_vec_2....");
    fflush(stdout);

    for (n = 0; n <= 3200; n += (n<100) ? 1 : n/3)
    {
        fmpq_poly_init2(P, n);
        fmpz_poly_init2(R, n);
        fmpz_poly_init2(S, n);

        fmpz_bernoulli_vec(P->den, P->coeffs, n);
        _fmpq_poly_set_length(P, n);

        fmpz_bernoulli_vec_2(R->coeffs, S->coeffs, n);
        _fmpz_poly_set_length(R, n);
        _fmpz_poly_set_length(S, n);

        for (i = 0; i < n; i++)
        {
            mpq_t x, y;
            mpq_init(x);
            mpq_init(y);

            fmpq_poly_get_coeff_mpq(x, P, i);
            fmpz_get_mpz(mpq_numref(y), R->coeffs + i);
            fmpz_get_mpz(mpq_denref(y), S->coeffs + i);

            if (!mpq_equal(x, y))
            {
                printf("ERROR: different results for %ld = \n", n);
                abort();
            }

            mpq_clear(y);
            mpq_clear(x);
        }


        fmpq_poly_clear(P);
        fmpz_poly_clear(R);
        fmpz_poly_clear(S);
    }

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
