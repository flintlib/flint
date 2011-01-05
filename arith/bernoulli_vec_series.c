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
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "fmpz_vec.h"
#include "arith.h"


void _bernoulli_vec_series_buhler(fmpq_poly_t B[4], long n)
{
    long i, j, k, l;
    fmpz_t t;

    fmpq_poly_t A[5];

    for (i = 0; i < 5; i++)
        fmpq_poly_init2(A[i], n);

    fmpz_init(t);

    fmpq_poly_set_coeff_si(A[0], 0, 6L);
    fmpq_poly_set_coeff_si(A[0], 1, -792L);
    fmpq_poly_set_coeff_si(A[1], 0, 20L);
    fmpq_poly_set_coeff_si(A[1], 1, -2704L);
    fmpq_poly_set_coeff_si(A[2], 0, -28L);
    fmpq_poly_set_coeff_si(A[2], 1, 3824L);
    fmpq_poly_set_coeff_si(A[3], 0, 96L);
    fmpq_poly_set_coeff_si(A[3], 1, -13056L);
    fmpq_poly_set_coeff_si(A[4], 0, 24L);
    fmpq_poly_set_coeff_si(A[4], 1, -3168L);

    /* Two-term recurrence */
    for (i = 0; i < 5; i++)
    {
        for (k = 2; k < n; k++)
        {
            fmpz_mul_si(A[i]->coeffs + k, A[i]->coeffs + k - 1, -136L);
            fmpz_submul_ui(A[i]->coeffs + k, A[i]->coeffs + k - 2, 16L);
        }
        _fmpq_poly_set_length(A[i], n);
    }

    /* Divide by factorials */
    for (i = 0; i < 5; i++)
    {
        fmpz_set_ui(t, 1UL);
        for (k = n-2; k >= 0; k--)
        {
            if (i == 4)
                j = 1;
            else
                j = 2*i;
            for (l = 4; l <= 11; l++)
                fmpz_mul_ui(t, t, 8*k + l + j);
            fmpz_mul(A[i]->coeffs + k, A[i]->coeffs + k, t);
        }
        if (i == 4)
            fmpz_mul_ui(t, t, 1*2*3*4);
        else
            for (k = 1; k <= 2*i+3; k++)
                fmpz_mul_ui(t, t, k);

        fmpz_set(fmpq_poly_denref(A[i]), t);
    }

    fmpq_poly_inv_series(A[4], A[4], n);

    for (i = 0; i < 4; i++)
    {
        fmpq_poly_mullow(B[i], A[i], A[4], n);
        fmpq_poly_clear(A[i]);
    }

    fmpq_poly_clear(A[4]);
    fmpz_clear(t);
}

void _fmpz_bernoulli_vec_series(fmpz_t den, fmpz * b, long n)
{
    fmpq_poly_t B[4];
    fmpz_t ratio_p[4];
    fmpz_t ratio_q[4];
    fmpz_t t;

    long i, k, m;

    m = FLINT_MAX(2, (n + 8 - 1) / 8);

    fmpz_primorial(den, n + 1);

    for (i = 0; i < 4; i++)
        fmpq_poly_init2(B[i], m);

    _bernoulli_vec_series_buhler(B, m);

    /* Set odd entries */
    if (n > 1)
        fmpz_divexact_si(b + 1, den, -2L);
    for (k = 3; k < n; k += 2)
        fmpz_zero(b + k);

    fmpz_init(t);

    /* Rescale */
    for (i = 0; i < 4; i++)
    {
        fmpz_init(ratio_p[i]);
        fmpz_init(ratio_q[i]);
        fmpz_gcd(t, B[i]->den, den);
        fmpz_divexact(ratio_p[i], B[i]->den, t);
        fmpz_divexact(ratio_q[i], den, t);
    }

    fmpz_set_ui(t, 1UL);
    for (k = 0; k < n; k += 2)
    {
        i = (k/2) % 4;
        fmpz_mul(b + k, B[i]->coeffs + (k/8), t);
        fmpz_mul(b + k, b + k, ratio_q[i]);
        fmpz_divexact(b + k, b + k, ratio_p[i]);
        fmpz_mul_ui(t, t, (k+1)*(k+2)/2);
    }

    fmpz_clear(t);

    for (i = 0; i < 4; i++)
    {
        fmpq_poly_clear(B[i]);
        fmpz_clear(ratio_p[i]);
        fmpz_clear(ratio_q[i]);
    }
}
