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
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_vec.h"
#include "arith.h"
#include "ulong_extras.h"


static void
_rising_factorial(fmpz * res, long a, long b, long trunc)
{
    const long span = b - a;

    switch (span)
    {
    case 0:
        fmpz_one(res);
        break;
    case 1:
        fmpz_set_ui(res, a);
        if (trunc > 1) fmpz_one(res+1);
        break;
    case 2:
        fmpz_set_ui(res, a);
        fmpz_mul_ui(res, res, a + 1UL);
        if (trunc > 1)
        {
            fmpz_set_ui(res+1, 2*a + 1UL);
            if (trunc > 2) fmpz_one(res+2);
        }
        break;
    case 3:
        fmpz_set_ui(res, a);
        fmpz_mul_ui(res, res, a + 1UL);
        fmpz_mul_ui(res, res, a + 2UL);
        if (trunc > 1)
        {
            fmpz_set_ui(res+1, 3*a);
            fmpz_mul_ui(res+1, res+1, a + 2UL);
            fmpz_add_ui(res+1, res+1, 2);
            if (trunc > 2)
            {
                fmpz_set_ui(res+2, 3*(a+1));
                if (trunc > 3)
                    fmpz_one(res+3);
            }
        }
        break;
    default:
        {
            const long mid    = (a + b) / 2;
            const int  chk    = (b - a + 1 < trunc);  /* i.e. nprod < trunc */
            const long nleft  = chk ? mid - a + 1 : trunc;
            const long nright = chk ? b - mid + 1 : trunc;

            fmpz *left  = _fmpz_vec_init(nleft);
            fmpz *right = _fmpz_vec_init(nright);

            _rising_factorial(left, a, mid, trunc);
            _rising_factorial(right, mid, b, trunc);

            if (chk)
                _fmpz_poly_mul(res, right, nright, left, nleft);
            else
                _fmpz_poly_mullow(res, left, nleft, right, nright, trunc);

            _fmpz_vec_clear(left, nleft);
            _fmpz_vec_clear(right, nright);
        }
    }
}

void
arith_stirling_number_1u(fmpz_t s, long n, long k)
{
    /* Various special cases
       TODO: factorials, binomial coefficients, harmonic numbers ... */
    if (k < 1)
    {
        fmpz_set_ui(s, (n == 0) & (k == 0));
    }
    if (k >= n)
    {
        fmpz_set_ui(s, n == k);
    }
    else
    {
        fmpz *tmp = _fmpz_vec_init(k+1);
        _rising_factorial(tmp, 0, n, k+1);
        fmpz_set(s, tmp+k);
        _fmpz_vec_clear(tmp, k+1);
    }
}

void
arith_stirling_number_1(fmpz_t s, long n, long k)
{
    arith_stirling_number_1u(s, n, k);
    if ((n + k) % 2)
        fmpz_neg(s, s);
}

void
arith_stirling_number_1u_vec(fmpz * row, long n, long klen)
{
    if (klen > 0)
        _rising_factorial(row, 0, n, klen);
}

void
arith_stirling_number_1_vec(fmpz * row, long n, long klen)
{
    long k;

    arith_stirling_number_1u_vec(row, n, klen);

    for (k = (n + 1) % 2; k < klen; k += 2)
        fmpz_neg(row + k, row + k);
}
