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

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_factor.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "arith.h"

static void
sum_of_two_squares(fmpz_t r, const fmpz_t n)
{
    fmpz_factor_t fac;
    len_t i;

    fmpz_factor_init(fac);
    fmpz_factor(fac, n);
    fmpz_one(r);

    for (i = 0; i < fac->num; i++)
    {
        const int res = fmpz_fdiv_ui(fac->p + i, 4);

        if (res == 1)
        {
            fac->exp[i]++;
            fmpz_mul_ui(r, r, fac->exp[i]);
        }
        else if (res == 3)
        {
            if (fac->exp[i] % 2)
            {
                fmpz_zero(r);
                break;
            }
        }
    }

    fmpz_mul_ui(r, r, 4);
    fmpz_factor_clear(fac);
}

static void
sum_of_four_squares(fmpz_t r, const fmpz_t n)
{
    const mp_bitcnt_t v = fmpz_val2(n);

    if (v == 0)
    {
        arith_divisor_sigma(r, n, 1);
        fmpz_mul_ui(r, r, 8);
    }
    else
    {
        fmpz_tdiv_q_2exp(r, n, v);
        arith_divisor_sigma(r, r, 1);
        fmpz_mul_ui(r, r, 24);
    }
}

static void
sum_of_squares_recursive(fmpz_t r, len_t k, ulong n)
{
    fmpz_t t, u;
    len_t i, j;

    fmpz_init(t);
    fmpz_init(u);
    fmpz_zero(r);

    for (i = j = 0; j <= n; i++)
    {
        fmpz_set_ui(u, n - j);
        arith_sum_of_squares(t, k - 1, u);

        if (j > 0)
            fmpz_mul_ui(t, t, 2);
        fmpz_add(r, r, t);

        j += 2 * i + 1;
    }

    fmpz_clear(t);
    fmpz_clear(u);
}

static void
sum_of_squares_series(fmpz_t r, ulong k, len_t n)
{
    fmpz * t;

    t = _fmpz_vec_init(n + 1);
    arith_sum_of_squares_vec(t, k, n + 1);
    fmpz_set(r, t + n);
    _fmpz_vec_clear(t, n + 1);
}

void
arith_sum_of_squares(fmpz_t r, ulong k, const fmpz_t n)
{
    if (fmpz_sgn(n) <= 0 || k == 0)
        fmpz_set_ui(r, fmpz_is_zero(n) != 0);
    else if (k == 1)
        fmpz_set_ui(r, 2 * (fmpz_is_square(n) != 0));
    else if (k == 2)
        sum_of_two_squares(r, n);
    else if (k == 4)
        sum_of_four_squares(r, n);
    else if (k == 3 || k == 5)
        sum_of_squares_recursive(r, k, fmpz_get_ui(n));
    else if (fmpz_fits_si(n))
        sum_of_squares_series(r, k, fmpz_get_ui(n));
    else
    {
        printf("Exception (arith_sum_of_squares). n is too large.\n");
        abort();
    }
}
