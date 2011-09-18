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
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "arith.h"
#include "mpfr.h"
#include "math.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"


/* The coefficients in 2^d * \prod_{i=1}^d (x - cos(a_i)) are
   easily bounded using the binomial theorem. */
static long
magnitude_bound(long d)
{
    long res;
    fmpz_t t;
    fmpz_init(t);
    fmpz_bin_uiui(t, d, d / 2);
    res = fmpz_bits(t);
    fmpz_clear(t);
    return FLINT_ABS(res) + d;
}

static void
fmpz_mul_or_div_2exp(fmpz_t x, fmpz_t y, long s)
{
    if (s >= 0)
        fmpz_mul_2exp(x, y, s);
    else
        fmpz_fdiv_q_2exp(x, y, -s);
}


/* Balanced product of linear factors (x+alpha_i) using
   fixed-point arithmetic with prec bits */
static void
balanced_product(fmpz * c, fmpz * alpha, long len, long prec)
{
    if (len == 1)
    {
        fmpz_set_ui(c + 1, 1UL);
        fmpz_mul_2exp(c + 1, c + 1, prec);
        fmpz_set(c, alpha);
    }
    else if (len == 2)
    {
        fmpz_mul(c, alpha, alpha + 1);
        fmpz_fdiv_q_2exp(c, c, prec);
        fmpz_add(c + 1, alpha, alpha + 1);
        fmpz_set_ui(c + 2, 1UL);
        fmpz_mul_2exp(c + 2, c + 2, prec);
    }
    else
    {
        fmpz *L, *R;
        long i, m;

        m = len / 2;
        L = _fmpz_vec_init(len + 2);
        R = L + m + 1;

        balanced_product(L, alpha, m, prec);
        balanced_product(R, alpha + m, len - m, prec);
        _fmpz_poly_mul(c, R, len - m + 1, L, m + 1);

        for (i = 0; i < len + 1; i++)
            fmpz_fdiv_q_2exp(c + i, c + i, prec);

        _fmpz_vec_clear(L, len + 2);
    }
}

void
_cyclotomic_cos_polynomial(fmpz * coeffs, long d, ulong n)
{
    long i, j, prec, exp;
    fmpz * alpha;
    fmpz_t half;
    mpfr_t t, u;

    prec = magnitude_bound(d) + 5 + FLINT_BIT_COUNT(d);

    alpha = _fmpz_vec_init(d);
    fmpz_init(half);
    mpfr_init2(t, prec);
    mpfr_init2(u, prec);

    fmpz_set_ui(half, 1UL);
    fmpz_mul_2exp(half, half, prec - 1);
    mpfr_const_pi(t, prec);
    mpfr_div_ui(t, t, n, MPFR_RNDN);

    for (i = j = 0; j < d; i++)
    {
        if (n_gcd(n, i) == 1)
        {
            mpfr_mul_ui(u, t, 2 * i, MPFR_RNDN);
            mpfr_cos(u, u, MPFR_RNDN);
            mpfr_neg(u, u, MPFR_RNDN);
            exp = mpfr_get_z_2exp(_fmpz_promote(alpha + j), u);
            _fmpz_demote_val(alpha + j);
            fmpz_mul_or_div_2exp(alpha + j, alpha + j, exp + prec);
            j++;
        }
    }

    balanced_product(coeffs, alpha, d, prec);

    /* Scale and round */
    for (i = 0; i < d + 1; i++)
    {
        long r = d;
        if ((n & (n - 1)) == 0)
            r--;
        fmpz_mul_2exp(coeffs + i, coeffs + i, r);
        fmpz_add(coeffs + i, coeffs + i, half);
        fmpz_fdiv_q_2exp(coeffs + i, coeffs + i, prec);
    }

    fmpz_clear(half);
    mpfr_clear(t);
    mpfr_clear(u);
    _fmpz_vec_clear(alpha, d);
}

void
cyclotomic_cos_polynomial(fmpz_poly_t poly, ulong n)
{
    if (n == 0)
    {
        fmpz_poly_set_ui(poly, 1UL);
    }
    else
    {
        long d = (n <= 2) ? 1 : n_euler_phi(n) / 2;

        fmpz_poly_fit_length(poly, d + 1);
        _cyclotomic_cos_polynomial(poly->coeffs, d, n);
        _fmpz_poly_set_length(poly, d + 1);
    }
}
