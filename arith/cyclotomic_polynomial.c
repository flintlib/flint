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

#include "arith.h"

void
_arith_cyclotomic_polynomial(fmpz * a, ulong n, mp_ptr factors,
                                        slong num_factors, ulong phi)
{
    slong i, k;
    int small;
    ulong D;

    D = phi / 2;

    /* Phi_p(x) = 1 + x + x^2 + ... + x^{p-1} */
    if (num_factors == 1)
    {
        for (i = 0; i <= D; i++)
            fmpz_one(a + i);
        return;
    }

    /* Phi_{2n}(x) = Phi_n(-x)*/
    if (factors[0] == UWORD(2))
    {
        _arith_cyclotomic_polynomial(a, n / 2, factors + 1,
            num_factors - 1, phi);
        for (i = 1; i <= D; i += 2)
            fmpz_neg(a + i, a + i);
        return;
    }

    fmpz_one(a);
    for (i = 1; i <= D; i++)
        fmpz_zero(a + i);

    /* Coefficients are guaranteed not to overflow an fmpz */
    small = (num_factors == 2) ||                  /* Always +1/0/-1*/
            (n < WORD(10163195)) ||                     /* At most 27 bits */
            (FLINT_BITS == 64 && n < WORD(169828113));  /* At most 60 bits */

    /* Iterate over all divisors of n */
    for (k = 0; k < (WORD(1) << num_factors); k++)
    {
        int mu;
        ulong d;

        mu = (num_factors & 1) ? -1 : 1;
        d = WORD(1);
        for (i = 0; i < num_factors; i++)
        {
            if ((k >> i) & 1)
            {
                d *= factors[i];
                mu = -mu;
            }
        }

        /* Multiply by (x^d - 1)^{\mu(n/d)} */
        if (small)
        {
            if (mu == 1)
                for (i = D; i >= d; i--) a[i] -= a[i - d];
            else
                for (i = d; i <= D; i++) a[i] += a[i - d];
        }
        else
        {
            if (mu == 1)
                for (i = D; i >= d; i--) fmpz_sub(a + i, a + i, a + i - d);
            else
                for (i = d; i <= D; i++) fmpz_add(a + i, a + i, a + i - d);
        }
    }
}

void
arith_cyclotomic_polynomial(fmpz_poly_t poly, ulong n)
{
    n_factor_t factors;
    slong i, j;
    ulong s, phi;

    if (n <= 2)
    {
        if (n == 0)
        {
            fmpz_poly_set_ui(poly, UWORD(1));
        }
        else
        {
            fmpz_poly_fit_length(poly, 2);
            fmpz_set_si(poly->coeffs, (n == 1) ? WORD(-1) : WORD(1));
            fmpz_set_si(poly->coeffs + 1, WORD(1));
            _fmpz_poly_set_length(poly, 2);
        }
        return;
    }

    /* Write n = q * s where q is squarefree, compute the factors of q,
      and compute phi(s) which determines the degree of the polynomial. */
    n_factor_init(&factors);
    n_factor(&factors, n, 1);
    s = phi = UWORD(1);
    for (i = 0; i < factors.num; i++)
    {
        phi *= factors.p[i] - 1;
        while (factors.exp[i] > 1)
        {
            s *= factors.p[i];
            factors.exp[i]--;
        }
    }

    fmpz_poly_fit_length(poly, phi * s + 1);

    /* Evaluate lower half of Phi_s(x) */
    _arith_cyclotomic_polynomial(poly->coeffs, n / s,
        factors.p, factors.num, phi);

    /* Palindromic extension */
    for (i = 0; i < (phi + 1) / 2; i++)
        fmpz_set(poly->coeffs + phi - i, poly->coeffs + i);

    /* Stretch */
    if (s != 1)
    {
        for (i = phi; i > 0; i--)
        {
            fmpz_set(poly->coeffs + i*s, poly->coeffs + i);
            for (j = 1; j < s; j++)
                fmpz_zero(poly->coeffs + i*s - j);
        }
    }

    _fmpz_poly_set_length(poly, phi * s + 1);
}
