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
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "arith.h"


void
_arith_cyclotomic_polynomial(fmpz * a, ulong n, mp_ptr factors,
                                        long num_factors, ulong phi)
{
    long i, k;
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
    if (factors[0] == 2UL)
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
            (n < 10163195L) ||                     /* At most 27 bits */
            (FLINT_BITS == 64 && n < 169828113L);  /* At most 60 bits */

    /* Iterate over all divisors of n */
    for (k = 0; k < (1L << num_factors); k++)
    {
        int mu;
        ulong d;

        mu = (num_factors & 1) ? -1 : 1;
        d = 1L;
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
    long i, j;
    ulong s, phi;

    if (n <= 2)
    {
        if (n == 0)
        {
            fmpz_poly_set_ui(poly, 1UL);
        }
        else
        {
            fmpz_poly_fit_length(poly, 2);
            fmpz_set_si(poly->coeffs, (n == 1) ? -1L : 1L);
            fmpz_set_si(poly->coeffs + 1, 1L);
            _fmpz_poly_set_length(poly, 2);
        }
        return;
    }

    /* Write n = q * s where q is squarefree, compute the factors of q,
      and compute phi(s) which determines the degree of the polynomial. */
    n_factor_init(&factors);
    n_factor(&factors, n, 1);
    s = phi = 1UL;
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
