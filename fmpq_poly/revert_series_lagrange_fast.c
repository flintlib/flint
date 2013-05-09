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

   Copyright (C) 2010 Sebastian Pancratz
   Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "ulong_extras.h"


/* pointer to (x/Q)^i */
#define Ri(ii) (R + (n-1)*((ii)-1))
#define Rdeni(ii) (Rden + ii - 1)

static void
_set_vec(fmpz * rnum, fmpz_t den,
                const fmpz * xnum, const fmpz * xden, len_t len)
{
    len_t j;
    fmpz_t t;
    fmpz_init(t);
    fmpz_one(den);

    for (j = 0; j < len; j++)
        fmpz_lcm(den, den, xden + j);

    for (j = 0; j < len; j++)
    {
        fmpz_divexact(t, den, xden + j);
        fmpz_mul(rnum + j, xnum + j, t);
    }

    fmpz_clear(t);
}

void
_fmpq_poly_revert_series_lagrange_fast(fmpz * Qinv, fmpz_t den,
                                    const fmpz * Q, const fmpz_t Qden, len_t n)
{
    len_t i, j, k, m;
    fmpz *R, *Rden, *S, *T, *dens, *tmp;
    fmpz_t Sden, Tden, t;

    if (fmpz_is_one(Qden) && (n > 1) && fmpz_is_pm1(Q + 1))
    {
        _fmpz_poly_revert_series(Qinv, Q, n);
        fmpz_one(den);
        return;
    }

    if (n <= 2)
    {
        fmpz_zero(Qinv);
        if (n == 2)
        {
            fmpz_set(Qinv + 1, Qden);
            fmpz_set(den, Q + 1);
            _fmpq_poly_canonicalise(Qinv, den, 2);
        }
        return;
    }

    m = n_sqrt(n);

    fmpz_init(t);
    dens = _fmpz_vec_init(n);
    R = _fmpz_vec_init((n - 1) * m);
    S = _fmpz_vec_init(n - 1);
    T = _fmpz_vec_init(n - 1);
    Rden = _fmpz_vec_init(m);
    fmpz_init(Sden);
    fmpz_init(Tden);

    fmpz_zero(Qinv);
    fmpz_one(dens);

    _fmpq_poly_inv_series(Ri(1), Rdeni(1), Q + 1, Qden, n - 1);
    _fmpq_poly_canonicalise(Ri(1), Rdeni(1), n - 1);

    for (i = 2; i <= m; i++)
    {
        _fmpq_poly_mullow(Ri(i), Rdeni(i), Ri(i-1), Rdeni(i-1), n - 1,
                Ri(1), Rdeni(1), n - 1, n - 1);
        _fmpq_poly_canonicalise(Ri(i), Rdeni(i), n - 1);
    }

    for (i = 1; i < m; i++)
    {
        fmpz_set(Qinv + i, Ri(i) + i - 1);
        fmpz_mul_ui(dens + i, Rdeni(i), i);
    }

    _fmpz_vec_set(S, Ri(m), n - 1);
    fmpz_set(Sden, Rdeni(m));

    for (i = m; i < n; i += m)
    {
        fmpz_set(Qinv + i, S + i - 1);
        fmpz_mul_ui(dens + i, Sden, i);

        for (j = 1; j < m && i + j < n; j++)
        {
            fmpz_mul(t, S + 0, Ri(j) + i + j - 1);

            for (k = 1; k <= i + j - 1; k++)
                fmpz_addmul(t, S + k, Ri(j) + i + j - 1 - k);

            fmpz_set(Qinv + i + j, t);
            fmpz_mul(dens + i + j, Sden, Rdeni(j));
            fmpz_mul_ui(dens + i + j, dens + i + j, i + j);
        }

        if (i + 1 < n)
        {
            _fmpq_poly_mullow(T, Tden, S, Sden, n - 1,
                Ri(m), Rdeni(m), n - 1, n - 1);
            _fmpq_poly_canonicalise(T, Tden, n - 1);
            fmpz_swap(Tden, Sden);
            tmp = S; S = T; T = tmp;
        }
    }

    _set_vec(Qinv, den, Qinv, dens, n);
    _fmpq_poly_canonicalise(Qinv, den, n);

    fmpz_clear(t);
    _fmpz_vec_clear(dens, n);
    _fmpz_vec_clear(R, (n - 1) * m);
    _fmpz_vec_clear(S, n - 1);
    _fmpz_vec_clear(T, n - 1);
    _fmpz_vec_clear(Rden, m);
    fmpz_clear(Sden);
    fmpz_clear(Tden);
}

void
fmpq_poly_revert_series_lagrange_fast(fmpq_poly_t res,
            const fmpq_poly_t poly, len_t n)
{
    fmpz *copy;
    int alloc;

    if (poly->length < 2 || !fmpz_is_zero(poly->coeffs)
                         || fmpz_is_zero(poly->coeffs + 1))
    {
        printf("Exception (fmpq_poly_revert_series_lagrange_fast). Input must \n"
               "have zero constant term and nonzero coefficient of x^1.\n");
        abort();
    }

    if (n < 2)
    {
        fmpq_poly_zero(res);
        return;
    }

    if (poly->length >= n)
    {
        copy = poly->coeffs;
        alloc = 0;
    }
    else
    {
        len_t i;
        copy = (fmpz *) flint_malloc(n * sizeof(fmpz));
        for (i = 0; i < poly->length; i++)
            copy[i] = poly->coeffs[i];
        for ( ; i < n; i++)
            copy[i] = 0;
        alloc = 1;
    }

    if (res != poly)
    {
        fmpq_poly_fit_length(res, n);
        _fmpq_poly_revert_series_lagrange_fast(res->coeffs,
                res->den, copy, poly->den, n);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        _fmpq_poly_revert_series_lagrange_fast(t->coeffs,
                t->den, copy, poly->den, n);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }

    _fmpq_poly_set_length(res, n);
    _fmpq_poly_normalise(res);

    if (alloc)
        flint_free(copy);
}
