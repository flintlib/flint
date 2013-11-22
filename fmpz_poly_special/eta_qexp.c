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

    Copyright (C) 2011, 2013 Fredrik Johansson

******************************************************************************/

#include "fmpz_poly.h"

/* e = 3 :

prod(k>=1, 1-x^k )^3 = sum(n>=0, (-1)^n*(2*n+1)*x^(n*(n+1)/2) ) (jacobi)

Sparse formulas are used when e = 1, 2, 3, 4, 6.

Otherwise the computation is reduced to one of those cases,
followed by exponentiation.

note: P4 = P1 * P3

products?
P4 = P2 ^ 2
P5 = P2 * P3
P7 = P1 * P6
P8 = P2 * P6
P9 = P3 * P6
P12 = P6 ^ 2
P24 = P6 ^ 4


*/

static void
_eta_one(fmpz * f, slong len)
{
    slong k, n;

    _fmpz_vec_zero(f, len);
    f[0] = 1;

    for (n = k = 1; n + 4*k + 2 < len; k += 2)
    {
        f[n] = -1;
        f[n + k] = -1;
        f[n + 3*k + 1] = 1;
        f[n + 4*k + 2] = 1;
        n += 6*k + 5;
    }

    if (n < len) f[n] = -1;
    if (n + k < len) f[n + k] = -1;
    if (n + 3*k + 1 < len) f[n + 3*k + 1] = 1;
}

void
_eta_two(fmpz * c, slong N)
{
    slong k1, n1, k2, n2;
    int s, t;

    _fmpz_vec_zero(c, N);

#if 0
    for (k1=0, n1=0, s=1; n1 < N; n1+=3*k1+1, k1++, s=-s)
        for (k2=0, n2=0, t=s; n1 + n2 < N; n2+=3*k2+1, k2++, t=-t)
            c[n1 + n2] += t;

    for (k1=1, n1=2, s=-1; n1 < N; n1+=3*k1+2, k1++, s=-s)
        for (k2=1, n2=2, t=-s; n1 + n2 < N; n2+=3*k2+2, k2++, t=-t)
            c[n1 + n2] += t;
#endif

    /* TODO: SIGN ALWAYS -1 AT START OF IT? */

    /* P^2 */
#if 1
    for (k1=0, n1=0, s=1; 2 * n1 < N; n1+=3*k1+1, k1++)
        c[2 * n1] += s;
    for (k1=0, n1=0, s=2; n1 < N; n1+=3*k1+1, k1++, s=-s)
        for (k2=k1+1, n2=n1+3*k1+1, t=-2; n1 + n2 < N; n2+=3*k2+1, k2++, t=-t)
            c[n1 + n2] += t;

    /* Q^2 */
    for (k1=1, n1=2, s=1; 2 * n1 < N; n1+=3*k1+2, k1++)
        c[2 * n1] += s;
    for (k1=1, n1=2, s=-2; n1 < N; n1+=3*k1+2, k1++, s=-s)
        for (k2=k1+1, n2=n1+3*k1+2, t=-2; n1 + n2 < N; n2+=3*k2+2, k2++, t=-t)
            c[n1 + n2] += t;
#endif

    /* 2PQ */
    for (k1=0, n1=0, s=2; n1 < N; n1+=3*k1+1, k1++, s=-s)
        for (k2=1, n2=2, t=-s; n1 + n2 < N; n2+=3*k2+2, k2++, t=-t)
            c[n1 + n2] += t;
}

static void
_eta_three(fmpz * f, slong len)
{
    slong k, n;

    _fmpz_vec_zero(f, len);
    f[0] = 1;

    for (n = k = 1; k < len; k += n)
    {
        f[k] = (n % 2) ? (-(2*n+1)) : (2*n+1);
        n++;
    }
}

/* TODO: this can be done! (1 * 3 multiplication) */
static void
_eta_four(fmpz * f, long slen)
{

}

/* TODO: use squaring! */
static void
_eta_six(fmpz * f, slong len)
{
    slong j, k, jv, kv;
    fmpz_t tmp;

    _fmpz_vec_zero(f, len);
    fmpz_init(tmp);

    for (j = jv = 0; jv < len; jv += ++j)
    {
        fmpz_set_ui(tmp, 2*j + 1);

        for (k = kv = 0; jv + kv < len; kv += ++k)
        {
            if ((j+k) & 1)
                fmpz_submul_ui(f + jv + kv, tmp, 2*k+1);
            else
                fmpz_addmul_ui(f + jv + kv, tmp, 2*k+1);
        }
    }

    fmpz_clear(tmp);
}

void
_fmpz_poly_eta_qexp(fmpz * f, slong e, slong len)
{
    if (e == 1)
    {
        _eta_one(f, len);
    }
    else if (e == 2)
    {
        _eta_two(f, len);
    }
    else if (e == 3)
    {
        _eta_three(f, len);
    }
    else if (e == 4)
    {
        _eta_four(f, len);
    }
    else if (e == 6)
    {
        _eta_six(f, len);
    }
    else
    {
        fmpz *t, *u;

        t = u = _fmpz_vec_init(len);

        if (e < 0)
            t = f;

        else if (e % 6 == 0)
        {
            _eta_six(t, len);
            e /= 6;
        }
        else if (e % 3 == 0)
        {
            _eta_three(t, len);
            e /= 3;
        }
        else if (e % 4 == 0)
        {
            _eta_four(t, len);
            e /= 4;
        }
        else if (e % 2 == 0)
        {
            _eta_two(t, len);
            e /= 2;
        }
        else
        {
            _eta_one(t, len);
        }

        if (e < 0)
        {
            _fmpz_poly_pow_trunc(u, t, -e, len);
            _fmpz_poly_inv_series(f, u, len);
        }
        else if (e == 2)
        {
            _fmpz_poly_sqrlow(f, t, len, len);
        }
        else if (e == 4)
        {
            _fmpz_poly_sqrlow(f, t, len, len);
            _fmpz_poly_sqrlow(t, f, len, len);
            _fmpz_vec_swap(f, t, len);
        }
        else
        {
            _fmpz_poly_pow_trunc(f, t, e, len);
        }

        _fmpz_vec_clear(u, len);
    }
}

void
fmpz_poly_eta_qexp(fmpz_poly_t f, slong e, slong n)
{
    if (n < 1)
    {
        fmpz_poly_zero(f);
    }
    else if (e == 0 || n == 1)
    {
        fmpz_poly_one(f);
    }
    else
    {
        fmpz_poly_fit_length(f, n);
        _fmpz_poly_eta_qexp(f->coeffs, e, n);
        _fmpz_poly_set_length(f, n);
        _fmpz_poly_normalise(f);
    }
}
