/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "fmpz_poly.h"

/* Bound coefficients using (x + u)^(2^n) and the binomial
   coefficients. TODO: this is about 2x too large... */
static slong __bound_prec(ulong n)
{
    slong i;
    double u, N;

    N = UWORD(1) << n;

    /* u = (sum of square roots)^(2^n) */
    u = 0;
    for (i = 0; i < n; i++)
        u += sqrt(n_nth_prime(1 + i));
    u = N * log(u) * 1.44269504088897;

    /* Central binomial coefficient C(N,N/2) < 2^N / sqrt(3*N/2) */
    u += N - 0.5*(n-1) - 0.792481250360578;  /* log(sqrt(3)) */

    return u;
}

void
_fmpz_poly_swinnerton_dyer(fmpz * T, ulong n)
{
    fmpz *square_roots, *tmp1, *tmp2, *tmp3;
    fmpz_t one;
    slong i, j, k, N;
    slong prec;

    if (n == 0)
    {
        fmpz_zero(T);
        fmpz_one(T + 1);
        return;
    }

    N = WORD(1) << n;
    prec = __bound_prec(n);

    fmpz_init(one);
    fmpz_one(one);
    fmpz_mul_2exp(one, one, prec);

    square_roots = _fmpz_vec_init(n);
    tmp1 = flint_malloc((N/2 + 1) * sizeof(fmpz));
    tmp2 = flint_malloc((N/2 + 1) * sizeof(fmpz));
    tmp3 = _fmpz_vec_init(N);

    for (i = 0; i < n; i++)
    {
        fmpz_set_ui(square_roots + i, n_nth_prime(i + 1));
        fmpz_mul_2exp(square_roots + i, square_roots + i, 2 * prec);
        fmpz_sqrt(square_roots + i, square_roots + i);
    }

    /* Build linear factors */
    for (i = 0; i < N; i++)
    {
        fmpz_zero(T + i);
        for (j = 0; j < n; j++)
        {
            if ((i >> j) & 1)
                fmpz_add(T + i, T + i, square_roots + j);
            else
                fmpz_sub(T + i, T + i, square_roots + j);
        }
    }

    /* For each level... */
    for (i = 0; i < n; i++)
    {
        slong stride = UWORD(1) << i;

        for (j = 0; j < N; j += 2*stride)
        {
            for (k = 0; k < stride; k++)
            {
                tmp1[k] = T[j + k];
                tmp2[k] = T[j + stride + k];
            }
            tmp1[stride] = *one;
            tmp2[stride] = *one;

            _fmpz_poly_mullow(tmp3, tmp1, stride + 1, tmp2, stride + 1, 2*stride);
            _fmpz_vec_scalar_fdiv_q_2exp(T + j, tmp3, 2*stride, prec);
        }
    }

    /* Round */
    fmpz_fdiv_q_2exp(one, one, 1);
    for (i = 0; i < N; i++)
        fmpz_add(T + i, T + i, one);

    _fmpz_vec_scalar_fdiv_q_2exp(T, T, N, prec);
    fmpz_one(T + (UWORD(1) << n));

    _fmpz_vec_clear(square_roots, n);
    flint_free(tmp1);
    flint_free(tmp2);
    _fmpz_vec_clear(tmp3, UWORD(1) << n);
    fmpz_clear(one);
}

void
fmpz_poly_swinnerton_dyer(fmpz_poly_t poly, ulong n)
{
    slong N = (WORD(1) << n);
    fmpz_poly_fit_length(poly, N + 1);
    _fmpz_poly_swinnerton_dyer(poly->coeffs, n);
    _fmpz_poly_set_length(poly, N + 1);
}

