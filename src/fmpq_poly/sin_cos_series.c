/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void
_fmpq_poly_sin_cos_series_basecase_can(fmpz * S, fmpz_t Sden,
    fmpz * C, fmpz_t Cden, const fmpz * A, const fmpz_t Aden, slong Alen, slong n, int can)
{
    fmpz_t t, u, v;
    slong j, k;

    Alen = FLINT_MIN(Alen, n);

    if (Alen == 1 || n == 1)
    {
        fmpz_zero(S);
        fmpz_one(C);
        _fmpz_vec_zero(S + 1, n - 1);
        _fmpz_vec_zero(C + 1, n - 1);
        fmpz_one(Sden);
        fmpz_one(Cden);
        return;
    }

    /* support aliasing */
    if (A == S || A == C)
    {
        fmpz * tmp = _fmpz_vec_init(Alen + 1);
        _fmpz_vec_set(tmp, A, Alen);
        fmpz_set(tmp + Alen, Aden);
        _fmpq_poly_sin_cos_series_basecase_can(S, Sden, C, Cden,
            tmp, tmp + Alen, Alen, n, can);
        _fmpz_vec_clear(tmp, Alen + 1);
        return;
    }

    fmpz_init(t);
    fmpz_init(u);
    fmpz_init(v);

    fmpz_fac_ui(t, n - 1);
    fmpz_pow_ui(v, Aden, n - 1);
    fmpz_mul(Sden, t, v);
    fmpz_set(C, Sden);
    fmpz_set(Cden, Sden);
    fmpz_zero(S);

    for (k = 1; k < n; k++)
    {
        fmpz_zero(t);
        fmpz_zero(u);

        for (j = 1; j < FLINT_MIN(Alen, k + 1); j++)
        {
            fmpz_mul_ui(v, A + j, j);
            fmpz_submul(t, v, S + k - j);
            fmpz_addmul(u, v, C + k - j);
        }

        fmpz_mul_ui(v, Aden, k);
        fmpz_divexact(C + k, t, v);
        fmpz_divexact(S + k, u, v);
    }

    if (can & 1) _fmpq_poly_canonicalise(S, Sden, n);
    if (can & 2) _fmpq_poly_canonicalise(C, Cden, n);

    fmpz_clear(t);
    fmpz_clear(u);
    fmpz_clear(v);
}

void
_fmpq_poly_sin_cos_series_basecase(fmpz * S, fmpz_t Sden,
    fmpz * C, fmpz_t Cden, const fmpz * A, const fmpz_t Aden, slong Alen, slong n)
{
    _fmpq_poly_sin_cos_series_basecase_can(S, Sden, C, Cden, A, Aden, Alen, n, 3);
}

void
_fmpq_poly_sin_cos_series_tangent(fmpz * S, fmpz_t Sden,
    fmpz * C, fmpz_t Cden, const fmpz * A, const fmpz_t Aden,
    slong Alen, slong n)
{
    fmpz * t;
    fmpz * u;
    fmpz_t tden;
    fmpz_t uden;

    Alen = FLINT_MIN(Alen, n);

    t = _fmpz_vec_init(n);
    u = _fmpz_vec_init(n);
    fmpz_init(tden);
    fmpz_init(uden);

    /* sin(x) = 2*tan(x/2)/(1+tan(x/2)^2) */
    /* cos(x) = (1-tan(x/2)^2)/(1+tan(x/2)^2) */

    /* t = tan(x/2), u = 1+tan(x/2)^2 */
    fmpz_mul_ui(uden, Aden, 2);
    _fmpq_poly_tan_series(t, tden, A, uden, Alen, n);
    _fmpq_poly_mullow(u, uden, t, tden, n, t, tden, n, n);
    fmpz_set(u, uden);
    _fmpq_poly_canonicalise(u, uden, n);

    /* C = 1/(1+tan(x/2))^2 */
    _fmpq_poly_inv_series(C, Cden, u, uden, n, n);

    /* S = sin(x)/2 */
    _fmpq_poly_mullow(S, Sden, t, tden, n, C, Cden, n, n);
    _fmpq_poly_canonicalise(S, Sden, n);

    /* u = sin(x)/2 * tan(x/2) */
    /* C = C - u */
    _fmpq_poly_mullow(u, uden, S, Sden, n, t, tden, n, n);
    _fmpq_poly_canonicalise(u, uden, n);
    _fmpq_poly_sub(C, Cden, C, Cden, n, u, uden, n);

    /* S = sin(x) */
    _fmpq_poly_scalar_mul_ui(S, Sden, S, Sden, n, 2);

    _fmpz_vec_clear(t, n);
    _fmpz_vec_clear(u, n);
    fmpz_clear(tden);
    fmpz_clear(uden);
}

void
_fmpq_poly_sin_cos_series(fmpz * S, fmpz_t Sden,
    fmpz * C, fmpz_t Cden, const fmpz * A, const fmpz_t Aden,
    slong Alen, slong n)
{
    if (Alen < 20 || n < 20)
        _fmpq_poly_sin_cos_series_basecase(S, Sden, C, Cden, A, Aden, Alen, n);
    else
        _fmpq_poly_sin_cos_series_tangent(S, Sden, C, Cden, A, Aden, Alen, n);
}

void
fmpq_poly_sin_cos_series(fmpq_poly_t res1, fmpq_poly_t res2, const fmpq_poly_t poly, slong n)
{
    if (n == 0)
    {
        fmpq_poly_zero(res1);
        fmpq_poly_zero(res2);
        return;
    }

    if (fmpq_poly_is_zero(poly) || n == 1)
    {
        fmpq_poly_zero(res1);
        fmpq_poly_one(res2);
        return;
    }

    if (!fmpz_is_zero(poly->coeffs))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_poly_sin_cos_series). Constant term != 0.\n");
    }

    fmpq_poly_fit_length(res1, n);
    fmpq_poly_fit_length(res2, n);
    _fmpq_poly_sin_cos_series(res1->coeffs, res1->den,
        res2->coeffs, res2->den, poly->coeffs, poly->den, poly->length, n);
    _fmpq_poly_set_length(res1, n);
    _fmpq_poly_normalise(res1);
    _fmpq_poly_set_length(res2, n);
    _fmpq_poly_normalise(res2);
}

