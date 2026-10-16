/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_mat.h"
#include "fmpq_poly.h"

void
_fmpq_poly_compose_series_brent_kung(fmpz * res, fmpz_t den, const fmpz * poly1,
        const fmpz_t den1, slong len1, const fmpz * poly2,
        const fmpz_t den2, slong len2, slong n)
{
    fmpz_mat_t A, B, C;
    fmpz *A_den, *t, *h;
    fmpz_t lcd, scale;
    fmpz_t C_den, tden, hden;
    slong i, m;

    if (n == 1)
    {
        fmpz_set(res, poly1);
        fmpz_set(den, den1);
        _fmpq_poly_canonicalise(res, den, 1);
        return;
    }

    m = n_sqrt(n) + 1;

    fmpz_mat_init(A, m, n);
    fmpz_mat_init(B, m, m);
    fmpz_mat_init(C, m, n);
    A_den = _fmpz_vec_init(m);

    fmpz_init(lcd);
    fmpz_init(scale);
    fmpz_init(C_den);
    fmpz_init(tden);
    fmpz_init(hden);
    h = _fmpz_vec_init(n);
    t = _fmpz_vec_init(n);

    /* Set rows of B to the segments of poly1 with common denominator den1. */
    for (i = 0; i < len1; i++)
        fmpz_set(fmpz_mat_entry(B, i / m, i % m), poly1 + i);
    /* Remark: if the poly is non-canonical e.g. due to being a truncation
       of a longer power series, it could be helpful to remove its content
       here. We could also consider removing content of B row by row. */

    /* Set rows of A to the numerators of powers of poly2 with corresponding
       denominators in A_den[i]. */
    fmpz_one(fmpz_mat_entry(A, 0, 0));
    fmpz_one(A_den + 0);
    _fmpz_vec_set(fmpz_mat_row(A, 1), poly2, len2);
    fmpz_set(A_den + 1, den2);
    /* Optional: may improve performance if poly2 is non-canonical e.g. due
       to being a truncation of a longer power series. */
    _fmpq_poly_canonicalise(fmpz_mat_row(A, 1), A_den + 1, n);

    for (i = 2; i < m; i++)
    {
        fmpz * Ai = fmpz_mat_row(A, i);
        fmpz * Ai1 = fmpz_mat_row(A, i - 1);
        fmpz * Ai_den = A_den + i;
        fmpz * Ai1_den = A_den + i - 1;

        _fmpq_poly_mullow(Ai, Ai_den, Ai1, Ai1_den, n, poly2, den2, len2, n);
        _fmpq_poly_canonicalise(Ai, Ai_den, n);
    }

    /* Compute h = poly2 ^ m */
    _fmpq_poly_mullow(h, hden, fmpz_mat_row(A, m - 1), A_den + m - 1, n, poly2, den2, len2, n);
    _fmpq_poly_canonicalise(h, hden, n);

    /* Matrix multiply C = B * A over the integers.
       B (m x m) has common denominator den1.
       A (m x n) has row denominator A_den[i].

       We could remove gcd of A columnwise before multiplying, but
       this is slower for small to moderate n where we want to use Brent-Kung
       and faster only for large n where we want to use Kinoshita-Li instead. */
    fmpz_one(lcd);
    for (i = 0; i < m; i++)
        fmpz_lcm(lcd, lcd, A_den + i);

    for (i = 0; i < m; i++)
    {
        fmpz_divexact(scale, lcd, A_den + i);
        if (!fmpz_is_one(scale))
            _fmpz_vec_scalar_mul_fmpz(fmpz_mat_row(A, i), fmpz_mat_row(A, i), n, scale);
    }

    fmpz_mat_mul(C, B, A);
    fmpz_mat_clear(A);
    _fmpz_vec_clear(A_den, m);
    fmpz_mat_clear(B);

    fmpz_mul(C_den, den1, lcd);

    _fmpz_vec_set(res, fmpz_mat_row(C, m - 1), n);
    fmpz_set(den, C_den);
    _fmpq_poly_canonicalise(res, den, n);
    /* Evaluate block composition using the Horner scheme */
    for (i = m - 2; i >= 0; i--)
    {
        _fmpq_poly_mullow(t, tden, res, den, n, h, hden, n, n);
        /* Remark: we could canonicalise t and/or fmpz_mat_row(C, i)
           here; in practice this seems to be slower at least for moderate n. */
        _fmpq_poly_add(res, den, t, tden, n, fmpz_mat_row(C, i), C_den, n);
    }

    _fmpq_poly_canonicalise(res, den, n);

    fmpz_mat_clear(C);

    _fmpz_vec_clear(t, n);
    _fmpz_vec_clear(h, n);
    fmpz_clear(lcd);
    fmpz_clear(scale);
    fmpz_clear(C_den);
    fmpz_clear(tden);
    fmpz_clear(hden);
}

void
fmpq_poly_compose_series_brent_kung(fmpq_poly_t res,
                    const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong lenr;

    if (len2 != 0 && !fmpz_is_zero(poly2->coeffs))
    {
        flint_throw(FLINT_ERROR, "(fmpq_poly_compose_series_brent_kung): Inner polynomial must have zero constant term.\n");
    }

    if (len1 == 0 || n == 0)
    {
        fmpq_poly_zero(res);
        return;
    }

    if (len2 == 0 || len1 == 1)
    {
        fmpq_poly_fit_length(res, 1);
        fmpz_set(res->coeffs, poly1->coeffs);
        fmpz_set(res->den, poly1->den);
        {
            fmpz_t d;
            fmpz_init(d);
            fmpz_gcd(d, res->coeffs, res->den);
            if (!fmpz_is_one(d))
            {
                fmpz_divexact(res->coeffs, res->coeffs, d);
                fmpz_divexact(res->den, res->den, d);
            }
            fmpz_clear(d);
        }
        _fmpq_poly_set_length(res, 1);
        _fmpq_poly_normalise(res);
        return;
    }

    lenr = FLINT_MIN((len1 - 1) * (len2 - 1) + 1, n);
    len1 = FLINT_MIN(len1, lenr);
    len2 = FLINT_MIN(len2, lenr);

    if ((res != poly1) && (res != poly2))
    {
        fmpq_poly_fit_length(res, lenr);
        _fmpq_poly_compose_series_brent_kung(res->coeffs, res->den,
                           poly1->coeffs, poly1->den, len1,
                           poly2->coeffs, poly2->den, len2, lenr);
        _fmpq_poly_set_length(res, lenr);
        _fmpq_poly_normalise(res);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, lenr);
        _fmpq_poly_compose_series_brent_kung(t->coeffs, t->den,
                           poly1->coeffs, poly1->den, len1,
                           poly2->coeffs, poly2->den, len2, lenr);
        _fmpq_poly_set_length(t, lenr);
        _fmpq_poly_normalise(t);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }
}
