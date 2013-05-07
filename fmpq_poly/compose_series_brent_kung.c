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
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "fmpq.h"
#include "fmpq_mat.h"
#include "ulong_extras.h"

static void
_fmpq_mat_get_row(fmpz * rnum, fmpz_t den, fmpq_mat_t A, long i)
{
    long j;
    fmpz_t t;
    fmpz_init(t);
    fmpz_one(den);

    for (j = 0; j < fmpq_mat_ncols(A); j++)
        fmpz_lcm(den, den, fmpq_mat_entry_den(A, i, j));

    for (j = 0; j < fmpq_mat_ncols(A); j++)
    {
        fmpz_divexact(t, den, fmpq_mat_entry_den(A, i, j));
        fmpz_mul(rnum + j, fmpq_mat_entry_num(A, i, j), t);
    }

    fmpz_clear(t);
}


void
_fmpq_poly_compose_series_brent_kung(fmpz * res, fmpz_t den, const fmpz * poly1,
        const fmpz_t den1, long len1, const fmpz * poly2,
        const fmpz_t den2, long len2, long n)
{
    fmpq_mat_t A, B, C;
    fmpz_t tden, uden, hden;
    fmpz *t, *u, *h, *swap;
    long i, j, m;

    if (fmpz_is_one(den2))
    {
        _fmpz_poly_compose_series(res, poly1, len1, poly2, len2, n);
        fmpz_set(den, den1);
        _fmpq_poly_canonicalise(res, den, n);
        return;
    }

    if (n == 1)
    {
        fmpz_set(res, poly1);
        fmpz_set(den, den1);
        _fmpq_poly_canonicalise(res, den, 1);
        return;
    }

    m = n_sqrt(n) + 1;

    fmpq_mat_init(A, m, n);
    fmpq_mat_init(B, m, m);
    fmpq_mat_init(C, m, n);

    fmpz_init(tden);
    fmpz_init(uden);
    fmpz_init(hden);
    h = _fmpz_vec_init(n);
    t = _fmpz_vec_init(n);
    u = _fmpz_vec_init(n);

    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < len1; i++)
    {
        fmpz_set(fmpq_mat_entry_num(B, i / m, i % m), poly1 + i);
        fmpz_set(fmpq_mat_entry_den(B, i / m, i % m), den1);
        fmpq_canonicalise(fmpq_mat_entry(B, i / m, i % m));
    }

    /* Set rows of A to powers of poly2 */
    fmpq_set_si(fmpq_mat_entry(A, 0, 0), 1L, 1L);

    for (i = 0; i < len2; i++)
    {
        fmpz_set(fmpq_mat_entry_num(A, 1, i), poly2 + i);
        fmpz_set(fmpq_mat_entry_den(A, 1, i), den2);
        fmpq_canonicalise(fmpq_mat_entry(A, 1, i));
    }

    _fmpz_vec_set(h, poly2, len2);
    fmpz_set(hden, den2);

    for (i = 2; i < m; i++)
    {
        _fmpq_poly_mullow(t, tden, h, hden, n, poly2, den2, len2, n);
        _fmpq_poly_canonicalise(t, tden, n);

        for (j = 0; j < n; j++)
        {
            fmpz_set(fmpq_mat_entry_num(A, i, j), t + j);
            fmpz_set(fmpq_mat_entry_den(A, i, j), tden);
            fmpq_canonicalise(fmpq_mat_entry(A, i, j));
        }
        swap = t; t = h; h = swap;
        fmpz_swap(hden, tden);
    }

    /* Compute h = poly2 ^ m */
    _fmpq_poly_mullow(t, tden, h, hden, n, poly2, den2, len2, n);
    _fmpq_poly_canonicalise(t, tden, n);
    swap = t; t = h; h = swap;
    fmpz_swap(hden, tden);

    /* Matrix multiply */
    fmpq_mat_mul(C, B, A);
    fmpq_mat_clear(A);
    fmpq_mat_clear(B);

    /* Evaluate block composition using the Horner scheme */
    _fmpq_mat_get_row(res, den, C, m - 1);

    for (i = m - 2; i >= 0; i--)
    {
        _fmpq_poly_mullow(t, tden, res, den, n, h, hden, n, n);
        /* we could canonicalise t here, but it does not seem to make
           much of a difference */
        _fmpq_mat_get_row(u, uden, C, i);
        _fmpq_poly_add(res, den, t, tden, n, u, uden, n);
    }

    _fmpq_poly_canonicalise(res, den, n);

    fmpq_mat_clear(C);

    _fmpz_vec_clear(t, n);
    _fmpz_vec_clear(u, n);
    _fmpz_vec_clear(h, n);
    fmpz_clear(tden);
    fmpz_clear(uden);
    fmpz_clear(hden);
}

void
fmpq_poly_compose_series_brent_kung(fmpq_poly_t res, 
                    const fmpq_poly_t poly1, const fmpq_poly_t poly2, long n)
{
    long len1 = poly1->length;
    long len2 = poly2->length;
    long lenr;

    if (len2 != 0 && !fmpz_is_zero(poly2->coeffs))
    {
        printf("Exception (fmpq_poly_compose_series_brent_kung). \n"
               "Inner polynomial must have zero constant term.\n");
        abort();
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
