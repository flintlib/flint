/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void _fmpz_poly_num_real_roots_sturm(slong * n_neg, slong * n_pos, const fmpz * pol, slong len)
{
    fmpz_t a, b, g, h;
    fmpz *A, *B, *W;
    slong lenA, lenB;
    int s, s0a, s0, s0b;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(g);
    fmpz_init(h);

    A = W = _fmpz_vec_init(2*len - 1);
    B = W + len;
    lenA = len;
    lenB = len - 1;

    _fmpz_poly_content(a, pol, lenA);
    _fmpz_vec_scalar_divexact_fmpz(A, pol, lenA, a);
    _fmpz_poly_derivative(B, A, lenA);
    _fmpz_poly_content(b, B, lenB);
    _fmpz_vec_scalar_divexact_fmpz(B, B, lenB, b);

    fmpz_one(g);
    fmpz_one(h);

    s0a = s0b = fmpz_sgn(A + (lenA - 1));
    if (lenA % 2 == 0)
        s0b = -s0b;
    s0 = fmpz_sgn(A);
    (*n_neg) = (*n_pos) = 0;
    while (1)
    {
        const slong delta = lenA - lenB;

        /* sign change at + infinity */
        s = fmpz_sgn(B + (lenB - 1));
        if (s != s0a)
        {
            (*n_pos)--;
            s0a = s;
        }
        /* sign change at - infinity */
        if (lenB % 2 == 0) s = -s;
        if (s != s0b)
        {
            (*n_neg)++;
            s0b = s;
        }
        /* sign change at 0 */
        s = fmpz_sgn(B);
        if (s && (s != s0))
        {
            (*n_neg)--;
            (*n_pos)++;
            s0 = s;
        }

        /* now compute the next element of the remainder sequence */
        _fmpz_poly_pseudo_rem_cohen(A, A, lenA, B, lenB);
        if ((fmpz_sgn(B + (lenB - 1)) > 0) || ((lenA ^ lenB) & 1))
        {
            _fmpz_vec_neg(A, A, lenA);
        }
        FMPZ_VEC_NORM(A, lenA);

        if (lenA <= 1)
            break;

        {
            fmpz *T;
            slong len;
            T = A, A = B, B = T, len = lenA, lenA = lenB, lenB = len;
        }

        if (delta == 1)
        {
            fmpz_mul(b, g, h);
            if (fmpz_sgn(b) < 0) fmpz_neg(b, b);
            _fmpz_vec_scalar_divexact_fmpz(B, B, lenB, b);
            fmpz_set(g, A + (lenA - 1));
            fmpz_set(h, g);
        }
        else
        {
            fmpz_pow_ui(a, h, delta);
            fmpz_mul(b, g, a);
            if (fmpz_sgn(b) < 0) fmpz_neg(b, b);
            _fmpz_vec_scalar_divexact_fmpz(B, B, lenB, b);
            fmpz_pow_ui(b, A + (lenA - 1), delta);
            fmpz_mul(g, h, b);
            fmpz_divexact(h, g, a);
            fmpz_set(g, A + (lenA - 1));
        }
    }
    if (lenA)
    {
        /* sign change at + infinity */
        s = fmpz_sgn(A);
        if (s != s0a)
            (*n_pos)--;
        /* sign change at - infinity */
        if (s != s0b)
            (*n_neg)++;
        /* sign change at 0 */
        if (s != s0)
        {
            (*n_neg)--;
            (*n_pos)++;
        }
    }

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(g);
    fmpz_clear(h);

    _fmpz_vec_clear(W, 2*len - 1);
}


slong fmpz_poly_num_real_roots_sturm(const fmpz_poly_t pol)
{
    slong i, len;
    slong n_neg = 0;
    slong n_pos = 0;

    if (fmpz_poly_is_zero(pol))
        flint_throw(FLINT_ERROR, "Zero polynomial in %s\n", __func__);

    for (i = 0; (i < pol->length) && fmpz_is_zero(pol->coeffs + i); i++);
    len = pol->length - i;

    if (len == 1)
        return i;
    else if (len == 2)
        return i + 1;
    else
    {
        _fmpz_poly_num_real_roots_sturm(&n_neg, &n_pos, i + pol->coeffs, len);
        return i + n_neg + n_pos;
    }
}
