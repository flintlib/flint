/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"


slong _fmpz_poly_num_real_roots_sturm(fmpz * pol, slong len)
{
    fmpz_t a, b, g, h;
    fmpz *A, *B, *W;
    slong lenA, lenB;
    int s, s0a, s0b, t;

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

    t = 0;
    s0a = s0b = fmpz_sgn(A + (lenA - 1));
    if (lenA % 2 == 0)
        s0b = -s0b;
    while (1)
    {
        const slong delta = lenA - lenB;

        /* sign change at + infinity */
        s = fmpz_sgn(B + (lenB - 1));
        if (s != s0a)
        {
            t--;
            s0a = s;
        }
        /* sign change at - infinity */
        if (lenB % 2 == 0) s = -s;
        if (s != s0b)
        {
            t++;
            s0b = s;
        }

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
            t--;
        /* sign change at - infinity */
        if (s != s0b)
            t++;
    }

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(g);
    fmpz_clear(h);

    _fmpz_vec_clear(W, 2*len - 1);

    return t;
}


slong fmpz_poly_num_real_roots_sturm(fmpz_poly_t pol)
{
    if (fmpz_poly_is_zero(pol))
    {
        printf("ERROR (fmpz_poly_num_real_roots_sturm): zero polynomial\n");
        flint_abort();
    }
    if (pol->length == 1)
        return 0;
    if (pol->length == 2)
        return 1;

    return _fmpz_poly_num_real_roots_sturm(pol->coeffs, pol->length);
}
