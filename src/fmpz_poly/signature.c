/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void _fmpz_poly_signature(slong * r1, slong * r2, const fmpz * poly, slong len)
{
    fmpz *A, *B, *f, *g, *h, *w;
    slong lenA, lenB;
    int s, t;

    if (len <= 2)
    {
        *r1 = (len == 2);
        *r2 = 0;
        return;
    }

    w = _fmpz_vec_init(2 * len + 2);
    A = w;
    B = w + len;
    lenA = len;
    lenB = lenA - 1;
    f = w + 2 * len - 1;
    g = w + 2 * len;
    h = w + 2 * len + 1;

    _fmpz_poly_primitive_part(A, poly, lenA);
    _fmpz_poly_derivative(B, A, lenA);
    _fmpz_poly_primitive_part(B, B, lenB);

    fmpz_one(g);
    fmpz_one(h);

    s = 1;
    t = (lenA & WORD(1)) ? -s : s;
    *r1 = 1;

    while (1)
	{
        slong delta = lenA - lenB;
        int sgnA;

        _fmpz_poly_pseudo_rem_cohen(A, A, lenA, B, lenB);

        lenA = lenB;
        FMPZ_VEC_NORM(A, lenA);

		if (lenA == 0)
		{
            _fmpz_vec_clear(w, 2 * len + 2);
			flint_throw(FLINT_ERROR, "(fmpz_poly_signature): Non-squarefree polynomial detected.\n");
		}

        if ((fmpz_sgn(B + (lenB - 1)) > 0) || (delta & WORD(1)))
            _fmpz_vec_neg(A, A, lenA);

        sgnA = fmpz_sgn(A + (lenA - 1));
		if (sgnA != s)
		{
			s = -s;
			(*r1)--;
        }
		if (sgnA != ((lenA & WORD(1)) ? t : -t))
		{
			t = -t;
			(*r1)++;
		}

        if (lenA == 1)
        {
            *r2 = ((len - 1) - *r1) / 2;

            _fmpz_vec_clear(w, 2 * len + 2);
            return;
        }
		else
		{
            {
                fmpz * temp = A;
                A = B;
                B = temp;
            }
            {
                slong temp = lenA;
                lenA = lenB;
                lenB = temp;
            }

            if (delta == 1)
            {
                fmpz_mul(f, g, h);
                _fmpz_vec_scalar_divexact_fmpz(B, B, lenB, f);
                fmpz_abs(g, A + (lenA - 1));
                fmpz_set(h, g);
            }
            else
            {
                fmpz_pow_ui(f, h, delta);
                fmpz_mul(f, f, g);
                _fmpz_vec_scalar_divexact_fmpz(B, B, lenB, f);
                fmpz_pow_ui(f, h, delta - 1);
                fmpz_pow_ui(g, A + (lenA - 1), delta);
                fmpz_abs(g, g);
                fmpz_divexact(h, g, f);
                fmpz_abs(g, A + (lenA - 1));
            }
		}
	}
}

void fmpz_poly_signature(slong * r1, slong * r2, const fmpz_poly_t poly)
{
    _fmpz_poly_signature(r1, r2, poly->coeffs, poly->length);
}

