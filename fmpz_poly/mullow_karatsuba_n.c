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

    Copyright (C) 2010 William Hart

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void 
_fmpz_poly_mullow_kara_recursive(fmpz * out, const fmpz * pol1,
                                 const fmpz * pol2, fmpz * temp, len_t len);

/*
   Multiplication using truncated karatsuba.
   
   Below length 7, classical truncated multiplication is always theoretically 
   faster, so we switch to that as the basecase.  Above that we use the 
   ordinary (left/right) karatsuba identity and recursively do one full
   karatsuba multiplication and two truncated karatsuba multiplications.
 */

void
_fmpz_poly_mullow_kara_recursive(fmpz * out, const fmpz * pol1,
                                 const fmpz * pol2, fmpz * temp, len_t len)
{
    len_t m1 = len / 2;
    len_t m2 = len - m1;
    int odd = (len & 1);

    if (len <= 6)
    {
        _fmpz_poly_mullow_classical(out, pol1, len, pol2, len, len);
        return;
    }

    _fmpz_vec_add(temp + m2, pol1, pol1 + m1, m1);
    if (odd)
        fmpz_set(temp + m2 + m1, pol1 + 2 * m1);

    _fmpz_vec_add(temp + 2 * m2, pol2, pol2 + m1, m1);
    if (odd)
        fmpz_set(temp + 2 * m2 + m1, pol2 + 2 * m1);

    _fmpz_poly_mul_karatsuba(out, pol1, m1, pol2, m1);
    fmpz_zero(out + 2 * m1 - 1);

    _fmpz_poly_mullow_kara_recursive(temp, temp + m2, temp + 2 * m2,
                                     temp + 3 * m2, m2);

    _fmpz_poly_mullow_kara_recursive(temp + m2, pol1 + m1, pol2 + m1,
                                     temp + 2 * m2, m2);

    _fmpz_vec_sub(temp, temp, out, m2);
    _fmpz_vec_sub(temp, temp, temp + m2, m2);

    if (odd)
        fmpz_set(out + 2 * m1, temp + m2);
    _fmpz_vec_add(out + m1, out + m1, temp, m2);
}

/* Assumes poly1 and poly2 are not length 0. */
void
_fmpz_poly_mullow_karatsuba_n(fmpz * res, const fmpz * poly1,
                              const fmpz * poly2, len_t n)
{
    fmpz *temp;
    len_t len, loglen = 0;

    if (n == 1)
    {
        fmpz_mul(res, poly1, poly2);
        return;
    }

    while ((1L << loglen) < n)
        loglen++;
    len = (1L << loglen);

    temp = _fmpz_vec_init(3 * len);

    _fmpz_poly_mullow_kara_recursive(res, poly1, poly2, temp, n);

    _fmpz_vec_clear(temp, 3 * len);
}

void
fmpz_poly_mullow_karatsuba_n(fmpz_poly_t res, const fmpz_poly_t poly1, 
                             const fmpz_poly_t poly2, len_t n)
{
    const len_t len1 = FLINT_MIN(poly1->length, n);
    const len_t len2 = FLINT_MIN(poly2->length, n);
    len_t i, lenr;

    int clear = 0;
    fmpz *copy1, *copy2;

    if (len1 == 0 || len2 == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    lenr = len1 + len2 - 1;
    if (n > lenr)
        n = lenr;

    if (len1 >= n)
        copy1 = poly1->coeffs;
    else
    {
        copy1 = (fmpz *) flint_malloc(n * sizeof(fmpz));
        for (i = 0; i < len1; i++)
            copy1[i] = poly1->coeffs[i];
        flint_mpn_zero((mp_ptr) copy1 + len1, n - len1);
        clear |= 1;
    }

    if (len2 >= n)
        copy2 = poly2->coeffs;
    else
    {
        copy2 = (fmpz *) flint_malloc(n * sizeof(fmpz));
        for (i = 0; i < len2; i++)
            copy2[i] = poly2->coeffs[i];
        flint_mpn_zero((mp_ptr) copy2 + len2, n - len2);
        clear |= 2;
    }

    if (res != poly1 && res != poly2)
    {
        fmpz_poly_fit_length(res, n);
        _fmpz_poly_mullow_karatsuba_n(res->coeffs, copy1, copy2, n);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        _fmpz_poly_mullow_karatsuba_n(t->coeffs, copy1, copy2, n);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
    _fmpz_poly_set_length(res, n);
    _fmpz_poly_normalise(res);

    if (clear & 1)
        flint_free(copy1);
    if (clear & 2)
        flint_free(copy2);
}
