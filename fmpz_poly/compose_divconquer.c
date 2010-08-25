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

******************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_compose_divconquer(fmpz * res, const fmpz * poly1, long len1, 
                                          const fmpz * poly2, long len2)
{
    long alloc, i, j, n;
    fmpz_poly_t *w, *half, *temp, pow2, pow2t;
    
	if (len1 == 1)
	{
        _fmpz_vec_copy(res, poly1, len1);
		return;
	}
    if (len2 == 1)
    {
        _fmpz_poly_evaluate(res, poly1, len1, poly2);
        return;
    }
	if (len1 == 2)
	{
		_fmpz_poly_compose_horner(res, poly1, len1, poly2, len2);
		return;
	}

    alloc = (len1 & 1L) ? len1 + 1 : len1;
    w     = (fmpz_poly_t *) malloc(alloc * sizeof(fmpz_poly_t));
    half  = w;
    temp  = w + ((len1 + 1) / 2);

    for (i = 0; i < (len1 + 1) / 2; i++)
        fmpz_poly_init2(temp[i], len2);
    for (i = 0; i < (len1 + 1) / 2; i++)
        fmpz_poly_init(half[i]);
    
    j = 0;
    for (i = 0; i < len1 / 2; i++)
    {
        if (poly1[j + 1] != 0L)
        {
            _fmpz_vec_scalar_mul_fmpz(temp[i]->coeffs, poly2, len2, poly1 + j + 1);
            fmpz_add(temp[i]->coeffs, temp[i]->coeffs, poly1 + j);
            temp[i]->length = len2;
        }
        else if (poly1[j] != 0L)
        {
            fmpz_set(temp[i]->coeffs, poly1 + j);
            temp[i]->length = 1;
        }
        j += 2;
    }
    if ((len1 & 1L))
    {
        if (poly1[j] != 0L)
        {
            fmpz_set(temp[i]->coeffs, poly1 + j);
            temp[i]->length = 1;
        }
    }
    
    fmpz_poly_init2(pow2, (len1 - 1) * (len2 - 1) + 1);
    fmpz_poly_init2(pow2t, (len1 - 1) * (len2 - 1) + 1);
    
    _fmpz_poly_mul(pow2->coeffs, poly2, len2, poly2, len2);
    _fmpz_poly_set_length(pow2, 2 * len2 - 1);
    
    for (n = (len1 + 1) / 2; n > 2; n = (n + 1) / 2)
    {
        for (i = 0; i < n / 2; i++)
        {
            fmpz_poly_mul(half[i], temp[2 * i + 1], pow2);
            fmpz_poly_add(half[i], half[i], temp[2 * i]);
        }

        if ((n & 1L))
        {
            fmpz_poly_set(half[i], temp[2 * i]);
        }
        
        fmpz_poly_mul(pow2t, pow2, pow2);
        fmpz_poly_swap(pow2, pow2t);
        
        for (i = 0; i < n; i++)
        {
            fmpz_poly_swap(half[i], temp[i]);
        }
    }

    _fmpz_poly_mul(res, pow2->coeffs, pow2->length, 
                        temp[1]->coeffs, temp[1]->length);
    _fmpz_vec_add(res, res, temp[0]->coeffs, temp[0]->length);

    fmpz_poly_clear(pow2);
    fmpz_poly_clear(pow2t);
    for (i = 0; i < alloc; i++)
        fmpz_poly_clear(w[i]);
    free(w);
}

void
fmpz_poly_compose_divconquer(fmpz_poly_t res, 
                             const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    const long len1 = poly1->length;
    const long len2 = poly2->length;
    long lenr;
    
    if (len1 == 0)
    {
        fmpz_poly_zero(res);
        return;
    }
    if (len1 == 1 || len2 == 0)
    {
        fmpz_poly_fit_length(res, 1);
        fmpz_set(res->coeffs, poly1->coeffs);
        _fmpz_poly_set_length(res, 1);
        _fmpz_poly_normalise(res);
        return;
    }
    
    lenr = (len1 - 1) * (len2 - 1) + 1;
    
    if ((res != poly1) && (res != poly2))
    {
        fmpz_poly_fit_length(res, lenr);
        _fmpz_poly_compose_divconquer(res->coeffs, poly1->coeffs, len1, 
                                                   poly2->coeffs, len2);
        _fmpz_poly_set_length(res, lenr);
        _fmpz_poly_normalise(res);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, lenr);
        _fmpz_poly_compose_divconquer(t->coeffs, poly1->coeffs, len1,
                                                 poly2->coeffs, len2);
        _fmpz_poly_set_length(t, lenr);
        _fmpz_poly_normalise(t);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
}
