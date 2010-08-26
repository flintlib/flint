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
    long i, j, n;
    long *alloc, talloc;
    long *half_len, *temp_len, pow2_len;
    fmpz *v, **half, **temp, *pow2, *pow2t;
    
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
    
    /*
       Run through the algorithm to figure out the allocations:
       
       -  alloc[i] is the maximum length of half[i] and temp[i]
       -  talloc is the sum of these values
     */
    
    alloc = (long *) malloc(((len1 + 1) / 2) * sizeof(long));
    
    for (i = 0; i < (len1 + 1) / 2; i++)
        alloc[i] = len2;
    pow2_len = 2 * len2 + 1;
    for (n = (len1 + 1) / 2; n > 2; n = (n + 1) / 2)
    {
        pow2_len += pow2_len - 1;
        for (i = 0; i < n / 2; i++)
            alloc[i] = alloc[2 * i + 1] + pow2_len - 1;
        if ((n & 1L))
            alloc[i] = alloc[2 * i];
        pow2_len += pow2_len - 1;
    }

    talloc = 0;
    for (i = 0; i < (len1 + 1) / 2; i++)
        talloc += alloc[i];

    /*
       Allocate
     */

    v    = _fmpz_vec_init(2 * talloc + 2 * pow2_len);
    
    half = (fmpz **) malloc(((len1 & 1L) ? len1 + 1 : len1) * sizeof(fmpz *));
    temp = half + ((len1 + 1) / 2);
    
    half[0] = v;
    temp[0] = v + talloc;
    
    half_len = (long *) calloc(2 * talloc, sizeof(long));
    temp_len = half_len + talloc;
    
    for (i = 1; i < (len1 + 1) / 2; i++)
    {
        half[i] = half[i-1] + alloc[i-1];
        temp[i] = temp[i-1] + alloc[i-1];
    }
    
    pow2  = v + 2 * talloc;
    pow2t = pow2 + pow2_len;

    /*
       Let's start the actual work
     */

    j = 0;
    for (i = 0; i < len1 / 2; i++)
    {
        if (poly1[j + 1] != 0L)
        {
            _fmpz_vec_scalar_mul_fmpz(temp[i], poly2, len2, poly1 + j + 1);
            fmpz_add(temp[i], temp[i], poly1 + j);
            temp_len[i] = len2;
        }
        else if (poly1[j] != 0L)
        {
            fmpz_set(temp[i], poly1 + j);
            temp_len[i] = 1;
        }
        j += 2;
    }
    if ((len1 & 1L))
    {
        if (poly1[j] != 0L)
        {
            fmpz_set(temp[i], poly1 + j);
            temp_len[i] = 1;
        }
    }
    
    _fmpz_poly_mul(pow2, poly2, len2, poly2, len2);
    pow2_len = 2 * len2 - 1;
    
    for (n = (len1 + 1) / 2; n > 2; n = (n + 1) / 2)
    {
        for (i = 0; i < n / 2; i++)
        {
            if (temp_len[2 * i + 1] > 0)
            {
                _fmpz_poly_mul(half[i], pow2, pow2_len, temp[2 * i + 1], temp_len[2 * i + 1]);
                half_len[i] = temp_len[2 * i + 1] + pow2_len - 1;
                
                _fmpz_poly_add(half[i], half[i], half_len[i], temp[2 * i], temp_len[2 * i]);
            }
            else
            {
                _fmpz_vec_copy(half[i], temp[2 * i], temp_len[2 * i]);
                half_len[i] = temp_len[2 * i];
            }
        }
        if ((n & 1L))
        {
            _fmpz_vec_copy(half[i], temp[2 * i], temp_len[2 * i]);
            half_len[i] = temp_len[2 * i];
        }

        _fmpz_poly_mul(pow2t, pow2, pow2_len, pow2, pow2_len);
        pow2_len += pow2_len - 1;
        {
            fmpz * t = pow2t;
            pow2t    = pow2;
            pow2     = t;
        }
        
        for (i = 0; i < n; i++)
        {
            {
                fmpz * t = half[i];
                half[i]  = temp[i];
                temp[i]  = t;
            }
            {
                long t = half_len[i];
                half_len[i] = temp_len[i];
                temp_len[i] = t;
            }
        }
    }
    
    _fmpz_poly_mul(res, pow2, pow2_len, temp[1], temp_len[1]);
    _fmpz_vec_add(res, res, temp[0], temp_len[0]);

    free(alloc);
    _fmpz_vec_clear(v, 2 * talloc + 2 * pow2_len);
    free(half);
    free(half_len);
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
