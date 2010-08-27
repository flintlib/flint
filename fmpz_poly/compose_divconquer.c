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
    long *halloc, *hlen, alloc, powlen;
    fmpz *v, **h, *pow, *temp;
    
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

    /* Run through the algorithm to determine the allocation sizes */
    
    halloc = (long *) malloc(((len1 + 1) / 2) * sizeof(long));
    for (i = 0; i < (len1 + 1) / 2; i++)
        halloc[i] = len2;
    
    powlen = 2 * len2 - 1;
    for (n = (len1 + 1) / 2; n > 2; n = (n + 1) / 2)
    {
        halloc[0] = FLINT_MAX(halloc[0], powlen + halloc[1] - 1);
        
        for (i = 1; i < n / 2; i++)
        {
            halloc[i] = FLINT_MAX(halloc[i], halloc[2 * i + 1] + powlen - 1);
            halloc[i] = FLINT_MAX(halloc[i], halloc[2 * i]);
        }

        if ((n & 1L))
            halloc[i] = FLINT_MAX(halloc[i], halloc[2 * i]);
        
        powlen += powlen - 1;
    }
    
    /* Complete the initialisation */
    
    alloc = 0;
    for (i = 0; i < (len1 + 1) / 2; i++)
        alloc += halloc[i];

    v = _fmpz_vec_init(alloc +  2 * powlen);
    h = (fmpz **) malloc(((len1 + 1) / 2) * sizeof(fmpz *));
    hlen = (long *) calloc(((len1 + 1) / 2), sizeof(long));
    h[0] = v;
    for (i = 0; i < (len1 - 1) / 2; i++)
        h[i + 1] = h[i] + halloc[i];
    pow  = v + alloc;
    temp = pow + powlen;
    
    /* Let's start the actual work */
    
    j = 0;
    for (i = 0; i < len1 / 2; i++)
    {
        if (poly1[j + 1] != 0L)
        {
            _fmpz_vec_scalar_mul_fmpz(h[i], poly2, len2, poly1 + j + 1);
            fmpz_add(h[i], h[i], poly1 + j);
            hlen[i] = len2;
        }
        else if (poly1[j] != 0L)
        {
            fmpz_set(h[i], poly1 + j);
            hlen[i] = 1;
        }
        j += 2;
    }
    if ((len1 & 1L))
    {
        if (poly1[j] != 0L)
        {
            fmpz_set(h[i], poly1 + j);
            hlen[i] = 1;
        }
    }
    
    _fmpz_poly_mul(pow, poly2, len2, poly2, len2);
    powlen = 2 * len2 - 1;
    
    for (n = (len1 + 1) / 2; n > 2; n = (n + 1) / 2)
    {
        if (hlen[1] > 0)
        {
            long templen = powlen + hlen[1] - 1;
            _fmpz_poly_mul(temp, pow, powlen, h[1], hlen[1]);
            _fmpz_poly_add(h[0], temp, templen, h[0], hlen[0]);
            hlen[0] = FLINT_MAX(hlen[0], templen);
        }
        
        for (i = 1; i < n / 2; i++)
        {
            if (hlen[2*i + 1] > 0)
            {
                _fmpz_poly_mul(h[i], pow, powlen, h[2*i + 1], hlen[2*i + 1]);
                hlen[i] = hlen[2*i + 1] + powlen - 1;
            }
            _fmpz_poly_add(h[i], h[i], hlen[i], h[2*i], hlen[2*i]);
            hlen[i] = FLINT_MAX(hlen[i], hlen[2*i]);
        }
        if ((n & 1L))
        {
            _fmpz_vec_copy(h[i], h[2*i], hlen[2*i]);
            hlen[i] = hlen[2*i];
        }
        
        _fmpz_poly_mul(temp, pow, powlen, pow, powlen);
        powlen += powlen - 1;
        {
            fmpz * t = pow;
            pow      = temp;
            temp     = t;
        }
    }

    _fmpz_poly_mul(res, pow, powlen, h[1], hlen[1]);
    _fmpz_vec_add(res, res, h[0], hlen[0]);
    
    free(halloc);
    _fmpz_vec_clear(v, alloc + 2 * powlen);
    free(h);
    free(hlen);
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
