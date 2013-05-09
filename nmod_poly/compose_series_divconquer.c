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

    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2010 William Hart

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

/*
    See fmpz_poly/compose_divconquer.c
 */

void
_nmod_poly_compose_series_divconquer(mp_ptr res, mp_srcptr poly1, len_t len1, 
                                                 mp_srcptr poly2, len_t len2, 
                                                 len_t N, nmod_t mod)
{
    len_t i, j, k, n;
    len_t *hlen, alloc, powlen;
    mp_ptr v, *h, pow, temp;
    
    if (len1 == 1)
    {
        res[0] = poly1[0];
        return;
    }
    if (len2 == 1)
    {
        res[0] = _nmod_poly_evaluate_nmod(poly1, len1, poly2[0], mod);
        return;
    }
    if (len1 == 2)
    {
        mp_limb_t t = poly1[0];
        _nmod_vec_scalar_mul_nmod(res, poly2, len2, poly1[1], mod);
        res[0] = n_addmod(res[0], t, mod.n);
        return;

    }

    /* Initialisation */
    
    hlen = (len_t *) flint_malloc(((len1 + 1) / 2) * sizeof(len_t));
    
    for (k = 1; (2 << k) < len1; k++) ;
    
    hlen[0] = hlen[1] = FLINT_MIN(N, ((1 << k) - 1) * (len2 - 1) + 1);
    for (i = k - 1; i > 0; i--)
    {
        len_t hi = (len1 + (1 << i) - 1) / (1 << i);
        len_t t  = FLINT_MIN(N, ((1 << i) - 1) * (len2 - 1) + 1);
        for (n = (hi + 1) / 2; n < hi; n++)
            hlen[n] = t;
    }
    powlen = FLINT_MIN(N, (1 << k) * (len2 - 1) + 1);
    
    alloc = 0;
    for (i = 0; i < (len1 + 1) / 2; i++)
        alloc += hlen[i];

    v = _nmod_vec_init(alloc +  2 * powlen);
    h = (mp_ptr *) flint_malloc(((len1 + 1) / 2) * sizeof(mp_ptr));
    h[0] = v;
    for (i = 0; i < (len1 - 1) / 2; i++)
    {
        h[i + 1] = h[i] + hlen[i];
        hlen[i]  = 0;
    }
    hlen[(len1 - 1) / 2] = 0;
    pow  = v + alloc;
    temp = pow + powlen;

    /* Let's start the actual work */
    
    for (i = 0, j = 0; i < len1 / 2; i++, j += 2)
    {
        if (poly1[j + 1] != 0L)
        {
            _nmod_vec_scalar_mul_nmod(h[i], poly2, len2, poly1[j + 1], mod);
            h[i][0] = n_addmod(h[i][0], poly1[j], mod.n);
            hlen[i] = len2;
        }
        else if (poly1[j] != 0L)
        {
            h[i][0] = poly1[j];
            hlen[i] = 1;
        }
    }
    if ((len1 & 1L))
    {
        if (poly1[j] != 0L)
        {
            h[i][0] = poly1[j];
            hlen[i] = 1;
        }
    }
    
    powlen = FLINT_MIN(N, 2 * len2 - 1);
    _nmod_poly_mullow(pow, poly2, len2, poly2, len2, powlen, mod);
    
    for (n = (len1 + 1) / 2; n > 2; n = (n + 1) / 2)
    {
        if (hlen[1] > 0)
        {
            len_t templen = FLINT_MIN(N, powlen + hlen[1] - 1);
            _nmod_poly_mullow(temp, pow, powlen, h[1], hlen[1], templen, mod);
            _nmod_poly_add(h[0], temp, templen, h[0], hlen[0], mod);
            hlen[0] = FLINT_MAX(hlen[0], templen);
        }
        
        for (i = 1; i < n / 2; i++)
        {
            if (hlen[2*i + 1] > 0)
            {
                hlen[i] = FLINT_MIN(N, hlen[2*i + 1] + powlen - 1);
                _nmod_poly_mullow(h[i], pow, powlen, h[2*i + 1], hlen[2*i + 1],
                    hlen[i], mod);
            }
            else
            {
                hlen[i] = 0;
            }
            _nmod_poly_add(h[i], h[i], hlen[i], h[2*i], hlen[2*i], mod);
            hlen[i] = FLINT_MAX(hlen[i], hlen[2*i]);
        }
        if ((n & 1L))
        {
            hlen[i] = FLINT_MIN(N, hlen[2*i]);
            flint_mpn_copyi(h[i], h[2*i], hlen[i]);
        }
        
        _nmod_poly_mullow(temp, pow, powlen, pow, powlen, 
                          FLINT_MIN(N, 2 * powlen - 1), mod);
        powlen = FLINT_MIN(N, 2 * powlen - 1);
        {
            mp_ptr t = pow;
            pow      = temp;
            temp     = t;
        }
    }

    _nmod_poly_mullow(res, pow, powlen, h[1], hlen[1], 
                      FLINT_MIN(N, powlen + hlen[1] - 1), mod);
    _nmod_vec_add(res, res, h[0], hlen[0], mod);
    
    _nmod_vec_clear(v);
    flint_free(h);
    flint_free(hlen);
}

void 
nmod_poly_compose_series_divconquer(nmod_poly_t res, 
    const nmod_poly_t poly1, const nmod_poly_t poly2, len_t N)
{
    const len_t len1 = poly1->length;
    const len_t len2 = FLINT_MIN(N, poly2->length);
    len_t lenr;
    
    if (len1 == 0 || N == 0)
    {
        nmod_poly_zero(res);
        return;
    }
    if (len1 == 1 || len2 == 0)
    {
        nmod_poly_set_coeff_ui(res, 0, poly1->coeffs[0]);
        nmod_poly_truncate(res, 1);
        return;
    }
    
    lenr = FLINT_MIN(N, (len1 - 1) * (len2 - 1) + 1);
    
    if (res != poly1 && res != poly2)
    {
        nmod_poly_fit_length(res, lenr);
        _nmod_poly_compose_series_divconquer(res->coeffs, 
            poly1->coeffs, len1, poly2->coeffs, len2, N, poly1->mod);
    }
    else
    {
        nmod_poly_t t;
        nmod_poly_init2(t, poly1->mod.n, lenr);
        _nmod_poly_compose_series_divconquer(t->coeffs, 
            poly1->coeffs, len1, poly2->coeffs, len2, N, poly1->mod);
        nmod_poly_swap(res, t);
        nmod_poly_clear(t);
    }

    res->length = lenr;
    _nmod_poly_normalise(res);
}
