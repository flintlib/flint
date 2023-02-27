/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_compose_divconquer(arb_ptr res, arb_srcptr poly1, slong len1,
                                          arb_srcptr poly2, slong len2, slong prec)
{
    slong i, j, k, n;
    slong *hlen, alloc, powlen;
    arb_ptr v, pow, temp;
    arb_ptr * h;

    if (len1 == 1)
    {
        arb_set(res, poly1);
        return;
    }
    if (len2 == 1)
    {
        _arb_poly_evaluate(res, poly1, len1, poly2, prec);
        return;
    }
    if (len1 == 2)
    {
        _arb_poly_compose_horner(res, poly1, len1, poly2, len2, prec);
        return;
    }

    /* Initialisation */
    
    hlen = (slong *) flint_malloc(((len1 + 1) / 2) * sizeof(slong));
    
    for (k = 1; (2 << k) < len1; k++) ;
    
    hlen[0] = hlen[1] = ((1 << k) - 1) * (len2 - 1) + 1;
    for (i = k - 1; i > 0; i--)
    {
        slong hi = (len1 + (1 << i) - 1) / (1 << i);
        for (n = (hi + 1) / 2; n < hi; n++)
            hlen[n] = ((1 << i) - 1) * (len2 - 1) + 1;
    }
    powlen = (1 << k) * (len2 - 1) + 1;
    
    alloc = 0;
    for (i = 0; i < (len1 + 1) / 2; i++)
        alloc += hlen[i];

    v = _arb_vec_init(alloc + 2 * powlen);
    h = (arb_ptr *) flint_malloc(((len1 + 1) / 2) * sizeof(arb_ptr));
    h[0] = v;
    for (i = 0; i < (len1 - 1) / 2; i++)
    {
        h[i + 1] = h[i] + hlen[i];
        hlen[i] = 0;
    }
    hlen[(len1 - 1) / 2] = 0;
    pow = v + alloc;
    temp = pow + powlen;
    
    /* Let's start the actual work */
    
    for (i = 0, j = 0; i < len1 / 2; i++, j += 2)
    {
        if (!arb_is_zero(poly1 + j + 1))
        {
            _arb_vec_scalar_mul(h[i], poly2, len2, poly1 + j + 1, prec);
            arb_add(h[i], h[i], poly1 + j, prec);
            hlen[i] = len2;
        }
        else if (!arb_is_zero(poly1 + j))
        {
            arb_set(h[i], poly1 + j);
            hlen[i] = 1;
        }
    }
    if ((len1 & WORD(1)))
    {
        if (!arb_is_zero(poly1 + j))
        {
            arb_set(h[i], poly1 + j);
            hlen[i] = 1;
        }
    }
    
    _arb_poly_mul(pow, poly2, len2, poly2, len2, prec);
    powlen = 2 * len2 - 1;
    
    for (n = (len1 + 1) / 2; n > 2; n = (n + 1) / 2)
    {
        if (hlen[1] > 0)
        {
            slong templen = powlen + hlen[1] - 1;
            _arb_poly_mul(temp, pow, powlen, h[1], hlen[1], prec);
            _arb_poly_add(h[0], temp, templen, h[0], hlen[0], prec);
            hlen[0] = FLINT_MAX(hlen[0], templen);
        }
        
        for (i = 1; i < n / 2; i++)
        {
            if (hlen[2*i + 1] > 0)
            {
                _arb_poly_mul(h[i], pow, powlen, h[2*i + 1], hlen[2*i + 1], prec);
                hlen[i] = hlen[2*i + 1] + powlen - 1;
            } else
                hlen[i] = 0;
            _arb_poly_add(h[i], h[i], hlen[i], h[2*i], hlen[2*i], prec);
            hlen[i] = FLINT_MAX(hlen[i], hlen[2*i]);
        }
        if ((n & WORD(1)))
        {
            _arb_vec_set(h[i], h[2*i], hlen[2*i]);
            hlen[i] = hlen[2*i];
        }
        
        _arb_poly_mul(temp, pow, powlen, pow, powlen, prec);
        powlen += powlen - 1;
        {
            arb_ptr t = pow;
            pow = temp;
            temp = t;
        }
    }

    _arb_poly_mul(res, pow, powlen, h[1], hlen[1], prec);
    _arb_vec_add(res, res, h[0], hlen[0], prec);
    
    _arb_vec_clear(v, alloc + 2 * powlen);
    flint_free(h);
    flint_free(hlen);
}

void
arb_poly_compose_divconquer(arb_poly_t res,
                             const arb_poly_t poly1, const arb_poly_t poly2, slong prec)
{
    const slong len1 = poly1->length;
    const slong len2 = poly2->length;
    
    if (len1 == 0)
    {
        arb_poly_zero(res);
    }
    else if (len1 == 1 || len2 == 0)
    {
        arb_poly_set_arb(res, poly1->coeffs);
    }
    else
    {
        const slong lenr = (len1 - 1) * (len2 - 1) + 1;
        
        if (res != poly1 && res != poly2)
        {
            arb_poly_fit_length(res, lenr);
            _arb_poly_compose_divconquer(res->coeffs, poly1->coeffs, len1,
                                                   poly2->coeffs, len2, prec);
        }
        else
        {
            arb_poly_t t;
            arb_poly_init2(t, lenr);
            _arb_poly_compose_divconquer(t->coeffs, poly1->coeffs, len1,
                                                 poly2->coeffs, len2, prec);
            arb_poly_swap(res, t);
            arb_poly_clear(t);
        }

        _arb_poly_set_length(res, lenr);
        _arb_poly_normalise(res);
    }
}
