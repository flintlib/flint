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

    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include "fq_poly.h"

void _fq_poly_compose_divconquer(fq_struct *rop, 
                                 const fq_struct *op1, long len1, 
                                 const fq_struct *op2, long len2, 
                                 const fq_ctx_t ctx)
{
    long i, j, k, n;
    long *hlen, alloc, powlen;
    fq_struct *v, **h, *pow, *temp;

    if (len1 <= 2 || len2 <= 1)
    {
        if (len1 == 1)
            fq_set(rop, op1);
        else if (len2 == 1)
            _fq_poly_evaluate_fq(rop, op1, len1, op2, ctx);
        else  /* len1 == 2 */
            _fq_poly_compose_horner(rop, op1, len1, op2, len2, ctx);
        return;
    }

    /* Initialisation */
    
    hlen = (long *) flint_malloc(((len1 + 1) / 2) * sizeof(long));

    k = FLINT_CLOG2(len1) - 1;

    hlen[0] = hlen[1] = ((1 << k) - 1) * (len2 - 1) + 1;
    for (i = k - 1; i > 0; i--)
    {
        long hi = (len1 + (1 << i) - 1) / (1 << i);
        for (n = (hi + 1) / 2; n < hi; n++)
            hlen[n] = ((1 << i) - 1) * (len2 - 1) + 1;
    }
    powlen = (1 << k) * (len2 - 1) + 1;
    
    alloc = 0;
    for (i = 0; i < (len1 + 1) / 2; i++)
        alloc += hlen[i];

    v = _fq_poly_init(alloc + 2 * powlen);
    h = (fq_struct **) flint_malloc(((len1 + 1) / 2) * sizeof(fq_struct *));
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
        if (!fq_is_zero(op1 + (j + 1)))
        {
            _fq_poly_scalar_mul_fq(h[i], op2, len2, op1 + j + 1, ctx);
            fq_add(h[i], h[i], op1 + j, ctx);
            hlen[i] = len2;
        }
        else if (!fq_is_zero(op1 + j))
        {
            fq_set(h[i], op1 + j);
            hlen[i] = 1;
        }
    }
    if ((len1 & 1L))
    {
        if (!fq_is_zero(op1 + j))
        {
            fq_set(h[i], op1 + j);
            hlen[i] = 1;
        }
    }
    
    _fq_poly_sqr(pow, op2, len2, ctx);
    powlen = 2 * len2 - 1;
    
    for (n = (len1 + 1) / 2; n > 2; n = (n + 1) / 2)
    {
        if (hlen[1] > 0)
        {
            long templen = powlen + hlen[1] - 1;
            _fq_poly_mul(temp, pow, powlen, h[1], hlen[1], ctx);
            _fq_poly_add(h[0], temp, templen, h[0], hlen[0], ctx);
            hlen[0] = FLINT_MAX(hlen[0], templen);
        }
        
        for (i = 1; i < n / 2; i++)
        {
            if (hlen[2*i + 1] > 0)
            {
                _fq_poly_mul(h[i], pow, powlen, h[2*i + 1], hlen[2*i + 1], ctx);
                hlen[i] = hlen[2*i + 1] + powlen - 1;
            } else
                hlen[i] = 0;
            _fq_poly_add(h[i], h[i], hlen[i], h[2*i], hlen[2*i], ctx);
            hlen[i] = FLINT_MAX(hlen[i], hlen[2*i]);
        }
        if ((n & 1L))
        {
            _fq_poly_set(h[i], h[2*i], hlen[2*i]);
            hlen[i] = hlen[2*i];
        }
        
        _fq_poly_sqr(temp, pow, powlen, ctx);
        powlen += powlen - 1;
        {
            fq_struct *t = pow;
            pow          = temp;
            temp         = t;
        }
    }

    _fq_poly_mul(rop, pow, powlen, h[1], hlen[1], ctx);
    _fq_poly_add(rop, rop, hlen[0], h[0], hlen[0], ctx);
    
    _fq_poly_clear(v, alloc + 2 * powlen);
    flint_free(h);
    flint_free(hlen);
}

void fq_poly_compose_divconquer(fq_poly_t rop, 
                                const fq_poly_t op1, const fq_poly_t op2, 
                                const fq_ctx_t ctx)
{
    const long len1 = op1->length;
    const long len2 = op2->length;
    const long lenr = (len1 - 1) * (len2 - 1) + 1;
    
    if (len1 == 0)
    {
        fq_poly_zero(rop);
    }
    else if (len1 == 1 || len2 == 0)
    {
        fq_poly_set_fq(rop, op1->coeffs + 0);
    }
    else if (rop != op1 && rop != op2)
    {
        fq_poly_fit_length(rop, lenr);
        _fq_poly_compose_divconquer(rop->coeffs, op1->coeffs, len1, 
                                                 op2->coeffs, len2, ctx);
        _fq_poly_set_length(rop, lenr);
        _fq_poly_normalise(rop);
    }
    else
    {
        fq_poly_t t;

        fq_poly_init2(t, lenr);
        _fq_poly_compose_divconquer(t->coeffs, op1->coeffs, len1,
                                               op2->coeffs, len2, ctx);
        _fq_poly_set_length(t, lenr);
        _fq_poly_normalise(t);
        fq_poly_swap(rop, t);
        fq_poly_clear(t);
    }
}

