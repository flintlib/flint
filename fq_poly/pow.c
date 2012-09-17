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

    Copyright (C) 2010, 2012 Sebastian Pancratz

******************************************************************************/

#include "fq_poly.h"

static fq_struct * __vec_init(long len)
{
    long i;
    fq_struct *v;

    v = flint_malloc(len * sizeof(fq_struct));
    for (i = 0; i < len; i++)
        fq_init(v + i);
    return v;
}

static void __vec_clear(fq_struct *v, long len)
{
    long i;

    for (i = 0; i < len; i++)
        fq_clear(v + i);
    flint_free(v);
}

void _fq_poly_pow(fq_struct *rop, const fq_struct *op, long len, ulong e, 
                                  const fq_ctx_t ctx)
{
    ulong bit = ~((~0UL) >> 1);
    long rlen;
    long alloc = (long) e * (len - 1) + 1;
    fq_struct *v = __vec_init(alloc);
    fq_struct *R, *S, *T;

    /*
       Set bits to the bitmask with a 1 one place lower than the msb of e
     */
    
    while ((bit & e) == 0UL)
        bit >>= 1;
    
    bit >>= 1;
    
    /*
       Trial run without any polynomial arithmetic to determine the parity 
       of the number of swaps;  then set R and S accordingly
     */
    
    {
        unsigned int swaps = 0U;
        ulong bit2 = bit;
        if ((bit2 & e))
            swaps = ~swaps;
        while (bit2 >>= 1)
            if ((bit2 & e) == 0UL)
                swaps = ~swaps;
        
        if (swaps == 0U)
        {
            R = rop;
            S = v;
        }
        else
        {
            R = v;
            S = rop;
        }
    }
    
    /*
       We unroll the first step of the loop, referring to {poly, len}
     */
    
    _fq_poly_sqr(R, op, len, ctx);
    rlen = 2 * len - 1;
    if ((bit & e))
    {
        _fq_poly_mul(S, R, rlen, op, len, ctx);
        rlen += len - 1;
        T = R;
        R = S;
        S = T;
    }
    
    while ((bit >>= 1))
    {
        if ((bit & e))
        {
            _fq_poly_sqr(S, R, rlen, ctx);
            rlen += rlen - 1;
            _fq_poly_mul(R, S, rlen, op, len, ctx);
            rlen += len - 1;
        }
        else
        {
            _fq_poly_sqr(S, R, rlen, ctx);
            rlen += rlen - 1;
            T = R;
            R = S;
            S = T;
        }
    }
    
    __vec_clear(v, alloc);
}

void fq_poly_pow(fq_poly_t rop, const fq_poly_t op, ulong e, 
                                const fq_ctx_t ctx)
{
    const long len = op->length;

    if ((len < 2) | (e < 3UL))
    {
        if (e == 0UL)
            fq_poly_one(rop);
        else if (len == 0)
            fq_poly_zero(rop);
        else if (len == 1)
        {
            fmpz_t f;
            fmpz_init_set_ui(f, e);

            fq_poly_fit_length(rop, 1);
            fq_pow(rop->coeffs + 0, op->coeffs + 0, f, ctx);
            _fq_poly_set_length(rop, 1);

            fmpz_clear(f);
        }
        else if (e == 1UL)
            fq_poly_set(rop, op);
        else  /* e == 2UL */
            fq_poly_sqr(rop, op, ctx);
    }
    else
    {
        const long rlen = (long) e * (len - 1) + 1;

        if (rop != op)
        {
            fq_poly_fit_length(rop, rlen);
            _fq_poly_pow(rop->coeffs, op->coeffs, len, e, ctx);
            _fq_poly_set_length(rop, rlen);
        }
        else
        {
            fq_poly_t t;
            fq_poly_init2(t, rlen);
            _fq_poly_pow(t->coeffs, op->coeffs, len, e, ctx);
            _fq_poly_set_length(t, rlen);
            fq_poly_swap(rop, t);
            fq_poly_clear(t);
        }
    }
}

