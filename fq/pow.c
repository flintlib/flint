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

    Copyright (C) 2011, 2012 Sebastian Pancratz 
 
******************************************************************************/

#include "fq.h"

void _fq_pow(fmpz *rop, const fmpz *op, long len, const fmpz_t e, 
             const fmpz *a, const long *j, long lena, const fmpz_t p)
{
    const long d = j[lena - 1];

    if (fmpz_is_zero(e))
    {
        fmpz_one(rop);
        _fmpz_vec_zero(rop + 1, 2 * d - 1 - 1);
    }
    else if (fmpz_is_one(e))
    {
        _fmpz_vec_set(rop, op, len);
        _fmpz_vec_zero(rop + len, 2 * d - 1 - len);
    }
    else
    {
        ulong bit;
        fmpz *v = _fmpz_vec_init(2 * d - 1);
        fmpz *R, *S, *T;

        _fmpz_vec_zero(rop, 2 * d - 1);

        /*
           Set bits to the bitmask with a 1 one place lower than the msb of e
         */

        bit = fmpz_bits(e) - 2;

        /*
           Trial run without any polynomial arithmetic to determine the parity 
           of the number of swaps;  then set R and S accordingly
         */
        
        {
            unsigned int swaps = 0U;
            ulong bit2 = bit;
            if (fmpz_tstbit(e, bit2))
                swaps = ~swaps;
            while (bit2--)
                if (!fmpz_tstbit(e, bit2))
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
           We unroll the first step of the loop, referring to {op, len}
         */

        _fmpz_poly_sqr(R, op, len);
        _fq_reduce(R, 2 * len - 1, a, j, lena, p);

        if (fmpz_tstbit(e, bit))
        {
            _fmpz_poly_mul(S, R, d, op, len);
            _fq_reduce(S, d + len - 1, a, j, lena, p);
            T = R;
            R = S;
            S = T;
        }

        while (bit--)
        {
            if (fmpz_tstbit(e, bit))
            {
                _fmpz_poly_sqr(S, R, d);
                _fq_reduce(S, 2 * d - 1, a, j, lena, p);
                _fmpz_poly_mul(R, S, d, op, len);
                _fq_reduce(R, d + len - 1, a, j, lena, p);
            }
            else
            {
                _fmpz_poly_sqr(S, R, d);
                _fq_reduce(S, 2 * d - 1, a, j, lena, p);
                T = R;
                R = S;
                S = T;
            }
        }

        _fmpz_vec_clear(v, 2 * d - 1);
    }
}

void fq_pow(fq_t rop, const fq_t op, const fmpz_t e, const fq_ctx_t ctx)
{
    if (fmpz_sgn(e) < 0)
    {
        printf("Exception (fq_pow).  e < 0.\n");
        abort();
    }

    if (fmpz_is_zero(e))
    {
        fq_one(rop);
    }
    else if (fq_is_zero(op))
    {
        fq_zero(rop);
    }
    else if (fmpz_is_one(e))
    {
        fq_set(rop, op);
    }
    else
    {
        const long d = fq_ctx_degree(ctx);
        fmpz *t;

        if (rop == op)
        {
            t = _fmpz_vec_init(2 * d - 1);
        }
        else
        {
            fmpz_poly_fit_length(rop, 2 * d - 1);
            t = rop->coeffs;
        }

        _fq_pow(t, op->coeffs, op->length, e, ctx->a, ctx->j, ctx->len, 
                fq_ctx_prime(ctx));

        if (rop == op)
        {
            _fmpz_vec_clear(rop->coeffs, rop->alloc);
            rop->coeffs = t;
            rop->alloc  = 2 * d - 1;
            rop->length = d;
        }
        else
        {
            _fmpz_poly_set_length(rop, d);
        }
        _fmpz_poly_normalise(rop);
    }
}

