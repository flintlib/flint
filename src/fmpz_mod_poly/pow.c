/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_pow(fmpz *res, const fmpz *poly, slong len, ulong e, 
                        const fmpz_t p)
{
    ulong bit = ~((~UWORD(0)) >> 1);
    slong rlen;
    slong alloc = (slong) e * (len - 1) + 1;
    fmpz *v = _fmpz_vec_init(alloc);
    fmpz *R, *S, *T;

    /*
       Set bits to the bitmask with a 1 one place lower than the msb of e
     */
    
    while ((bit & e) == UWORD(0))
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
            if ((bit2 & e) == UWORD(0))
                swaps = ~swaps;
        
        if (swaps == 0U)
        {
            R = res;
            S = v;
        }
        else
        {
            R = v;
            S = res;
        }
    }
    
    /*
       We unroll the first step of the loop, referring to {poly, len}
     */
    
    _fmpz_mod_poly_sqr(R, poly, len, p);
    rlen = 2 * len - 1;
    if ((bit & e))
    {
        _fmpz_mod_poly_mul(S, R, rlen, poly, len, p);
        rlen += len - 1;
        T = R;
        R = S;
        S = T;
    }
    
    while ((bit >>= 1))
    {
        if ((bit & e))
        {
            _fmpz_mod_poly_sqr(S, R, rlen, p);
            rlen += rlen - 1;
            _fmpz_mod_poly_mul(R, S, rlen, poly, len, p);
            rlen += len - 1;
        }
        else
        {
            _fmpz_mod_poly_sqr(S, R, rlen, p);
            rlen += rlen - 1;
            T = R;
            R = S;
            S = T;
        }
    }
    
    _fmpz_vec_clear(v, alloc);
}

void fmpz_mod_poly_pow(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, ulong e,
                                                      const fmpz_mod_ctx_t ctx)
{
    const slong len = op->length;
    slong rlen;

    if ((len < 2) || (e < UWORD(3)))
    {
        if (e == UWORD(0))
            fmpz_mod_poly_set_ui(rop, 1, ctx);
        else if (len == 0)
            fmpz_mod_poly_zero(rop, ctx);
        else if (len == 1)
        {
            fmpz_mod_poly_fit_length(rop, 1, ctx);
            fmpz_powm_ui(rop->coeffs, op->coeffs, e, fmpz_mod_ctx_modulus(ctx));
            _fmpz_mod_poly_set_length(rop, 1);
            _fmpz_mod_poly_normalise(rop);
        }
        else if (e == UWORD(1))
            fmpz_mod_poly_set(rop, op, ctx);
        else  /* e == UWORD(2) */
            fmpz_mod_poly_sqr(rop, op, ctx);
        return;
    }

    rlen = (slong) e * (len - 1) + 1;

    if (rop != op)
    {
        fmpz_mod_poly_fit_length(rop, rlen, ctx);
        _fmpz_mod_poly_pow(rop->coeffs, op->coeffs, len, e, fmpz_mod_ctx_modulus(ctx));
        _fmpz_mod_poly_set_length(rop, rlen);
    }
    else
    {
        fmpz *t = _fmpz_vec_init(rlen);

        _fmpz_mod_poly_pow(t, op->coeffs, len, e, fmpz_mod_ctx_modulus(ctx));

        _fmpz_vec_clear(rop->coeffs, rop->alloc);
        rop->coeffs = t;
        rop->alloc  = rlen;
        rop->length = rlen;
    }

    _fmpz_mod_poly_normalise(rop);
}

