/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void
_nmod_poly_pow_trunc_binexp(mp_ptr res, mp_srcptr poly, 
                                ulong e, slong trunc, nmod_t mod)
{
    ulong bit = ~((~UWORD(0)) >> 1);
    mp_ptr v = _nmod_vec_init(trunc);
    mp_ptr R, S, T;

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
    
    _nmod_poly_mullow(R, poly, trunc, poly, trunc, trunc, mod);
    if ((bit & e))
    {
        _nmod_poly_mullow(S, R, trunc, poly, trunc, trunc, mod);
        T = R;
        R = S;
        S = T;
    }
    
    while ((bit >>= 1))
    {
        if ((bit & e))
        {
            _nmod_poly_mullow(S, R, trunc, R, trunc, trunc, mod);
            _nmod_poly_mullow(R, S, trunc, poly, trunc, trunc, mod);
        }
        else
        {
            _nmod_poly_mullow(S, R, trunc, R, trunc, trunc, mod);
            T = R;
            R = S;
            S = T;
        }
    }
    
    _nmod_vec_clear(v);
}

void
nmod_poly_pow_trunc_binexp(nmod_poly_t res, 
                           const nmod_poly_t poly, ulong e, slong trunc)
{
    const slong len = poly->length;
    mp_ptr p;
    int pcopy = 0;

    if (len < 2 || e < UWORD(3) || trunc == 0)
    {
        if (len == 0 || trunc == 0)
            nmod_poly_zero(res);
        else if (len == 1)
        {
            nmod_poly_fit_length(res, 1);
            res->coeffs[0] = n_powmod2_ui_preinv(poly->coeffs[0], e,
                poly->mod.n, poly->mod.ninv);
            res->length = 1;
            _nmod_poly_normalise(res);
        }
        else if (e == UWORD(0))
        {
            nmod_poly_set_coeff_ui(res, 0, UWORD(1));
            res->length = 1;
            _nmod_poly_normalise(res);
        }
        else if (e == UWORD(1))
        {
            nmod_poly_set(res, poly);
            nmod_poly_truncate(res, trunc);
        }
        else  /* e == UWORD(2) */
            nmod_poly_mullow(res, poly, poly, trunc);

        return;
    }

    if (poly->length < trunc)
    {
        p = _nmod_vec_init(trunc);
        flint_mpn_copyi(p, poly->coeffs, poly->length);
        flint_mpn_zero(p + poly->length, trunc - poly->length);
        pcopy = 1;
    } else
        p = poly->coeffs;

    if (res != poly || pcopy)
    {
        nmod_poly_fit_length(res, trunc);
        _nmod_poly_pow_trunc_binexp(res->coeffs, p, e, trunc, poly->mod);
    }
    else
    {
        nmod_poly_t t;
        nmod_poly_init2(t, poly->mod.n, trunc);
        _nmod_poly_pow_trunc_binexp(t->coeffs, p, e, trunc, poly->mod);
        nmod_poly_swap(res, t);
        nmod_poly_clear(t);
    }

    if (pcopy)
        _nmod_vec_clear(p);

    res->length = trunc;
    _nmod_poly_normalise(res);
}
