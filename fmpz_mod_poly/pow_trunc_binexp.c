/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2010 William Hart
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

void
_fmpz_mod_poly_pow_trunc_binexp(fmpz * res, const fmpz * poly,
                                ulong e, slong trunc, const fmpz_t p)
{
    ulong bit = ~((~UWORD(0)) >> 1);
    fmpz * v = _fmpz_vec_init(trunc);
    fmpz * R, * S, * T;

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

    _fmpz_mod_poly_mullow(R, poly, trunc, poly, trunc, p, trunc);
    if ((bit & e))
    {
        _fmpz_mod_poly_mullow(S, R, trunc, poly, trunc, p, trunc);
        T = R;
        R = S;
        S = T;
    }

    while ((bit >>= 1))
    {
        if ((bit & e))
        {
            _fmpz_mod_poly_mullow(S, R, trunc, R, trunc, p, trunc);
            _fmpz_mod_poly_mullow(R, S, trunc, poly, trunc, p, trunc);
        }
        else
        {
            _fmpz_mod_poly_mullow(S, R, trunc, R, trunc, p, trunc);
            T = R;
            R = S;
            S = T;
        }
    }

    _fmpz_vec_clear(v, trunc);
}

void
fmpz_mod_poly_pow_trunc_binexp(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly,
                                ulong e, slong trunc, const fmpz_mod_ctx_t ctx)
{
    const slong len = poly->length;
    fmpz * q;
    int qcopy = 0;

    if (len < 2 || e < UWORD(3) || trunc == 0)
    {
        if (len == 0 || trunc == 0)
        {
            fmpz_mod_poly_zero(res, ctx);
        }
        else if (len == 1)
        {
            fmpz_mod_poly_fit_length(res, 1, ctx);
            fmpz_powm_ui(res->coeffs, poly->coeffs, e, fmpz_mod_ctx_modulus(ctx));
            _fmpz_mod_poly_set_length(res, 1);
            _fmpz_mod_poly_normalise(res);
        }
        else if (e == UWORD(0))
        {
            fmpz_mod_poly_set_coeff_ui(res, 0, UWORD(1), ctx);
            _fmpz_mod_poly_set_length(res, 1);
            _fmpz_mod_poly_normalise(res);
        }
        else if (e == UWORD(1))
        {
            fmpz_mod_poly_set(res, poly, ctx);
            fmpz_mod_poly_truncate(res, trunc, ctx);
        }
        else  /* e == UWORD(2) */
            fmpz_mod_poly_mullow(res, poly, poly, trunc, ctx);

        return;
    }

    if (poly->length < trunc)
    {
        q = _fmpz_vec_init(trunc);
        _fmpz_vec_set(q, poly->coeffs, poly->length);
        _fmpz_vec_zero(q + poly->length, trunc - poly->length);
        qcopy = 1;
    } else
        q = poly->coeffs;

    if (res != poly || qcopy)
    {
        fmpz_mod_poly_fit_length(res, trunc, ctx);
        _fmpz_mod_poly_pow_trunc_binexp(res->coeffs, q, e, trunc,
                                                    fmpz_mod_ctx_modulus(ctx));
    }
    else
    {
        fmpz_mod_poly_t t;
        fmpz_mod_poly_init2(t, trunc, ctx);
        _fmpz_mod_poly_pow_trunc_binexp(t->coeffs, q, e, trunc,
                                                    fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_poly_swap(res, t, ctx);
        fmpz_mod_poly_clear(t, ctx);
    }

    if (qcopy)
        _fmpz_vec_clear(q, trunc);

    _fmpz_mod_poly_set_length(res, trunc);
    _fmpz_mod_poly_normalise(res);
}

