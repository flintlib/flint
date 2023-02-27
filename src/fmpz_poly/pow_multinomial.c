/*
    Copyright (C) 2010 Sebastian Pancratz

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

void
_fmpz_poly_pow_multinomial(fmpz * res, const fmpz * poly, slong len, ulong e)
{
    slong k, low, rlen;
    fmpz_t d, t;
    fmpz * P;
    
    rlen = (slong) e * (len - WORD(1)) + WORD(1);
    _fmpz_vec_zero(res, rlen);
    
    for (low = WORD(0); poly[low] == WORD(0); low++) ;
    if (low == WORD(0))
    {
        P = (fmpz *) poly;
    }
    else
    {
        P = (fmpz *) poly + low;
        len  -= low;
        res  += (slong) e * low;
        rlen -= (slong) e * low;
    }
    
    fmpz_init(d);
    fmpz_init(t);
    
    fmpz_pow_ui(res, P, e);
    
    for (k = 1; k < rlen; k++)
    {
        slong i, u = -k;
        for (i = 1; i <= FLINT_MIN(k, len - 1); i++)
        {
            fmpz_mul(t, P + i, res + (k - i));
            u += (slong) e + 1;
            if (u >= 0)
                fmpz_addmul_ui(res + k, t, (ulong) u);
            else
                fmpz_submul_ui(res + k, t, - ((ulong) u));
        }
        fmpz_add(d, d, P);
        fmpz_divexact(res + k, res + k, d);
    }
    
    fmpz_clear(d);
    fmpz_clear(t);
}

void
fmpz_poly_pow_multinomial(fmpz_poly_t res, const fmpz_poly_t poly, ulong e)
{
    const slong len = poly->length;
    slong rlen;

    if ((len < 2) | (e < UWORD(3)))
    {
        if (e == UWORD(0))
            fmpz_poly_set_ui(res, 1);
        else if (len == 0)
            fmpz_poly_zero(res);
        else if (len == 1)
        {
            fmpz_poly_fit_length(res, 1);
            fmpz_pow_ui(res->coeffs, poly->coeffs, e);
            _fmpz_poly_set_length(res, 1);
        }
        else if (e == UWORD(1))
            fmpz_poly_set(res, poly);
        else  /* e == UWORD(2) */
            fmpz_poly_sqr(res, poly);
        return;
    }
    
    rlen = (slong) e * (len - 1) + 1;

    if (res != poly)
    {
        fmpz_poly_fit_length(res, rlen);
        _fmpz_poly_pow_multinomial(res->coeffs, poly->coeffs, len, e);
        _fmpz_poly_set_length(res, rlen);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, rlen);
        _fmpz_poly_pow_multinomial(t->coeffs, poly->coeffs, len, e);
        _fmpz_poly_set_length(t, rlen);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
}
