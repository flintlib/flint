/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void
_nmod_poly_pow(mp_ptr res, mp_srcptr poly, slong len, ulong e, nmod_t mod)
{
    _nmod_poly_pow_binexp(res, poly, len, e, mod);
}

void
nmod_poly_pow(nmod_poly_t res, const nmod_poly_t poly, ulong e)
{
    const slong len = poly->length;
    slong rlen;

    if ((len < 2) | (e < UWORD(3)))
    {
        if (len == 0)
        {
            if (e == 0)
                nmod_poly_one(res);
            else
                nmod_poly_zero(res);
        }
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
            nmod_poly_set(res, poly);
        else  /* e == UWORD(2) */
            nmod_poly_mul(res, poly, poly);

        return;
    }

    rlen = (slong) e * (len - 1) + 1;

    if (res != poly)
    {
        nmod_poly_fit_length(res, rlen);
        _nmod_poly_pow(res->coeffs, poly->coeffs, len, e, poly->mod);
    }
    else
    {
        nmod_poly_t t;
        nmod_poly_init2(t, poly->mod.n, rlen);
        _nmod_poly_pow(t->coeffs, poly->coeffs, len, e, poly->mod);
        nmod_poly_swap(res, t);
        nmod_poly_clear(t);
    }

    res->length = rlen;
    _nmod_poly_normalise(res);
}
