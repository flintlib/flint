/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

ulong
_nmod_poly_evaluate_nmod(ulong_srcptr poly, slong len, ulong c, nmod_t mod)
{
    slong m;
    ulong val;

    if (len == 0)
        return 0;

    if (len == 1 || c == 0)
        return poly[0];

    m = len - 1;
    
    val = poly[m];
    m--;

    for ( ; m >= 0; m--)
    {
        val = n_mulmod2_preinv(val, c, mod.n, mod.ninv);
        val = n_addmod(val, poly[m], mod.n);
    }

    return val;
}

ulong
nmod_poly_evaluate_nmod(const nmod_poly_t poly, ulong c)
{
    return _nmod_poly_evaluate_nmod(poly->coeffs, poly->length, c, poly->mod);
}

