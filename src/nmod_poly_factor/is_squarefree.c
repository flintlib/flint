/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "mpn_extras.h"
#include "ulong_extras.h"

int
_nmod_poly_is_squarefree(mp_srcptr f, slong len, nmod_t mod)
{
    mp_ptr fd, g;
    slong dlen;
    int res;

    if (len <= 2)
        return len != 0;

    fd = flint_malloc(sizeof(mp_limb_t) * 2 * (len - 1));
    g = fd + len - 1;

    _nmod_poly_derivative(fd, f, len, mod);
    dlen = len - 1;
    MPN_NORM(fd, dlen);

    if (dlen)
        res = (_nmod_poly_gcd(g, f, len, fd, dlen, mod) == 1);
    else
        res = 0;   /* gcd(f, 0) = f, and len(f) > 2 */

    flint_free(fd);
    return res;
}

int
nmod_poly_is_squarefree(const nmod_poly_t f)
{
    return _nmod_poly_is_squarefree(f->coeffs, f->length, f->mod);
}
