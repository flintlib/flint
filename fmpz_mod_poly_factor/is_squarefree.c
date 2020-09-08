/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

#include <string.h>

int
_fmpz_mod_poly_is_squarefree(const fmpz * f, slong len, const fmpz_t p)
{
    fmpz * fd, * g;
    fmpz_t invd;
    slong dlen;
    int res;

    if (len <= 2)
        return len != 0;

    fd = _fmpz_vec_init(2 * (len - 1));
    g = fd + len - 1;

    _fmpz_mod_poly_derivative(fd, f, len, p);
    dlen = len - 1;
    FMPZ_VEC_NORM(fd, dlen);

    if (dlen)
    {
        fmpz_init(invd);
        fmpz_invmod(invd, fd + dlen - 1, p);
        res = (_fmpz_mod_poly_gcd(g, f, len, fd, dlen, invd, p) == 1);
        fmpz_clear(invd);
    }
    else
        res = 0;   /* gcd(f, 0) = f, and len(f) > 2 */

    _fmpz_vec_clear(fd, 2 * (len - 1));
    return res;
}

int fmpz_mod_poly_is_squarefree(const fmpz_mod_poly_t f,
                                                      const fmpz_mod_ctx_t ctx)
{
    return _fmpz_mod_poly_is_squarefree(f->coeffs, f->length,
                                                    fmpz_mod_ctx_modulus(ctx));
}
