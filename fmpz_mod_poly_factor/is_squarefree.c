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

    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

******************************************************************************/

#include "fmpz_mod_poly_factor.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

#include <string.h>

int
_fmpz_mod_poly_is_squarefree(const fmpz * f, len_t len, const fmpz_t p)
{
    fmpz * fd, * g;
    fmpz_t invd;
    len_t dlen;
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

int fmpz_mod_poly_is_squarefree(const fmpz_mod_poly_t f)
{
    return _fmpz_mod_poly_is_squarefree(f->coeffs, f->length, &f->p);
}
