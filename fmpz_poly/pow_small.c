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

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_pow_small(fmpz * res, const fmpz * poly, long len, ulong e)
{
    switch (e)
    {
        case 0:
            fmpz_one(res);
            break;
        case 1:
            _fmpz_vec_set(res, poly, len);
            break;
        case 2:
            _fmpz_poly_sqr(res, poly, len);
            break;
        case 3:
        {
            long alloc = 2 * len - 1;
            fmpz *t = _fmpz_vec_init(alloc);
            _fmpz_poly_sqr(t, poly, len);
            _fmpz_poly_mul(res, t, alloc, poly, len);
            _fmpz_vec_clear(t, alloc);
            break;
        }
        case 4:
        {
            long alloc = 2 * len - 1;
            fmpz *t = _fmpz_vec_init(alloc);
            _fmpz_poly_sqr(t, poly, len);
            _fmpz_poly_sqr(res, t, alloc);
            _fmpz_vec_clear(t, alloc);
            break;
        }
    }
}

