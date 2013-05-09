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

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"

void
_fmpz_vec_prod(fmpz_t res, const fmpz * vec, len_t len)
{
    if (len <= 1)
    {
        if (len == 1)
            fmpz_set(res, vec);
        else
            fmpz_one(res);
    }
    else if (len <= 3)
    {
        len_t i;

        fmpz_mul(res, vec, vec + 1);
        for (i = 2; i < len; i++)
            fmpz_mul(res, res, vec + i);
    }
    else
    {
        len_t m = len / 2;
        fmpz_t tmp;
        fmpz_init(tmp);
        _fmpz_vec_prod(res, vec, m);
        _fmpz_vec_prod(tmp, vec + m, len - m);
        fmpz_mul(res, res, tmp);
        fmpz_clear(tmp);
    }
}
