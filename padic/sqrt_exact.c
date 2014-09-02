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

    Copyright (C) 2011 Jan Tuitman
    Copyright (C) 2011, 2012 Sebastian Pancratz
    Copyright (C) 2014 William Hart

******************************************************************************/

#include "padic.h"

int padic_sqrt_exact(padic_t rop, const padic_t op)
{
    fmpz_t r;
    int res = 1;

    if (padic_is_zero(op))
    {
        padic_zero(rop);
        return 1;
    }
    if (padic_val(op) & WORD(1))
    {
        return 0;
    }

    padic_val(rop) = padic_val(op) / 2;

    fmpz_init(r);

    fmpz_sqrtrem(padic_unit(rop), r, padic_unit(op));

    res = fmpz_is_zero(r);

    fmpz_clear(r);

    return res;
}

