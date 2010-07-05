/*============================================================================

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

===============================================================================*/
/****************************************************************************

    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2010 William Hart

*****************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"

void _fmpz_vec_content(fmpz_t res, const fmpz * vec, ulong len)
{
	ulong i;
    fmpz_t temp;
    
    if (len == 0)
    {
        fmpz_zero(res);
        return;
    }
    
    fmpz_init(temp);
    fmpz_abs(temp, vec);
    for (i = 1; (*temp != 1) && (i < len); i++)
        fmpz_gcd(temp, temp, vec + i);
    fmpz_swap(res, temp);
    fmpz_clear(temp);
}

