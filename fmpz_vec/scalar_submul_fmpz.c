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

    Copyright (C) 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz
   
******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"

void
_fmpz_vec_scalar_submul_fmpz(fmpz * vec1, const fmpz * vec2, long len2,
                             const fmpz_t x)
{
    fmpz c = *x;

    if (!COEFF_IS_MPZ(c))
    {
        if (c == 0)
            return;
        else if (c == 1)
            _fmpz_vec_sub(vec1, vec1, vec2, len2);
        else if (c == -1)
            _fmpz_vec_add(vec1, vec1, vec2, len2);
        else
            _fmpz_vec_scalar_submul_si(vec1, vec2, len2, c);
    }
    else
    {
        long i;
        for (i = 0; i < len2; i++)
            fmpz_submul(vec1 + i, vec2 + i, x);
    }
}
