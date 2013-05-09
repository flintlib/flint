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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_factor.h"


void
fmpz_factor_expand_iterative(fmpz_t n, const fmpz_factor_t factor)
{
    len_t i;
    fmpz_t tmp;

    fmpz_set_si(n, factor->sign);

    fmpz_init(tmp);
    for (i = 0; i < factor->num; i++)
    {
        fmpz_pow_ui(tmp, factor->p + i, factor->exp[i]);
        fmpz_mul(n, n, tmp);
    }
    fmpz_clear(tmp);
}
