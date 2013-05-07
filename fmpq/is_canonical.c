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
#include "fmpq.h"
#include "ulong_extras.h"


int _fmpq_is_canonical(const fmpz_t num, const fmpz_t den)
{
    fmpz_t u;
    int result;

    if (fmpz_is_one(den)) return 1;
    if (fmpz_sgn(den) <= 0) return 0;
    if (fmpz_is_zero(num)) return fmpz_is_one(den);

    fmpz_init(u);
    fmpz_gcd(u, num, den);
    result = fmpz_is_one(u);
    fmpz_clear(u);
    return result;
}

int fmpq_is_canonical(const fmpq_t x)
{
    return _fmpq_is_canonical(&x->num, &x->den);
}
