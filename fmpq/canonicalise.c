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


void _fmpq_canonicalise(fmpz_t num, fmpz_t den)
{
    fmpz_t u;

    if (*den == 1L)
        return;

    if (*num == 0L)
    {
        fmpz_one(den);
        return;
    }

    fmpz_init(u);
    fmpz_gcd(u, num, den);

    if (!fmpz_is_one(u))
    {
        fmpz_divexact(num, num, u);
        fmpz_divexact(den, den, u);
    }

    fmpz_clear(u);

    if (fmpz_sgn(den) < 0)
    {
        fmpz_neg(num, num);
        fmpz_neg(den, den);
    }
}

void fmpq_canonicalise(fmpq_t res)
{
    _fmpq_canonicalise(&res->num, &res->den);
}
