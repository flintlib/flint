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


void
_fmpq_next_signed_minimal(fmpz_t rnum, fmpz_t rden,
    const fmpz_t num, const fmpz_t den)
{
    if (fmpz_sgn(num) > 0)
    {
        fmpz_neg(rnum, num);
        fmpz_set(rden, den);
    }
    else
    {
        fmpz_neg(rnum, num);
        _fmpq_next_minimal(rnum, rden, rnum, den);
    }
}

void
fmpq_next_signed_minimal(fmpq_t res, const fmpq_t x)
{
    _fmpq_next_signed_minimal(&res->num, &res->den, &x->num, &x->den);
}
