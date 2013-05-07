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
_fmpq_randbits(fmpz_t num, fmpz_t den, flint_rand_t state, mp_bitcnt_t bits)
{
    fmpz_randbits(num, state, bits);

    do {
        fmpz_randbits(den, state, bits);
    } while (fmpz_is_zero(den));

    _fmpq_canonicalise(num, den);
}

void fmpq_randbits(fmpq_t res, flint_rand_t state, mp_bitcnt_t bits)
{
    _fmpq_randbits(&res->num, &res->den, state, bits);
}
