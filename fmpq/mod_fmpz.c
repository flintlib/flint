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


int
_fmpq_mod_fmpz(fmpz_t res, const fmpz_t num, const fmpz_t den, const fmpz_t mod)
{
    int result;
    fmpz_t tmp;

    fmpz_init(tmp);
    result = fmpz_invmod(tmp, den, mod);
    fmpz_mul(tmp, tmp, num);
    fmpz_mod(res, tmp, mod);
    fmpz_clear(tmp);

    return result;
}

int fmpq_mod_fmpz(fmpz_t res, const fmpq_t x, const fmpz_t mod)
{
    return _fmpq_mod_fmpz(res, &x->num, &x->den, mod);
}

