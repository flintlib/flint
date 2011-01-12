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

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "ulong_extras.h"


void
_fmpq_mul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den,
            const fmpz_t op2num, const fmpz_t op2den)
{
    fmpz_mul(rnum, op1num, op2num);
    fmpz_mul(rden, op1den, op2den);
    _fmpq_canonicalise(rnum, rden);
}

void fmpq_mul(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
{
    _fmpq_mul(&res->num, &res->den, &op1->num, &op1->den,
            &op2->num, &op2->den);
}
