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
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "ulong_extras.h"

void
_fmpq_div(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den,
                                    const fmpz_t op2num, const fmpz_t op2den)
{
    fmpz_t t, u;

    fmpz_init(t);
    fmpz_init(u);
    fmpz_set(t, op2den);
    fmpz_set(u, op2num);

    _fmpq_mul(rnum, rden, op1num, op1den, t, u);

    fmpz_clear(t);
    fmpz_clear(u);

    if (fmpz_sgn(rden) < 0)
    {
        fmpz_neg(rnum, rnum);
        fmpz_neg(rden, rden);
    }
}

void fmpq_div(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
{
    if (fmpq_is_zero(op2))
    {
        printf("Exception (fmpq_div). Division by zero.\n");
        abort();
    }

    _fmpq_div(fmpq_numref(res), fmpq_denref(res),
              fmpq_numref(op1), fmpq_denref(op1),
              fmpq_numref(op2), fmpq_denref(op2));
}

