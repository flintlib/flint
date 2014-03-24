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

#include "arith.h"
#include "fmpq.h"

void _arith_harmonic_number(fmpz_t num, fmpz_t den, slong n)
{
    if (n <= 0)
    {
        fmpz_zero(num);
        fmpz_one(den);
    }
    else
    {
        _fmpq_harmonic_ui(num, den, n);
    }
}

void arith_harmonic_number(fmpq_t x, slong n)
{
    _arith_harmonic_number(fmpq_numref(x), fmpq_denref(x), n);
}
