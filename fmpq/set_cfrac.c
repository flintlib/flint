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

void
_fmpq_set_cfrac(fmpz_t p, fmpz_t q, const fmpz * c, long n)
{
    fmpz_t t, u;
    long i;

    fmpz_set(p, c);
    fmpz_set_ui(q, 1);

    fmpz_init(t);
    fmpz_init(u);
    fmpz_set_ui(t, 1);

    for (i = 1; i < n; i++)
    {
        fmpz_addmul(t, c + i, p);
        fmpz_addmul(u, c + i, q);
        fmpz_swap(t, p);
        fmpz_swap(u, q);
    }

    fmpz_clear(t);
    fmpz_clear(u);
}

void
fmpq_set_cfrac(fmpq_t x, const fmpz * c, long n)
{
    _fmpq_set_cfrac(fmpq_numref(x), fmpq_denref(x), c, n);
}

