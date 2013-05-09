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

#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "ulong_extras.h"

void _fmpq_set_si(fmpz_t rnum, fmpz_t rden, len_t p, ulong q)
{
    if (q == 1 || p == 0)
    {
        fmpz_set_si(rnum, p);
        fmpz_one(rden);
    }
    else
    {
        ulong r = n_gcd_full(p < 0 ? (-(ulong) p) : (ulong) p, q);

        if (p < 0)
        {
            fmpz_set_ui(rnum, (-(ulong) p) / r);
            fmpz_neg(rnum, rnum);
        }
        else
            fmpz_set_si(rnum, p / r);

        fmpz_set_ui(rden, q / r);
    }
}

void fmpq_set_si(fmpq_t res, len_t p, ulong q)
{
    _fmpq_set_si(&res->num, &res->den, p, q);
}
