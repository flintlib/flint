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
_fmpq_next_calkin_wilf(fmpz_t rnum, fmpz_t rden,
    const fmpz_t num, const fmpz_t den)
{
    fmpz n, d;

    n = *num;
    d = *den;

    if (!COEFF_IS_MPZ(n) && !COEFF_IS_MPZ(d))
    {
        /* This does not overflow, as the larger part at most doubles */
        fmpz_set_ui(rnum, d);
        fmpz_set_ui(rden, d*(n / d) + d - (n % d));
    }
    else
    {
        fmpz_t q, r, t;

        fmpz_init(q);
        fmpz_init(r);
        fmpz_init(t);

        fmpz_fdiv_qr(q, r, num, den);
        fmpz_set(rnum, den);
        fmpz_mul(t, q, den);
        fmpz_add(rden, t, den);
        fmpz_sub(rden, rden, r);

        fmpz_clear(q);
        fmpz_clear(r);
        fmpz_clear(t);
    }
}

void
fmpq_next_calkin_wilf(fmpq_t res, const fmpq_t x)
{
    _fmpq_next_calkin_wilf(&res->num, &res->den, &x->num, &x->den);
}
