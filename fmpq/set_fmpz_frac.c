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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "ulong_extras.h"

void
fmpq_set_fmpz_frac(fmpq_t res, const fmpz_t p, const fmpz_t q)
{
    if (fmpz_is_zero(p))
    {
        fmpq_zero(res);
    }
    else if (fmpz_is_pm1(q) || fmpz_is_pm1(p))
    {
        if (fmpz_sgn(q) < 0)
        {
            fmpz_neg(fmpq_numref(res), p);
            fmpz_neg(fmpq_denref(res), q);
        }
        else
        {
            fmpz_set(fmpq_numref(res), p);
            fmpz_set(fmpq_denref(res), q);
        }
    }
    else
    {
        fmpz_t t;

        fmpz_init(t);
        fmpz_gcd(t, p, q);

        if (fmpz_is_one(t))
        {
            fmpz_set(fmpq_numref(res), p);
            fmpz_set(fmpq_denref(res), q);
        }
        else
        {
            fmpz_divexact(fmpq_numref(res), p, t);
            fmpz_divexact(fmpq_denref(res), q, t);
        }

        if (fmpz_sgn(fmpq_denref(res)) < 0)
        {
            fmpz_neg(fmpq_numref(res), fmpq_numref(res));
            fmpz_neg(fmpq_denref(res), fmpq_denref(res));
        }

        fmpz_clear(t);
    }
}
