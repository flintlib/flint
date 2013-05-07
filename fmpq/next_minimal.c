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
_fmpq_next_minimal(fmpz_t rnum, fmpz_t rden,
                        const fmpz_t num, const fmpz_t den)
{
    fmpz p, q;

    p = *num;
    q = *den;

    if (!COEFF_IS_MPZ(p) && !COEFF_IS_MPZ(q))
    {
        if (p < q && p)
        {
            *rnum = q;
            *rden = p;
            return;
        }

        while (q < p)
        {
            q++;
            if (n_gcd(p, q) == 1)
            {
                *rnum = q;
                *rden = p;
                return;
            }
        }

        *rnum = 1;
        fmpz_set_ui(rden, p + 1);
    }
    else
    {
        fmpz_t t;

        if (fmpz_cmp(num, den) < 0)
        {
            fmpz_set(rnum, num);
            fmpz_set(rden, den);
            fmpz_swap(rnum, rden);
            return;
        }

        fmpz_init(t);
        fmpz_set(rnum, num);
        fmpz_set(rden, den);

        while (fmpz_cmp(rden, rnum) < 0)
        {
            fmpz_add_ui(rden, rden, 1UL);
            fmpz_gcd(t, rden, rnum);
            if (fmpz_is_one(t))
            {
                fmpz_swap(rnum, rden);
                fmpz_clear(t);
                return;
            }
        }

        fmpz_add_ui(rden, rden, 1UL);
        fmpz_one(rnum);
        fmpz_clear(t);
    }
}

void
fmpq_next_minimal(fmpq_t res, const fmpq_t x)
{
    _fmpq_next_minimal(&res->num, &res->den, &x->num, &x->den);
}
