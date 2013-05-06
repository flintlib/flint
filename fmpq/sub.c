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
_fmpq_sub(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q,
            const fmpz_t r, const fmpz_t s)
{
    fmpz_t g, a, b, t, u;

    /* Same denominator */
    if (fmpz_equal(q, s))
    {
        fmpz_sub(rnum, p, r);

        /* Both are integers */
        if (fmpz_is_one(q))
        {
            fmpz_set(rden, q);
        }
        else
        {
            fmpz_init(g);
            fmpz_gcd(g, rnum, q);

            if (fmpz_is_one(g))
            {
                fmpz_set(rden, q);
            }
            else
            {
                fmpz_divexact(rnum, rnum, g);
                fmpz_divexact(rden, q, g);
            }
            fmpz_clear(g);
        }
        return;
    }

    /* p/q is an integer */
    if (fmpz_is_one(q))
    {
        fmpz_init(t);
        fmpz_mul(t, p, s);
        fmpz_sub(rnum, t, r);
        fmpz_set(rden, s);
        fmpz_clear(t);
        return;
    }

    /* r/s is an integer */
    if (fmpz_is_one(s))
    {
        fmpz_init(t);
        fmpz_mul(t, r, q);
        fmpz_sub(rnum, p, t);
        fmpz_set(rden, q);
        fmpz_clear(t);
        return;
    }

    fmpz_init(g);
    fmpz_gcd(g, q, s);

    if (fmpz_is_one(g))
    {
        fmpz_init(t);
        fmpz_init(u);

        fmpz_mul(t, p, s);
        fmpz_mul(u, q, r);
        fmpz_sub(rnum, t, u);
        fmpz_mul(rden, q, s);

        fmpz_clear(t);
        fmpz_clear(u);
    }
    else
    {
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(t);
        fmpz_init(u);

        fmpz_divexact(a, q, g);
        fmpz_divexact(b, s, g);

        fmpz_mul(t, p, b);
        fmpz_mul(u, r, a);
        fmpz_sub(rnum, t, u);

        fmpz_gcd(t, rnum, g);

        if (fmpz_is_one(t))
        {
            fmpz_mul(rden, q, b);
        }
        else
        {
            fmpz_divexact(rnum, rnum, t);
            fmpz_divexact(g, q, t);
            fmpz_mul(rden, g, b);
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(t);
        fmpz_clear(u);
    }

    fmpz_clear(g);
}

void fmpq_sub(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
{
    _fmpq_sub(fmpq_numref(res), fmpq_denref(res),
              fmpq_numref(op1), fmpq_denref(op1),
              fmpq_numref(op2), fmpq_denref(op2));
}
