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
#include "arith.h"

void
arith_dedekind_sum(fmpq_t s, const fmpz_t h, const fmpz_t k)
{
    if (fmpz_cmp_ui(k, 2UL) <= 0 || fmpz_is_zero(h) || fmpz_equal(h, k))
    {
        fmpq_zero(s);
    }
    else if (fmpz_sgn(h) < 0)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_neg(t, h);
        arith_dedekind_sum(s, t, k);
        fmpq_neg(s, s);
        fmpz_clear(t);
    }
    else
    {
        fmpz_t t, u, q;

        fmpz_init(t);
        fmpz_init(u);
        fmpz_init(q);

        fmpz_gcd(q, h, k);
        fmpz_divexact(t, h, q);
        fmpz_divexact(u, k, q);

        if (fmpz_cmp(t, u) > 0)
        {
            fmpq_t r;
            fmpq_init(r);

            /* r = (1 + h(h-3k) + k^2) / (12hk) */
            fmpz_mul_ui(fmpq_numref(r), u, 3UL);
            fmpz_sub(fmpq_numref(r), t, fmpq_numref(r));
            fmpz_mul(fmpq_numref(r), fmpq_numref(r), t);
            fmpz_addmul(fmpq_numref(r), u, u);
            fmpz_add_ui(fmpq_numref(r), fmpq_numref(r), 1UL);
            fmpz_mul(fmpq_denref(r), t, u);
            fmpz_mul_ui(fmpq_denref(r), fmpq_denref(r), 12UL);
            fmpq_canonicalise(r);
            arith_dedekind_sum_coprime(s, u, t);
            fmpq_sub(s, r, s);

            fmpq_clear(r);
        }
        else
        {
            arith_dedekind_sum_coprime(s, t, u);
        }

        fmpz_clear(t);
        fmpz_clear(u);
        fmpz_clear(q);
    }
}
