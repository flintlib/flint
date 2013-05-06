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
arith_dedekind_sum_naive(fmpq_t s, const fmpz_t h, const fmpz_t k)
{
    fmpz_t i, j, q1, r1, q2, r2;

    if (fmpz_is_zero(k))
    {
        fmpq_zero(s);
        return;
    }

    fmpz_init(i);
    fmpz_init(j);
    fmpz_init(q1);
    fmpz_init(r1);
    fmpz_init(q2);
    fmpz_init(r2);

    fmpz_zero(fmpq_numref(s));

    for (fmpz_one(i); fmpz_cmp(i, k) < 0; fmpz_add_ui(i, i, 1))
    {
        fmpz_fdiv_qr(q1, r1, i, k);
        if (fmpz_is_zero(r1))
            continue;

        fmpz_mul(j, h, i);
        fmpz_fdiv_qr(q2, r2, j, k);
        if (fmpz_is_zero(r2))
            continue;

        fmpz_mul(q1, q1, k);
        fmpz_sub(q1, i, q1);
        fmpz_mul_ui(q1, q1, 2);
        fmpz_sub(q1, q1, k);

        fmpz_mul(q2, q2, k);
        fmpz_sub(q2, j, q2);
        fmpz_mul_ui(q2, q2, 2);
        fmpz_sub(q2, q2, k);

        fmpz_addmul(fmpq_numref(s), q1, q2);
    }

    fmpz_mul(fmpq_denref(s), k, k);
    fmpz_mul_ui(fmpq_denref(s), fmpq_denref(s), 4);
    fmpq_canonicalise(s);

    fmpz_clear(i);
    fmpz_clear(j);
    fmpz_clear(q1);
    fmpz_clear(r1);
    fmpz_clear(q2);
    fmpz_clear(r2);
}
