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
#include "ulong_extras.h"
#include "arith.h"

#define p (fmpq_numref(s))
#define q (fmpq_denref(s))

void
dedekind_sum_coprime_large(fmpq_t s, const fmpz_t h, const fmpz_t k)
{
    fmpz_t a, b, t, u;
    int sign;

    if (fmpz_cmp_ui(k, 2UL) <= 0)
    {
        fmpq_zero(s);
        return;
    }

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(t);
    fmpz_init(u);

    fmpz_set(a, k);
    fmpz_set(b, h);

    sign = 1;

    fmpz_zero(p);
    fmpz_set_ui(q, 1UL);

    while (!fmpz_is_zero(b))
    {
        fmpz_mul(t, a, a);
        fmpz_addmul(t, b, b);
        fmpz_add_ui(t, t, 1UL);
        fmpz_mul(u, a, b);

        if (sign == -1)
            fmpz_neg(t, t);

        _fmpq_add(p, q, p, q, t, u);

        fmpz_mod(t, a, b);
        fmpz_swap(t, a);
        fmpz_swap(a, b);

        sign = -sign;
    }

    fmpz_mul_ui(q, q, 12UL);

    if (sign < 0)
    {
        fmpz_mul_2exp(p, p, 2UL);
        fmpz_sub(p, p, q);
        fmpz_mul_2exp(q, q, 2UL);
    }

    _fmpq_canonicalise(p, q);

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(t);
    fmpz_clear(u);
}
