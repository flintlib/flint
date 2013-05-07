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

void
_fmpq_bsplit_sum_pq(fmpz_t P, fmpz_t Q, fmpz_t T,
    const fmpq * pq, long n1, long n2)
{
    if (n2 - n1 <= 0)
    {
        fmpz_zero(P);
        fmpz_one(Q);
    }
    else if (n2 - n1 == 1)
    {
        fmpz_set(P, fmpq_numref(pq + n1));
        fmpz_set(Q, fmpq_denref(pq + n1));
        fmpz_set(T, P);
    }
    else
    {
        long m = (n1 + n2) / 2;

        fmpz_t P2, Q2, T2;

        fmpz_init(P2);
        fmpz_init(Q2);
        fmpz_init(T2);

        _fmpq_bsplit_sum_pq(P,  Q,  T,  pq, n1, m);
        _fmpq_bsplit_sum_pq(P2, Q2, T2, pq, m, n2);

        fmpz_mul(T, T, Q2);

        fmpz_mul(T2, T2, P);
        fmpz_add(T, T, T2);

        fmpz_mul(P, P, P2);
        fmpz_mul(Q, Q, Q2);

        fmpz_clear(P2);
        fmpz_clear(Q2);
        fmpz_clear(T2);
    }
}

void
fmpq_bsplit_sum_pq(fmpq_bsplit_t s, const fmpq * pq, long n1, long n2)
{
    _fmpq_bsplit_sum_pq(s->P, s->Q, s->T, pq, n1, n2);
}
