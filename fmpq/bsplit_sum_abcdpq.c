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
fmpq_bsplit_sum_abcdpq(fmpq_bsplit_t s,
        const fmpq * ab, const fmpq * cd, const fmpq * pq, len_t n1, len_t n2)
{
    if (n2 <= n1)
    {
        return;
    }
    if (n1 == n2 - 1)
    {
        fmpz_set(s->P, fmpq_numref(pq + n1));
        fmpz_set(s->Q, fmpq_denref(pq + n1));
        fmpz_set(s->B, fmpq_denref(ab + n1));
        fmpz_set(s->C, fmpq_numref(cd + n1));
        fmpz_set(s->D, fmpq_denref(cd + n1));
        fmpz_mul(s->T, fmpq_numref(ab + n1), fmpq_numref(pq + n1));
        fmpz_mul(s->V, fmpq_numref(ab + n1), fmpq_numref(cd + n1));
        fmpz_mul(s->V, s->V, fmpq_numref(pq + n1));
    }
    else
    {
        fmpq_bsplit_t L, R;
        fmpz_t t, u, v;

        len_t m = (n1 + n2) / 2;

        fmpq_bsplit_init(L);
        fmpq_bsplit_init(R);

        fmpq_bsplit_sum_abcdpq(L, ab, cd, pq, n1, m);
        fmpq_bsplit_sum_abcdpq(R, ab, cd, pq, m, n2);

        fmpz_init(t);
        fmpz_init(u);
        fmpz_init(v);

        fmpz_mul(s->P, L->P, R->P);
        fmpz_mul(s->Q, L->Q, R->Q);
        fmpz_mul(s->B, L->B, R->B);
        fmpz_mul(s->D, L->D, R->D);

        /* T = LB LP RT + RB RQ LT*/
        fmpz_mul(u, L->B, L->P);
        fmpz_mul(t, u, R->T);
        fmpz_mul(v, R->B, R->Q);
        fmpz_mul(s->T, v, L->T);
        fmpz_add(s->T, s->T, t);

        /* C = LC RD RC LD */
        fmpz_mul(s->C, L->C, R->D);
        fmpz_addmul(s->C, R->C, L->D);

        /* V = RD (RB RQ LV + LC LB LP RT) + LD LB LP RV */
        fmpz_mul(u, u, R->V);
        fmpz_mul(u, u, L->D);
        fmpz_mul(v, v, L->V);
        fmpz_addmul(v, t, L->C);
        fmpz_mul(v, v, R->D);
        fmpz_add(s->V, u, v);

        fmpz_clear(t);
        fmpz_clear(u);
        fmpz_clear(v);

        fmpq_bsplit_clear(L);
        fmpq_bsplit_clear(R);
    }
}
