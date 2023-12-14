/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

/*
    Find the fractions directly below and above a1/q1 in the Farey sequence of
    order Q:

     a0     a1     a2
    ---- < ---- < ----
     q0     q1     q2

    The index v satisfies

               Q + q0      q2 + q0     a2 + a0
    v = floor(--------) = --------- = ---------
                 q1           q1          a1
*/
void fmpq_farey_neighbors(fmpq_t left, fmpq_t right,
                                           const fmpq_t mid_, const fmpz_t Q_)
{
    fmpz_t Q, t;
    fmpq_t mid;

    /* find left denominator */
    if (fmpz_sgn(fmpq_denref(mid_)) <= 0
        || fmpz_cmp(fmpq_denref(mid_), Q_) > 0
        || !fmpz_invmod(fmpq_denref(left), fmpq_numref(mid_), fmpq_denref(mid_)))
        flint_throw(FLINT_ERROR, "(%s): bad input\n", __func__);

    /* simple handling of aliasing */
    fmpz_init_set(fmpq_numref(mid), fmpq_numref(mid_));
    fmpz_init_set(fmpq_denref(mid), fmpq_denref(mid_));
    fmpz_init_set(Q, Q_);
    fmpz_init(t);

    fmpz_sub(t, Q, fmpq_denref(left));
    fmpz_mod(t, t, fmpq_denref(mid));
    fmpz_sub(fmpq_denref(left), Q, t);

    /* find left numerator */
    fmpz_mul(t, fmpq_numref(mid), fmpq_denref(left));
    fmpz_sub_ui(t, t, 1);
    fmpz_divexact(fmpq_numref(left), t, fmpq_denref(mid));

    /* find index t */
    fmpz_add(t, Q, fmpq_denref(left));
    fmpz_fdiv_q(t, t, fmpq_denref(mid));

    /* find right denominator */
    fmpz_mul(fmpq_denref(mid), fmpq_denref(mid), t);
    fmpz_sub(fmpq_denref(right), fmpq_denref(mid), fmpq_denref(left));

    /* find right numerator */
    fmpz_mul(fmpq_numref(mid), fmpq_numref(mid), t);
    fmpz_sub(fmpq_numref(right), fmpq_numref(mid), fmpq_numref(left));

    fmpz_clear(fmpq_numref(mid));
    fmpz_clear(fmpq_denref(mid));
    fmpz_clear(Q);
    fmpz_clear(t);
}
