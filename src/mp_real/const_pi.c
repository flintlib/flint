/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mp_real.h"
#include "impl.h"

/* pi by the Chudnovsky series, through the generic binary splitting
   (y-cruncher's formula file for pi):

       pi = S / sqrt(10005),
       S = ((c_Q + c_P sum_{k>=1} P(k)/Q(k) prod_{j<k} R(j)/Q(j)) / c_D)^(-1),

   c_P = 1, c_Q = 13591409, c_D = 4270934400,
   R(k) = (6k - 5)(2k - 1)(6k - 1), P(k) = (13591409 + 545140134 k) R(k),
   Q(k) = -10939058860032000 k^3 = -(640320^3 / 24) k^3.
   The series is computed at one limb more than the result, since
   mp_real_hypgeom_series gives FLINT_BITS (n - 1) accurate bits. */

static const int64_t pi_P[] = { -67957045, -2100495856, 23608573992,
    -57896553024, 39250089648 };
static const int64_t pi_Q[] = { 0, 0, 0, -10939058860032000 };
static const int64_t pi_R[] = { -5, 46, -108, 72 };

void
_mp_real_const_pi_compute(mp_real_t pi, slong n)
{
    mp_real_t t;

    mp_real_init(t);
    mp_real_hypgeom_series_int64(pi, -1, 1, 13591409, 4270934400,
        pi_P, 5, pi_Q, 4, pi_R, 4, n + 1);
    mp_real_rsqrt_ui(t, 10005, n + 1);
    mp_real_mul(pi, pi, t, n);
    mp_real_clear(t);
}

void
_mp_real_const_pi4_compute(mp_real_t pi, slong n)
{
    _mp_real_const_pi_compute(pi, n);
    mp_real_mul_2exp_si(pi, pi, -2);
}

/* 2/pi = 1/(2 (pi/4)), for the argument reduction of the trigonometric
   functions */
void
_mp_real_const_2_div_pi_compute(mp_real_t res, slong n)
{
    mp_real_t p, one;

    mp_real_init(p);
    mp_real_init(one);
    _mp_real_const_pi4_compute(p, n + 1);
    mp_real_mul_2exp_si(p, p, 1);
    mp_real_set_ui(one, 1);
    mp_real_div(res, one, p, n);
    mp_real_clear(p);
    mp_real_clear(one);
}
