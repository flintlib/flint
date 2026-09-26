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

/* log 2 by Zuniga's d = 2 Ramanujan-type series [Zun2025, Eq. 18],
   11.92 bits per term,

       log 2 = (1/2160) sum_{k>=0} A(k) prod_{j=1}^{k} p(j)/q(j),
       A(k) = 1497 + 1794 k,  p(k) = k (2k - 1),
       q(k) = 7776 k^2 + 7776 k + 1080 = 216 (6k + 1)(6k + 5),

   through the generic binary splitting: c_Q = A(0) = 1497, c_P = 1,
   c_D = 2160, P = A p = 3588 k^3 + 1200 k^2 - 1497 k, Q = q, R = p.
   The series is computed at one limb more than the result, since
   mp_real_hypgeom_series gives FLINT_BITS (n - 1) accurate bits. */

static const int64_t log2_P[] = { 0, -1497, 1200, 3588 };
static const int64_t log2_Q[] = { 1080, 7776, 7776 };
static const int64_t log2_R[] = { 0, -1, 2 };

void
_mp_real_const_log2_compute(mp_real_t res, slong n)
{
    mp_real_hypgeom_series_int64(res, 1, 1, 1497, 2160,
        log2_P, 4, log2_Q, 3, log2_R, 3, n + 1);
}
