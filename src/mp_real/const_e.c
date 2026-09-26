/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include "flint.h"
#include "mp_real.h"
#include "impl.h"

/* e = sum_{k >= 0} 1/k! = 1 + sum_{k >= 1} (1/k) prod_{j < k} (1/j):
   P = 1, Q = k, R = 1 with the k = 0 term as coefQ. */
void
_mp_real_const_e_compute(mp_real_t res, slong n)
{
    static const int64_t P[] = { INT64_C(1) };
    static const int64_t Q[] = { INT64_C(0), INT64_C(1) };
    static const int64_t R[] = { INT64_C(1) };

    mp_real_hypgeom_series_int64(res, 1, 1, 1, 1, P, 1, Q, 2, R, 1, n);
}
