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

/* Catalan's constant by Pilehrood's short series (2010),

       G = (1/64) sum_{k >= 0} 256^k (580k^2 - 184k + 15) (2k)!^3 (3k)!^2
                               / (k^3 (2k-1) (6k)!^2),

   in the hypergeometric form of y-cruncher's formula file
   (https://hal.inria.fr/hal-00990465/document).  Measured against the
   other files for Catalan at 10^5 and 10^6 bits, this is the fastest
   (Zuniga's 2023 series is 10% behind, Guillera's 2019 series 45%). */
void
_mp_real_const_catalan_compute(mp_real_t res, slong n)
{
    static const int64_t P[] = { INT64_C(15), INT64_C(-184), INT64_C(580) };
    static const int64_t Q[] = { INT64_C(225), INT64_C(-3240),
        INT64_C(14904), INT64_C(-23328), INT64_C(11664) };
    static const int64_t R[] = { INT64_C(0), INT64_C(0), INT64_C(0),
        INT64_C(-32), INT64_C(64) };

    mp_real_hypgeom_series_int64(res, 1, 1, 0, 2, P, 3, Q, 5, R, 5, n);
}
