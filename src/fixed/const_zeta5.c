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
#include "fixed.h"

/* zeta(5) by Zhi-Wei Sun's identity (2025),

       zeta(5) = (3 S + 56 pi^2 zeta(3)) / 540,
       S = sum_{n>=1} (-1/27)^n (560n^4 - 640n^3 + 408n^2 - 136n + 17)
           (1)_n (1)_n / ((1/3)_n (2/3)_n n^5 (2n-1)^4),

   S in the hypergeometric form of y-cruncher's formula file
   (https://mathoverflow.net/questions/486353).  The cost is 10.9 bits
   per term for S plus zeta(3) and pi, against 32.3 for Y. Zhao's
   single series (measured 3.6 times slower at 10^5 bits). */
/* the three independent parts: S, zeta(3) and pi */
typedef struct
{
    fball_struct * s, * z, * p;
    slong wp;
}
zeta5_work;

static void
_zeta5_worker(slong i, void * arg)
{
    zeta5_work * w = (zeta5_work *) arg;
    static const int64_t P[] = { INT64_C(17), INT64_C(-136), INT64_C(408),
        INT64_C(-640), INT64_C(560) };
    static const int64_t Q[] = { INT64_C(0), INT64_C(0), INT64_C(0),
        INT64_C(-6), INT64_C(75), INT64_C(-387), INT64_C(1056),
        INT64_C(-1608), INT64_C(1296), INT64_C(-432) };
    static const int64_t R[] = { INT64_C(0), INT64_C(0), INT64_C(0),
        INT64_C(0), INT64_C(0), INT64_C(1), INT64_C(-8), INT64_C(24),
        INT64_C(-32), INT64_C(16) };

    if (i == 0)
        fball_hypgeom_series_int64(w->s, 1, 3, 0, 1, P, 5, Q, 10, R, 10,
            w->wp);
    else if (i == 1)
        fball_const_zeta3(w->z, w->wp);
    else
        fball_const_pi_chudnovsky(w->p, w->wp);
}

void
fball_const_zeta5(fball_t res, slong n)
{
    slong wp = n + 2;
    fball_t s, z, p, t;
    zeta5_work w;

    fball_init(s);
    fball_init(z);
    fball_init(p);
    fball_init(t);

    w.s = s;
    w.z = z;
    w.p = p;
    w.wp = wp;
    _fixed_parallel_tasks(_zeta5_worker, &w, 3);

    /* (3 S + 56 pi^2 zeta(3)) / 540 */
    fball_mul(p, p, p, wp);
    fball_mul_ui(p, p, 56, wp);
    fball_mul(p, p, z, wp);
    fball_add(s, s, p, wp);
    fball_set_ui(t, 540);
    fball_div(res, s, t, n);

    fball_clear(s);
    fball_clear(z);
    fball_clear(p);
    fball_clear(t);
}
