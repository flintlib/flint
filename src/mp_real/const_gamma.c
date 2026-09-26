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

/* Gamma(1/3) by Guillera's series (2023),

       Gamma(1/3) = (810^(1/4) pi / X)^(1/3),
       X = 3200 sum_{n>=1} (-9/64000)^n n^2 (9108n - 8829)
           (1/12)_n (5/12)_n / ((1)_n (1)_n (12n-7)(12n-11)),

   X in the hypergeometric form of y-cruncher's formula file
   (https://math.stackexchange.com/questions/4793437).  Measured at
   10^6 bits this is a shade faster than Zuniga's 2024 series and
   Brown's 2011 series (both within 10%), and the series is 4 times
   faster than the one the old arb code used. */
/* the series and pi are independent */
typedef struct
{
    mp_real_struct * x, * p;
    slong wp;
    int quarter;
}
gamma_work;

static void
_gamma_worker(slong i, void * arg)
{
    gamma_work * w = (gamma_work *) arg;
    static const int64_t P3[] = { INT64_C(0), INT64_C(0), INT64_C(-8829),
        INT64_C(9108) };
    static const int64_t Q3[] = { INT64_C(0), INT64_C(0), INT64_C(1024000) };
    static const int64_t R3[] = { INT64_C(-77), INT64_C(216), INT64_C(-144) };
    static const int64_t P4[] = { INT64_C(0), INT64_C(0),
        INT64_C(110094756835840), INT64_C(-440362836049920),
        INT64_C(440352218808320) };
    static const int64_t Q4[] = { INT64_C(0), INT64_C(0),
        INT64_C(11008380780544), INT64_C(-44033523122176),
        INT64_C(44033523122176) };
    static const int64_t R4[] = { INT64_C(3465), INT64_C(-35136),
        INT64_C(114176), INT64_C(-147456), INT64_C(65536) };

    if (i == 1)
        _mp_real_const_pi_compute(w->p, w->wp);
    else if (w->quarter == 2)
    {
        /* agm(1, sqrt 2) */
        mp_real_t t;
        mp_real_init(t);
        mp_real_set_ui(t, 2);
        mp_real_sqrt(t, t, w->wp);
        mp_real_set_ui(w->x, 1);
        mp_real_agm(w->x, w->x, t, w->wp);
        mp_real_clear(t);
    }
    else if (w->quarter)
        mp_real_hypgeom_series_int64(w->x, -1, 1, 0, 1, P4, 5, Q4, 5, R4, 5,
            w->wp);
    else
        mp_real_hypgeom_series_int64(w->x, 1, 3200, 0, 1, P3, 4, Q3, 3, R3, 3,
            w->wp);
}

void
_mp_real_const_gamma_1_3_compute(mp_real_t res, slong n)
{
    slong wp = n + 2;
    mp_real_t x, t, p;
    gamma_work w;

    mp_real_init(x);
    mp_real_init(t);
    mp_real_init(p);

    w.x = x;
    w.p = p;
    w.wp = wp;
    w.quarter = 0;
    _mp_real_parallel_tasks(_gamma_worker, &w, 2);

    /* 810^(1/4) pi */
    mp_real_set_ui(t, 810);
    mp_real_sqrt(t, t, wp);
    mp_real_sqrt(t, t, wp);
    mp_real_mul(t, t, p, wp);

    mp_real_div(t, t, x, wp);
    mp_real_root_ui(res, t, 3, n);

    mp_real_clear(x);
    mp_real_clear(t);
    mp_real_clear(p);
}

/* Gamma(1/4) from the lemniscate constant by Ebisu's series (2016),

       Gamma(1/4) = (pi^6 / (322 S^4))^(1/8),
       S = 1 / sum_{n>=1} n^2 (440352218808320n^2 - 440362836049920n
           + 110094756835840) / (11008380780544 n^2 (2n-1)^2)
           prod ... ,

   S in the hypergeometric form of y-cruncher's formula file (the
   eighth root is three square roots).  Ebisu's and Zuniga's 2023-x
   series measured within 1% of each other and both about 3 times
   faster than the AGM the old arb code used. */
/* from this many limbs Gamma(1/4) = sqrt((2 pi)^(3/2) / agm(1, sqrt 2)),
   the AGM (O(M(n) log n)) in parallel with pi, rather than the series
   (O(M(n) log^2 n)): measured equal from 10^5 to 4 10^6 bits on one
   thread and 4% faster at 10^7 bits, where on two threads it is the
   1.8 s of the AGM against the 2.2 s of the series alongside pi */
#ifndef MP_REAL_GAMMA_1_4_AGM_CUTOFF
#define MP_REAL_GAMMA_1_4_AGM_CUTOFF 65536
#endif

void
_mp_real_const_gamma_1_4_compute(mp_real_t res, slong n)
{
    slong wp = n + 2;
    mp_real_t s, t, p;
    gamma_work w;

    mp_real_init(s);
    mp_real_init(t);
    mp_real_init(p);

    w.x = s;
    w.p = p;
    w.wp = wp;
    w.quarter = (n >= MP_REAL_GAMMA_1_4_AGM_CUTOFF) ? 2 : 1;
    _mp_real_parallel_tasks(_gamma_worker, &w, 2);

    if (w.quarter == 2)
    {
        /* (2 pi)^(3/2) / agm(1, sqrt 2) */
        mp_real_mul_2exp_si(p, p, 1);
        mp_real_sqrt(t, p, wp);
        mp_real_mul(t, t, p, wp);
        mp_real_div(t, t, s, wp);
        mp_real_sqrt(res, t, n);
        mp_real_clear(s);
        mp_real_clear(t);
        mp_real_clear(p);
        return;
    }

    /* S^4 (the series has power -1) */
    mp_real_mul(s, s, s, wp);
    mp_real_mul(s, s, s, wp);
    mp_real_mul_ui(s, s, 322, wp);

    /* pi^6 */
    mp_real_mul(res, p, p, wp);
    mp_real_mul(t, res, p, wp);
    mp_real_mul(t, t, t, wp);

    mp_real_div(t, t, s, wp);
    mp_real_sqrt(t, t, wp);
    mp_real_sqrt(t, t, wp);
    mp_real_sqrt(res, t, n);

    mp_real_clear(s);
    mp_real_clear(t);
    mp_real_clear(p);
}
