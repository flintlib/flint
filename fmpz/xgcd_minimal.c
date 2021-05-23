/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

void
fmpz_xgcd_minimal(fmpz_t d, fmpz_t a, fmpz_t b, const fmpz_t f, const fmpz_t g)
{
    fmpz_t kc; /* integer transforming new coefficients via ceil */
    fmpz_t kf; /* integer transforming new coefficients via floor */
    fmpz_t ca; /* ceil transformed `a` */
    fmpz_t cb; /* ceil transformed `b` */
    fmpz_t fa; /* floor transformed `a` */
    fmpz_t fb; /* floor transformed `b` */
    fmpz_t tmp;

    fmpz_xgcd(d, a, b, f, g);

    if (fmpz_is_zero(f) || fmpz_is_zero(g))
        return;

    fmpz_init(kc);
    fmpz_init(kf);
    fmpz_init_set(ca, a);
    fmpz_init_set(cb, b);
    fmpz_init_set(fa, a);
    fmpz_init_set(fb, b);
    fmpz_init(tmp);

    fmpz_mul(tmp, a, d);
    fmpz_cdiv_q(kc, tmp, g);
    fmpz_fdiv_q(kf, tmp, g);
    fmpz_divexact(tmp, g, d);
    fmpz_submul(ca, kc, tmp); /* ca = a - kc g / d */
    fmpz_submul(fa, kf, tmp); /* fa = a - kf g / d */
    fmpz_divexact(tmp, f, d);
    fmpz_addmul(cb, kc, tmp); /* cb = b + kc f / d */
    fmpz_addmul(fb, kf, tmp); /* fb = b + kf f / d */

    fmpz_abs(kc, ca); /* using kf and kc since they won't be of use anymore */
    fmpz_abs(tmp, cb);
    fmpz_add(kc, kc, tmp); /* kc = |ca| + |cb| */
    fmpz_abs(kf, fa);
    fmpz_abs(tmp, fb);
    fmpz_add(kf, kf, tmp); /* kf = |fa| + |fb| */

    if (fmpz_cmp(kf, kc) >= 0)
    {
        /* we want to use ceil */
        fmpz_abs(kf, a);
        fmpz_abs(tmp, b);
        fmpz_add(kf, kf, tmp); /* kf = |a| + |b| */
        if (fmpz_cmp(kf, kc) > 0)
        {
            fmpz_set(a, ca);
            fmpz_set(b, cb);
        }
    }
    else
    {
        /* we want to use floor */
        fmpz_abs(kc, a);
        fmpz_abs(tmp, b);
        fmpz_add(kc, kc, tmp); /* kc = |a| + |b| */
        if (fmpz_cmp(kc, kf) > 0)
        {
            fmpz_set(a, fa);
            fmpz_set(b, fb);
        }
    }

    fmpz_clear(kc);
    fmpz_clear(kf);
    fmpz_clear(ca);
    fmpz_clear(cb);
    fmpz_clear(fa);
    fmpz_clear(fb);
    fmpz_clear(tmp);
}
