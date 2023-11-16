/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"

void
_arb_hypgeom_gamma_stirling_term_bounds(slong * bound, const mag_t zinv, slong N)
{
    mag_t b, u;
    slong n;

    mag_init(b);
    mag_init(u);

    /* bound[0] = WORD_MAX;  -- should not be used */

    /* first term 1/(12z) */
    mag_set(b, zinv);
    mag_div_ui(b, b, 12);
    bound[1] = MAG_EXP(b);

    /* u = 1/(2 pi z)^2 */
    mag_const_pi_lower(u);
    mag_mul_2exp_si(u, u, 1);
    mag_inv(u, u);
    mag_mul(u, u, zinv);
    mag_mul(u, u, u);

    /* zeta(2n) 2 (2n-2)! / (2pi)^(2n) / z^(2n-1) */
    /* ratio bounded by (2n-2)(2n-3)/(2 pi z)^2 */

    for (n = 2; n < N; n++)
    {
        mag_mul_ui(b, b, (2*n-2) * (2*n-3));
        mag_mul(b, b, u);
        bound[n] = MAG_EXP(b);
    }

    mag_clear(b);
    mag_clear(u);
}

