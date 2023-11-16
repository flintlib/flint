/*
    Copyright (C) 2019 D.H.J Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

static void
_mag_div_ui_ui(mag_t res, ulong a, ulong b)
{
    mag_set_ui(res, a);
    mag_div_ui(res, res, b);
}

ulong
acb_dirichlet_turing_method_bound(const fmpz_t p)
{
    ulong result;
    arb_t t;
    fmpz_t k;
    mag_t m, b1, b2, c;

    arb_init(t);
    fmpz_init(k);
    mag_init(m);
    mag_init(b1);
    mag_init(b2);
    mag_init(c);

    acb_dirichlet_gram_point(t, p, NULL, NULL, FLINT_MAX(8, fmpz_bits(p)));
    arb_get_mag(m, t);
    mag_log(m, m);

    /* b1 = 0.0061*log(gram(p))^2 + 0.08*log(gram(p)) */
    _mag_div_ui_ui(c, 61, 10000);
    mag_mul(b1, c, m);
    _mag_div_ui_ui(c, 8, 100);
    mag_add(b1, b1, c);
    mag_mul(b1, b1, m);

    /* b2 = 0.0031*log(gram(p))^2 + 0.11*log(gram(p)) */
    _mag_div_ui_ui(c, 31, 10000);
    mag_mul(b2, c, m);
    _mag_div_ui_ui(c, 11, 100);
    mag_add(b2, b2, c);
    mag_mul(b2, b2, m);

    /* result = ceil(min(b1, b2)) */
    mag_min(m, b1, b2);
    mag_get_fmpz(k, m);
    result = fmpz_get_ui(k);

    arb_clear(t);
    fmpz_clear(k);
    mag_clear(m);
    mag_clear(b1);
    mag_clear(b2);
    mag_clear(c);

    return result;
}
