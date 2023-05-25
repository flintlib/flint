/*
    Copyright (C) 2022 Raoul Bourquin

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "qqbar.h"

void
qqbar_tribonacci_constant(qqbar_t res)
{
    /* Subexpressions */
    arb_t r33, r33p, r33m;

    fmpz_poly_zero(QQBAR_POLY(res));
    fmpz_poly_set_coeff_si(QQBAR_POLY(res), 3, 1);
    fmpz_poly_set_coeff_si(QQBAR_POLY(res), 2, -1);
    fmpz_poly_set_coeff_si(QQBAR_POLY(res), 1, -1);
    fmpz_poly_set_coeff_si(QQBAR_POLY(res), 0, -1);

    /* Init */
    arb_init(r33);
    arb_init(r33p);
    arb_init(r33m);

    /* r33 := 3*sqrt(33) */
    arb_sqrt_ui(r33, 33, QQBAR_DEFAULT_PREC);
    arb_mul_ui(r33, r33, 3, QQBAR_DEFAULT_PREC);

    /* r33p := cbrt(19 + r33) */
    arb_add_ui(r33p, r33, 19, QQBAR_DEFAULT_PREC);
    arb_root_ui(r33p, r33p, 3, QQBAR_DEFAULT_PREC);

    /* r33m := cbrt(19 - r33) */
    arb_sub_si(r33m, r33, 19, QQBAR_DEFAULT_PREC);
    arb_neg(r33m, r33m);
    arb_root_ui(r33m, r33m, 3, QQBAR_DEFAULT_PREC);

    /* res := 1 */
    arb_one(acb_realref(QQBAR_ENCLOSURE(res)));

    /* res += r33p */
    arb_add(acb_realref(QQBAR_ENCLOSURE(res)), acb_realref(QQBAR_ENCLOSURE(res)), r33p, QQBAR_DEFAULT_PREC);

    /* res += r33m */
    arb_add(acb_realref(QQBAR_ENCLOSURE(res)), acb_realref(QQBAR_ENCLOSURE(res)), r33m, QQBAR_DEFAULT_PREC);

    /* res /= 3 */
    arb_div_ui(acb_realref(QQBAR_ENCLOSURE(res)), acb_realref(QQBAR_ENCLOSURE(res)), 3, QQBAR_DEFAULT_PREC);

    /* zero imag part */
    arb_zero(acb_imagref(QQBAR_ENCLOSURE(res)));

    /* Free */
    arb_clear(r33);
    arb_clear(r33p);
    arb_clear(r33m);
}
