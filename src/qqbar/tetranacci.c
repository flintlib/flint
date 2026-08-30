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
qqbar_tetranacci_constant(qqbar_t res)
{
    /* Subexpressions */
    arb_t t1, t2, l1, p1;

    fmpz_poly_zero(QQBAR_POLY(res));
    fmpz_poly_set_coeff_si(QQBAR_POLY(res), 4, 1);
    fmpz_poly_set_coeff_si(QQBAR_POLY(res), 3, -1);
    fmpz_poly_set_coeff_si(QQBAR_POLY(res), 2, -1);
    fmpz_poly_set_coeff_si(QQBAR_POLY(res), 1, -1);
    fmpz_poly_set_coeff_si(QQBAR_POLY(res), 0, -1);

    /* Init */
    arb_init(t1);
    arb_init(t2);
    arb_init(l1);
    arb_init(p1);

    /* t1 := 3*sqrt(1689) */
    arb_sqrt_ui(t1, 1689, QQBAR_DEFAULT_PREC);
    arb_mul_ui(t1, t1, 3, QQBAR_DEFAULT_PREC);

    /* l1 := (cbrt(t1 - 65) - cbrt(t1 + 65)) / (12 * cbrt(2)) */
    arb_zero(l1);

    /* First part of l1 */
    arb_sub_ui(t2, t1, 65, QQBAR_DEFAULT_PREC);
    arb_root_ui(t2, t2, 3, QQBAR_DEFAULT_PREC);

    arb_add(l1, l1, t2, QQBAR_DEFAULT_PREC);

    /* Second part of l1 */
    arb_add_ui(t2, t1, 65, QQBAR_DEFAULT_PREC);
    arb_root_ui(t2, t2, 3, QQBAR_DEFAULT_PREC);

    arb_sub(l1, l1, t2, QQBAR_DEFAULT_PREC);

    /* Denominator of l1 */
    arb_set_ui(t2, 2);
    arb_root_ui(t2, t2, 3, QQBAR_DEFAULT_PREC);
    arb_mul_ui(t2, t2, 12, QQBAR_DEFAULT_PREC);

    /* Combine l1 */
    arb_div(l1, l1, t2, QQBAR_DEFAULT_PREC);

    /* p1 := sqrt(l1 + 11/48) */
    arb_set_ui(p1, 11);
    arb_div_ui(p1, p1, 48, QQBAR_DEFAULT_PREC);

    arb_add(p1, p1, l1, QQBAR_DEFAULT_PREC);
    arb_sqrt(p1, p1, QQBAR_DEFAULT_PREC);

    /* t2 := 1/4 - p1 */
    arb_one(t2);
    arb_mul_2exp_si(t2, t2, -2);
    arb_add(t2, t2, p1, QQBAR_DEFAULT_PREC);

    /* Invert p1 */
    arb_inv(p1, p1, QQBAR_DEFAULT_PREC);

    /* Final result */
    arb_one(acb_realref(QQBAR_ENCLOSURE(res)));
    arb_div_ui(acb_realref(QQBAR_ENCLOSURE(res)), acb_realref(QQBAR_ENCLOSURE(res)), 6, QQBAR_DEFAULT_PREC);

    arb_set_ui(t1, 7);
    arb_div_ui(t1, t1, 24, QQBAR_DEFAULT_PREC);
    arb_mul(t1, t1, p1, QQBAR_DEFAULT_PREC);
    arb_add(acb_realref(QQBAR_ENCLOSURE(res)), acb_realref(QQBAR_ENCLOSURE(res)), t1, QQBAR_DEFAULT_PREC);

    arb_sqr(t1, t2, QQBAR_DEFAULT_PREC);
    arb_add(acb_realref(QQBAR_ENCLOSURE(res)), acb_realref(QQBAR_ENCLOSURE(res)), t1, QQBAR_DEFAULT_PREC);

    arb_mul_ui(t1, t2, 2, QQBAR_DEFAULT_PREC);
    arb_mul(t1, t1, l1, QQBAR_DEFAULT_PREC);
    arb_mul(t1, t1, p1, QQBAR_DEFAULT_PREC);
    arb_sub(acb_realref(QQBAR_ENCLOSURE(res)), acb_realref(QQBAR_ENCLOSURE(res)), t1, QQBAR_DEFAULT_PREC);

    arb_sqrt(acb_realref(QQBAR_ENCLOSURE(res)), acb_realref(QQBAR_ENCLOSURE(res)), QQBAR_DEFAULT_PREC);

    arb_add(acb_realref(QQBAR_ENCLOSURE(res)), acb_realref(QQBAR_ENCLOSURE(res)), t2, QQBAR_DEFAULT_PREC);

    /* zero imag part */
    arb_zero(acb_imagref(QQBAR_ENCLOSURE(res)));

    /* Free */
    arb_clear(t1);
    arb_clear(t2);
    arb_clear(l1);
    arb_clear(p1);
}
