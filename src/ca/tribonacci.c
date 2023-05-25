/*
    Copyright (C) 2022 Raoul Bourquin

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_tribonacci_constant(ca_t res, ca_ctx_t ctx)
{
    /* Subexpressions */
    ca_t third, r33, r33p, r33m;

    /* Init */
    ca_init(third, ctx);
    ca_init(r33, ctx);
    ca_init(r33p, ctx);
    ca_init(r33m, ctx);

    /* third := 1/3 */
    ca_one(third, ctx);
    ca_div_ui(third, third, 3, ctx);

    /* r33 := 3*sqrt(33) */
    ca_sqrt_ui(r33, 33, ctx);
    ca_mul_ui(r33, r33, 3, ctx);

    /* r33p := cbrt(19 + r33) */
    ca_add_ui(r33p, r33, 19, ctx);
    /* Todo: now suitable root function in calcium yet */
    /* ca_root_ui(r33p, r33p, 3, ctx); */
    ca_pow(r33p, r33p, third, ctx);

    /* r33m := cbrt(19 - r33) */
    ca_sub_si(r33m, r33, 19, ctx);
    ca_neg(r33m, r33m, ctx);
    /* Todo: now suitable root function in calcium yet */
    /* ca_root_ui(r33m, r33m, 3, ctx); */
    ca_pow(r33m, r33m, third, ctx);

    /* res := 1 */
    ca_one(res, ctx);
    /* res += r33p */
    ca_add(res, res, r33p, ctx);
    /* res += r33m */
    ca_add(res, res, r33m, ctx);
    /* res /= 3 */
    ca_div_ui(res, res, 3, ctx);

    /* Free */
    ca_clear(third, ctx);
    ca_clear(r33, ctx);
    ca_clear(r33p, ctx);
    ca_clear(r33m, ctx);
}
