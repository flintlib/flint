/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "arb.h"
#include "hypgeom.h"

void
arb_const_apery_eval(arb_t s, slong prec)
{
    hypgeom_t series;
    arb_t t;
    arb_init(t);
    hypgeom_init(series);

    fmpz_poly_set_str(series->A, "12  2539221463380 55800350430619 543035311757517 3094818610007883 11495362203169095 29262452925092202 52160137207884216 65180430489299744 56019281176488240 31567339140195744 10506477648460032 1565994397644288");
    fmpz_poly_set_str(series->B, "1  1");
    fmpz_poly_set_str(series->P, "15  0 0 0 0 0 -30 691 -6781 37374 -127976 283232 -406224 364896 -186624 41472");
    fmpz_poly_set_str(series->Q, "15  -44008272000 -2334151436400 -53522442803340 -703273183134030 -5931859745397870 -34140867105175650 -139058868850409430 -409481300311614720 -880500176512163280 -1382139595517666400 -1565294958171053280 -1244539247650560000 -658690593528960000 -208277254886400000 -29753893555200000");

    prec += 4 + FLINT_CLOG2(prec);
    arb_hypgeom_infsum(s, t, series, prec, prec);
    arb_mul_ui(t, t, 1031443875, prec);
    arb_mul_2exp_si(t, t, 11);
    arb_div(s, s, t, prec);

    hypgeom_clear(series);
    arb_clear(t);
}

ARB_DEF_CACHED_CONSTANT(arb_const_apery, arb_const_apery_eval)

