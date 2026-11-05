/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_vec.h"
#include "fmpq_mat.h"
#include "acb.h"
#include "impl.h"

/* The NaCl Madelung constant M = -Z_{Z^3,1}(0; (1/2,1/2,1/2)), by the shell
   algorithm (acb_mzeta_epstein_gram_shell): the direct side is the
   alternating sum of r3(N) Gamma(1/2, lambda N) values, the reciprocal side
   the power series in q with coefficients r3(8k + 3)/(8k + 3), whose counts
   come from the t2 divisor sieve at high precision (the r3_odd handler of the
   shell counts). */
void
acb_mzeta_madelung_nacl(arb_t res, slong prec)
{
    fmpq_mat_t G;
    fmpq * u, * v;
    acb_t nu, z;
    slong i;

    fmpq_mat_init(G, 3, 3);
    fmpq_mat_one(G);
    u = _fmpq_vec_init(3);
    v = _fmpq_vec_init(3);
    for (i = 0; i < 3; i++)
        fmpq_set_si(v + i, 1, 2);
    acb_init(nu);
    acb_init(z);
    acb_one(nu);

    if (acb_mzeta_epstein_gram_shell(z, nu, G, u, v, 0, FLINT_MAX(prec, 2)))
        arb_neg(res, acb_realref(z));
    else
        arb_indeterminate(res);

    fmpq_mat_clear(G);
    _fmpq_vec_clear(u, 3);
    _fmpq_vec_clear(v, 3);
    acb_clear(nu);
    acb_clear(z);
}
