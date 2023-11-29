/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"
#include "acb_poly.h"

void
acb_dirichlet_hurwitz_precomp_eval(acb_t res,
    const acb_dirichlet_hurwitz_precomp_t pre, ulong p, ulong q, slong prec)
{
    slong i;
    acb_t a, t;

    if (p > q)
    {
        flint_throw(FLINT_ERROR, "hurwitz_precomp_eval: require p <= n\n");
    }

    if (pre->A == 0)
    {
        acb_init(a);
        acb_set_ui(a, p);
        acb_div_ui(a, a, q, prec);

        if (pre->deflate == 0)
            acb_hurwitz_zeta(res, &pre->s, a, prec);
        else
            _acb_poly_zeta_cpx_series(res, &pre->s, a, 1, 1, prec);

        acb_clear(a);
        return;
    }

    acb_init(a);
    acb_init(t);

    /* todo: avoid integer overflows below */
    if (p == q)
        i = pre->N - 1;
    else
        i = (pre->N * p) / q;
    acb_set_si(a, 2 * pre->N * p - q * (2 * i + 1));
    acb_div_ui(a, a, 2 * q * pre->N, prec);

    /* compute zeta(s,A+p/q) */
    _acb_poly_evaluate(res, pre->coeffs + i * pre->K, pre->K, a, prec);

    /* error bound */
    if (acb_is_real(&pre->s))
        arb_add_error_mag(acb_realref(res), &pre->err);
    else
        acb_add_error_mag(res, &pre->err);

    /* zeta(s,p/q) = (p/q)^-s + ... + (p/q+A-1)^-s zeta(s,A+p/q) */
    for (i = 0; i < pre->A; i++)
    {
        acb_set_ui(a, p);
        acb_div_ui(a, a, q, prec);
        acb_add_ui(a, a, i, prec);
        acb_neg(t, &pre->s);
        acb_pow(a, a, t, prec);
        acb_add(res, res, a, prec);
    }

    acb_clear(a);
    acb_clear(t);
}

