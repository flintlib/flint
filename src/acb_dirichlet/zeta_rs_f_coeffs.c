/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "acb_poly.h"
#include "acb_dirichlet.h"

void
acb_dirichlet_zeta_rs_f_coeffs(acb_ptr c, const arb_t p, slong N, slong prec)
{
    arb_ptr R, I, T, X;
    slong i, len;

    R = _arb_vec_init(N);
    I = _arb_vec_init(N);
    T = _arb_vec_init(N);
    X = _arb_vec_init(2);

    arb_set(X, p);
    arb_one(X + 1);

    /* I, R = sin,cos(pi*(X^2/2 + 3/8)) */
    len = FLINT_MIN(N, 3);
    _arb_poly_mullow(T, X, 2, X, 2, len, prec);
    _arb_vec_scalar_mul_2exp_si(T, T, len, -1);
    arb_set_d(R, 0.375);
    arb_add(T, T, R, prec);
    _arb_poly_sin_cos_pi_series(I, R, T, len, N, prec);

    /* I -= cos(pi*x/2) * sqrt(2) */
    _arb_vec_scalar_mul_2exp_si(X, X, 2, -1);
    _arb_poly_cos_pi_series(T, X, 2, N, prec);
    arb_sqrt_ui((arb_ptr) c, 2, prec);
    _arb_vec_scalar_mul(T, T, N, (arb_ptr) c, prec);
    _arb_vec_sub(I, I, T, N, prec);
    _arb_vec_scalar_mul_2exp_si(X, X, 2, 1);

    /* T = 1 / (2 cos(pi*x)) */
    _arb_poly_cos_pi_series(T, X, 2, N, prec);
    _arb_vec_scalar_mul_2exp_si(T, T, N, 1);
    _arb_poly_inv_series((arb_ptr) c, T, N, N, prec);
    _arb_vec_swap(T, (arb_ptr) c, N);

    /* R, I *= T */
    _arb_poly_mullow((arb_ptr) c, R, N, T, N, N, prec);
    _arb_vec_swap(R, (arb_ptr) c, N);
    _arb_poly_mullow((arb_ptr) c, I, N, T, N, N, prec);
    _arb_vec_swap(I, (arb_ptr) c, N);

    for (i = 0; i < N; i++)
    {
        arb_swap(acb_realref(c + i), R + i);
        arb_swap(acb_imagref(c + i), I + i);
    }

    _acb_poly_inv_borel_transform(c, c, N, prec);

    _arb_vec_clear(R, N);
    _arb_vec_clear(I, N);
    _arb_vec_clear(T, N);
    _arb_vec_clear(X, 2);
}

