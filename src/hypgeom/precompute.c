/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "mag.h"
#include "hypgeom.h"

/* Compute a pure ratio P2(k)/Q2(k) for the term
   A(k)/B(k) * [P(1)P(2)...P(k)] / [Q(1)Q(2)...Q(k)] */
void
hypgeom_standardize(fmpz_poly_t P2, fmpz_poly_t Q2,
    const fmpz_poly_t A, const fmpz_poly_t B,
    const fmpz_poly_t P, const fmpz_poly_t Q)
{
    fmpz_t s;
    fmpz_poly_t T;

    fmpz_init(s);
    fmpz_poly_init(T);

    fmpz_set_si(s, -WORD(1));

    /* P = A * B(k-1) * P */
    fmpz_poly_taylor_shift(T, B, s);
    fmpz_poly_mul(P2, A, T);
    fmpz_poly_mul(P2, P2, P);

    /* Q = B * A(k-1) * Q */
    fmpz_poly_taylor_shift(T, A, s);
    fmpz_poly_mul(Q2, B, T);
    fmpz_poly_mul(Q2, Q2, Q);

    fmpz_clear(s);
    fmpz_poly_clear(T);
}


/* quotient of absolute value, rounded up */
static inline void
fmpz_cdiv_abs_q(fmpz_t q, const fmpz_t x, const fmpz_t y)
{
    if (fmpz_sgn(x) == fmpz_sgn(y))
    {
        fmpz_cdiv_q(q, x, y);
    }
    else
    {
        fmpz_fdiv_q(q, x, y);
        fmpz_neg(q, q);
    }
}

slong
hypgeom_root_norm(const fmpz_poly_t P)
{
    slong res, i, p;
    fmpz_t t, A;

    fmpz_init(A);
    fmpz_init(t);

    p = fmpz_poly_degree(P);

    fmpz_zero(A);

    for (i = 1; i <= p; i++)
    {
        fmpz_cdiv_abs_q(t, P->coeffs + p - i, P->coeffs + p);
        fmpz_root(t, t, i);
        fmpz_add_ui(t, t, 1);

        if (fmpz_cmp(t, A) > 0)
            fmpz_swap(t, A);
    }

    if (!fmpz_fits_si(A))
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    res = fmpz_get_si(A);

    fmpz_clear(A);
    fmpz_clear(t);

    return res;
}


static inline void
fmpz_poly_evaluate_si(fmpz_t y, const fmpz_poly_t poly, slong x)
{
    fmpz_set_si(y, x);
    fmpz_poly_evaluate_fmpz(y, poly, y);
}

void
_hypgeom_precompute(hypgeom_t hyp, const fmpz_poly_t P, const fmpz_poly_t Q)
{
    slong k;

    fmpz_t t;
    fmpz_init(t);

    hyp->r = fmpz_poly_degree(Q) - fmpz_poly_degree(P);
    hyp->boundC = hypgeom_root_norm(P);
    hyp->boundD = hypgeom_root_norm(Q);
    hyp->boundK = 1 + FLINT_MAX(hyp->boundC, 2 * hyp->boundD);

    mag_one(hyp->MK);

    for (k = 1; k <= hyp->boundK; k++)
    {
        fmpz_poly_evaluate_si(t, P, k);
        mag_mul_fmpz(hyp->MK, hyp->MK, t);

        fmpz_poly_evaluate_si(t, Q, k);
        mag_div_fmpz(hyp->MK, hyp->MK, t);
    }

    fmpz_clear(t);
}

void
hypgeom_precompute(hypgeom_t hyp)
{
    if (fmpz_poly_is_one(hyp->A) && fmpz_poly_is_one(hyp->B))
    {
        _hypgeom_precompute(hyp, hyp->P, hyp->Q);
    }
    else
    {
        fmpz_poly_t P2, Q2;
        fmpz_poly_init(P2);
        fmpz_poly_init(Q2);

        hypgeom_standardize(P2, Q2, hyp->A, hyp->B, hyp->P, hyp->Q);
        _hypgeom_precompute(hyp, P2, Q2);

        {
            fmpz_t t;
            fmpz_init(t);

            fmpz_poly_evaluate_si(t, hyp->A, 0);
            mag_mul_fmpz(hyp->MK, hyp->MK, t);

            fmpz_poly_evaluate_si(t, hyp->B, 0);
            mag_div_fmpz(hyp->MK, hyp->MK, t);

            fmpz_clear(t);
        }

        fmpz_poly_clear(P2);
        fmpz_poly_clear(Q2);
    }
}
