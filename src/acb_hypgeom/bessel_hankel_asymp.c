/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Asymptotic expansions of the Hankel functions H1, H2, following the
   structure of acb_hypgeom_bessel_j_asymp: with

       omega = z - (2 nu + 1) pi / 4,
       U1 = U*(nu + 1/2, 2 nu + 1, -2iz),   U2 = U*(..., +2iz),

   the J asymptotic expansion is J = (2 pi z)^{-1/2} (e^{i omega} U1 +
   e^{-i omega} U2), and for Re(z) > 0 the Hankel functions are the
   individual terms:

       H1 = 2 (2 pi z)^{-1/2} e^{+i omega} U1,
       H2 = 2 (2 pi z)^{-1/2} e^{-i omega} U2. */

#include "acb.h"
#include "acb_hypgeom.h"

/* e^{+-i omega} and (2 pi z)^{-1/2} */
static void
hankel_asymp_prefactors(acb_t Ap, acb_t Am, acb_t C,
    const acb_t nu, const acb_t z, slong prec)
{
    /* the caller has already checked Re(z) > 0 */
    acb_t t, u;

    acb_init(t);
    acb_init(u);

    /* -(2nu+1)/4 * pi + z */
    acb_mul_2exp_si(t, nu, 1);
    acb_add_ui(t, t, 1, prec);
    acb_mul_2exp_si(t, t, -2);
    acb_neg(t, t);
    acb_const_pi(u, prec);
    acb_mul(t, t, u, prec);
    acb_add(t, t, z, prec);
    acb_mul_onei(t, t);
    acb_exp_invexp(Ap, Am, t, prec);

    /* (2 pi z)^(-1/2) */
    acb_const_pi(C, prec);
    acb_mul_2exp_si(C, C, 1);
    acb_mul(C, C, z, prec);
    acb_rsqrt(C, C, prec);

    acb_clear(t);
    acb_clear(u);
}

static void
bessel_h_asymp(acb_t res, const acb_t nu, const acb_t z, int which, slong prec)
{
    acb_t A1, A2, C, U, s, t, u;

    if (acb_is_finite(nu) && !acb_is_finite(z) &&
        acb_is_real(z) && !acb_contains_zero(z))
    {
        acb_zero(res);
        return;
    }

    /* The identification of the individual terms of the J expansion
       with H1 and H2 requires the principal branches in the prefactors
       to align with the Hankel functions' own expansions, which is the
       case for Re(z) > 0.  Outside this sector the terms pick up
       phases (their sum still gives J); a phase-corrected extension is
       possible but not implemented, so be conservative. */
    if (!arb_is_positive(acb_realref(z)))
    {
        acb_indeterminate(res);
        return;
    }

    acb_init(A1); acb_init(A2); acb_init(C);
    acb_init(U); acb_init(s); acb_init(t); acb_init(u);

    hankel_asymp_prefactors(A1, A2, C, nu, z, prec);

    if (!acb_is_finite(A1) || !acb_is_finite(A2) || !acb_is_finite(C))
    {
        acb_indeterminate(res);
    }
    else
    {
        /* s = 1/2 + nu, t = 1 + 2 nu */
        acb_one(s);
        acb_mul_2exp_si(s, s, -1);
        acb_add(s, s, nu, prec);
        acb_mul_2exp_si(t, nu, 1);
        acb_add_ui(t, t, 1, prec);

        acb_mul_onei(u, z);
        acb_mul_2exp_si(u, u, 1);
        if (which == 1)
            acb_neg(u, u);          /* H1: U*(s, t, -2iz) */
        acb_hypgeom_u_asymp(U, s, t, u, -1, prec);

        acb_mul(res, (which == 1) ? A1 : A2, U, prec);
        acb_mul(res, res, C, prec);
        acb_mul_2exp_si(res, res, 1);
    }

    acb_clear(A1); acb_clear(A2); acb_clear(C);
    acb_clear(U); acb_clear(s); acb_clear(t); acb_clear(u);
}

void
acb_hypgeom_bessel_h1_asymp(acb_t res, const acb_t nu, const acb_t z, slong prec)
{
    bessel_h_asymp(res, nu, z, 1, prec);
}

void
acb_hypgeom_bessel_h2_asymp(acb_t res, const acb_t nu, const acb_t z, slong prec)
{
    bessel_h_asymp(res, nu, z, 2, prec);
}
