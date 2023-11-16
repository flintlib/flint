/*
    Copyright (C) 2019 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

/*
Claim: the error is bounded by 1/64 if n <= 1 and (1/64) (log(n)/n) if n >= 2.

A crude lower bound for g_n is 2 pi exp(W(n)), or 8*n/log(n) for n >= 8.

We want to solve pi n = -t/2 log(2 pi/t) - t/2 - pi/8 + epsilon for t (= g_n).
Using (47) in Brent [https://arxiv.org/abs/1609.03682], |epsilon| <= 1/(8 t) for for t >= 2.
Also, for x >= 3, |f'(x)| < 0.5 where f(x) = exp(W(x)).

Assume n >= 9, so that (n+1/8)/e >= 3.35. Then inverting gives

t = 2 pi exp[W( [pi n - epsilon + pi/8] / (pi e) ) + 1]
  = 2 pi e exp[W((n+1/8)/e  - epsilon / (pi e))]
  = 2 pi e exp[W((n+1/8)/e)] + epsilon2, |epsilon2| <= 1/(8 t) <= (1/64) (log(n)/n)

One can check 0 <= n <= 8 separately.
*/
static void
gram_point_initial(arb_t x, const fmpz_t n, slong prec)
{
    arb_t pi, e;
    mag_t b;

    arb_init(pi);
    arb_init(e);
    mag_init(b);

    arb_const_pi(pi, prec);
    arb_const_e(e, prec);

    /* x = 2*pi*exp(1 + W((n+1/8)/e)) */
    arb_one(x);
    arb_mul_2exp_si(x, x, -3);
    arb_add_fmpz(x, x, n, prec);
    arb_div(x, x, e, prec);
    arb_lambertw(x, x, 0, prec);
    arb_add_ui(x, x, 1, prec);
    arb_exp(x, x, prec);
    arb_mul(x, x, pi, prec);
    arb_mul_2exp_si(x, x, 1);

    if (fmpz_cmp_ui(n, 1) <= 0)
    {
        mag_set_ui_2exp_si(b, 1, -6);
    }
    else
    {
        mag_set_fmpz(b, n);
        mag_log(b, b);
        mag_div_fmpz(b, b, n);
        mag_mul_2exp_si(b, b, -6);
    }

    arb_add_error_mag(x, b);

    arb_clear(pi);
    arb_clear(e);
    mag_clear(b);
}

void
acb_dirichlet_gram_point(arb_t res, const fmpz_t n, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
{
    slong asymp_accuracy;

    /* Only implemented for n >= -1 and Riemann zeta. */
    if (fmpz_cmp_si(n, -1) < 0 || G != NULL || chi != NULL)
    {
        arb_indeterminate(res);
        return;
    }

    asymp_accuracy = 2 * fmpz_bits(n);
    asymp_accuracy = FLINT_MIN(asymp_accuracy, prec);

    gram_point_initial(res, n, asymp_accuracy + 20);
    asymp_accuracy = arb_rel_accuracy_bits(res);

    if (asymp_accuracy < prec)
    {
        acb_struct tmp[2];
        arb_t f, fprime, root;
        mag_t C, r;
        slong * steps;
        slong wp, step;

        acb_init(tmp);
        acb_init(tmp + 1);
        arb_init(f);
        arb_init(fprime);
        arb_init(root);
        mag_init(C);
        mag_init(r);
        steps = flint_malloc(sizeof(slong) * FLINT_BITS);

        step = 0;
        steps[step] = prec * 1.05 + 10;

        while (steps[step] / 2 > asymp_accuracy)
        {
            steps[step + 1] = steps[step] / 2;
            step++;
        }

        arb_set(root, res);

        /* theta''(x) <= C = 1/x, x >= 1 */
        arb_get_mag_lower(C, root);
        if (mag_cmp_2exp_si(C, 0) >= 0)
            mag_inv(C, C);
        else
            mag_inf(C);

        arb_set(root, res);

        for ( ; step >= 0; step--)
        {
            wp = steps[step] + 10;
            wp = FLINT_MAX(wp, arb_rel_accuracy_bits(root) + 10);

            /* store radius, set root to the midpoint */
            mag_set(r, arb_radref(root));
            mag_zero(arb_radref(root));

            acb_set_arb(tmp, root);
            acb_dirichlet_hardy_theta(tmp, tmp, NULL, NULL, 2, wp);
            arb_set(f, acb_realref(tmp));
            arb_const_pi(acb_imagref(tmp), wp);
            arb_submul_fmpz(f, acb_imagref(tmp), n, wp);

            arb_set(fprime, acb_realref(tmp + 1));

            /* f'([m+/-r]) = f'(m) +/- f''([m +/- r]) * r */
            mag_mul(r, C, r);
            arb_add_error_mag(fprime, r);
            arb_div(f, f, fprime, wp);
            arb_sub(root, root, f, wp);

            /* Verify inclusion so that C is still valid. */
            if (!arb_contains(res, root))
            {
                flint_printf("unexpected: no containment computing Gram point\n");
                arb_set(root, res);
                break;
            }
        }

        arb_set(res, root);

        acb_clear(tmp);
        acb_clear(tmp + 1);
        arb_clear(f);
        arb_clear(fprime);
        arb_clear(root);
        mag_clear(C);
        mag_clear(r);
        flint_free(steps);
    }

    arb_set_round(res, res, prec);
}
