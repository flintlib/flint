/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_calc.h"
#include "acb_dirichlet.h"

static void
integral_tail(mag_t bound, const acb_t z, const acb_t log_z, const acb_t s, const acb_t a, const arb_t R, slong prec)
{
    arb_t C, s1;
    mag_t t;

    arb_init(C);
    arb_init(s1);
    mag_init(t);

    /* re(s) - 1 */
    arb_sub_ui(s1, acb_realref(s), 1, prec);

    /* C = re(a) - max(0, re(s) - 1) / R */
    arb_nonnegative_part(C, s1);
    arb_div(C, C, R, prec);
    arb_sub(C, acb_realref(a), C, prec);

    /* C > 0 */
    if (arb_is_positive(C))
    {
        /* |log(z)| + 1 */
        acb_get_mag(bound, log_z);
        mag_add_ui(bound, bound, 1);
        /* R */
        arb_get_mag_lower(t, R);

        /* R > |log(z)| + 1 */
        if (mag_cmp(t, bound) > 0)
        {
            /* 2 * R^(re(s)-1) * C^(-1) * exp(-re(a)*R)  */
            arb_pow(s1, R, s1, prec);
            arb_div(C, s1, C, prec);
            arb_mul_2exp_si(C, C, 1);
            arb_mul(s1, acb_realref(a), R, prec);
            arb_neg(s1, s1);
            arb_exp(s1, s1, prec);
            arb_mul(C, C, s1, prec);
            arb_get_mag(bound, C);
        }
        else
        {
            mag_inf(bound);
        }
    }
    else
    {
        mag_inf(bound);
    }

    arb_clear(C);
    arb_clear(s1);
    mag_clear(t);
}

static int
_integrand(acb_ptr res, const acb_t t, void * param, slong order, int negate_power, slong prec)
{
    acb_srcptr z, s, a;
    acb_t u, v;

    if (order > 1)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);  /* Would be needed for Taylor method. */

    z = ((acb_srcptr)(param)) + 0;
    s = ((acb_srcptr)(param)) + 1;
    a = ((acb_srcptr)(param)) + 2;

    acb_init(u);
    acb_init(v);

    acb_neg(u, t);
    acb_exp(u, u, prec);
    acb_mul(u, u, z, prec);
    acb_sub_ui(u, u, 1, prec);
    acb_neg(u, u);

    if (acb_contains_zero(u))
    {
        acb_indeterminate(res);
    }
    else
    {
        /* t^(s-1) * exp(-a*t) = exp((s-1)*log(t) - a*t) */
        /* (-t)^(s-1) * exp(-a*t) = exp((s-1)*log(t) - a*t) */

        acb_sub_ui(v, s, 1, prec);

        if (acb_is_int(s))
        {
            if (negate_power)
            {
                acb_neg(res, t);
                acb_pow(v, res, v, prec);
            }
            else
            {
                acb_pow(v, t, v, prec);
            }

            acb_div(u, v, u, prec);
            acb_mul(v, a, t, prec);
            acb_neg(v, v);
            acb_exp(v, v, prec);
            acb_mul(res, u, v, prec);
        }
        else
        {
            if (negate_power)
            {
                acb_neg(res, t);
                acb_log_analytic(res, res, order != 0, prec);
            }
            else
            {
                acb_log_analytic(res, t, order != 0, prec);
            }

            acb_mul(res, res, v, prec);
            acb_submul(res, a, t, prec);
            acb_exp(res, res, prec);
            acb_div(res, res, u, prec);
        }
    }

    acb_clear(u);
    acb_clear(v);

    return 0;
}

static int
integrand(acb_ptr res, const acb_t t, void * param, slong order, slong prec)
{
    return _integrand(res, t, param, order, 0, prec);
}

static int
integrand2(acb_ptr res, const acb_t t, void * param, slong order, slong prec)
{
    return _integrand(res, t, param, order, 1, prec);
}

void
_acb_dirichlet_lerch_phi_integral(acb_t res, const acb_t z, const acb_t s, const acb_t a, slong prec)
{
    acb_t log_z, t, u, v, w, zero, N;
    acb_ptr param;
    mag_t abs_tol, tail_bound;
    slong i, rel_goal;
    acb_calc_integrate_opt_t options;
    int is_real;
    mag_t log_z_im_lower;
    acb_t xa, xb;

    acb_init(log_z);
    acb_init(t);
    acb_init(u);
    acb_init(v);
    acb_init(w);
    acb_init(zero);
    acb_init(N);
    mag_init(abs_tol);
    mag_init(tail_bound);
    mag_init(log_z_im_lower);
    acb_init(xa);
    acb_init(xb);
    param = _acb_vec_init(3);

    /* compute either the principal log or +/ 2 pi i
       (does not matter as long as the imaginary part is not too far off) */
    if (arb_is_positive(acb_realref(z)) || !arb_contains_zero(acb_imagref(z)))
    {
        acb_log(log_z, z, prec);
    }
    else
    {
        acb_neg(log_z, z);
        acb_log(log_z, log_z, prec);

        acb_const_pi(t, prec);
        acb_mul_2exp_si(t, t, 1);
        acb_mul_onei(t, t);

        if (arf_sgn(arb_midref(acb_realref(z))) >= 0)
            acb_add(log_z, log_z, t, prec);
        else
            acb_sub(log_z, log_z, t, prec);
    }

    arb_one(acb_realref(t));
    is_real = acb_is_real(z) && acb_is_real(s) && acb_is_real(a) &&
        arb_is_positive(acb_realref(a)) && arb_lt(acb_realref(z), acb_realref(t));

    acb_set(param + 0, z);
    acb_set(param + 1, s);
    acb_set(param + 2, a);

    acb_one(N);
    mag_inf(tail_bound);

    /* todo: good relative magnitude */
    acb_one(t);
    acb_get_mag(abs_tol, t);
    mag_mul_2exp_si(abs_tol, abs_tol, -prec);
    acb_zero(t);

    for (i = 1; i < prec; i++)
    {
        acb_one(N);
        acb_mul_2exp_si(N, N, i);
        integral_tail(tail_bound, z, log_z, s, a, acb_realref(N), 53);
        if (mag_cmp(tail_bound, abs_tol) < 0)
            break;
    }

    rel_goal = prec;
    acb_calc_integrate_opt_init(options);

    arb_get_mag_lower(log_z_im_lower, acb_imagref(log_z));

    /* todo: with several integrals, add the errors to tolerances */

    if (acb_is_int(s) && arb_is_positive(acb_realref(s)))
    {
        /* integer s > 0: use direct integral */

        /* no need to avoid pole near path */
        if (arb_is_negative(acb_realref(log_z)) || mag_cmp_2exp_si(log_z_im_lower, -2) > 0)
        {
            acb_zero(xa);
            acb_set(xb, N);
            acb_calc_integrate(t, integrand, param, xa, xb, rel_goal, abs_tol, options, prec);
        }
        /* take detour through upper plane */
        else if (arb_is_nonpositive(acb_imagref(log_z)))
        {
            acb_zero(xa);
            acb_onei(xb);
            acb_calc_integrate(t, integrand, param, xa, xb, rel_goal, abs_tol, options, prec);

            acb_swap(xa, xb);
            acb_add(xb, xa, N, prec);
            acb_calc_integrate(u, integrand, param, xa, xb, rel_goal, abs_tol, options, prec);
            acb_add(t, t, u, prec);

            acb_swap(xa, xb);
            acb_set(xb, N);
            acb_calc_integrate(u, integrand, param, xa, xb, rel_goal, abs_tol, options, prec);
            acb_add(t, t, u, prec);
        }
        /* take detour through lower plane */
        else if (arb_is_positive(acb_imagref(log_z)))
        {
            acb_zero(xa);
            acb_onei(xb);
            acb_conj(xb, xb);
            acb_calc_integrate(t, integrand, param, xa, xb, rel_goal, abs_tol, options, prec);

            acb_swap(xa, xb);
            acb_add(xb, xa, N, prec);
            acb_calc_integrate(u, integrand, param, xa, xb, rel_goal, abs_tol, options, prec);
            acb_add(t, t, u, prec);

            acb_swap(xa, xb);
            acb_set(xb, N);
            acb_calc_integrate(u, integrand, param, xa, xb, rel_goal, abs_tol, options, prec);
            acb_add(t, t, u, prec);
        }
        else
        {
            /* todo: compute union of both branches? */
            acb_indeterminate(t);
        }

        acb_add_error_mag(t, tail_bound);

        if (is_real && acb_is_finite(t))
            arb_zero(acb_imagref(t));

        acb_rgamma(u, s, prec);
        acb_mul(res, t, u, prec);
    }
    else
    {
        arb_t left, right, bottom, top;
        acb_t residue;
        mag_t rm, im;

        arb_init(left);
        arb_init(right);
        arb_init(bottom);
        arb_init(top);
        acb_init(residue);
        mag_init(rm);
        mag_init(im);

        arb_get_mag_lower(rm, acb_realref(log_z));
        arb_get_mag_lower(im, acb_imagref(log_z));

        if (arb_is_negative(acb_realref(log_z)) && mag_cmp_2exp_si(rm, -1) > 0)
        {
            /* re(log(z)) < -0.5 - pole certainly excluded */
            /* left = min(|re(log(z))| / 2, 1) */
            arf_set_mag(arb_midref(left), rm);
            arf_mul_2exp_si(arb_midref(left), arb_midref(left), -1);
            if (arf_cmpabs_2exp_si(arb_midref(left), 0) > 0)
                arb_one(left);

            arb_set(right, left);
            arb_set(top, left);
            arb_set(bottom, left);
        }
        else if (mag_cmp_2exp_si(im, -1) > 0)
        {
            /* im(log(z)) > 0.5 - pole certainly excluded */
            arf_set_mag(arb_midref(left), im);
            arf_mul_2exp_si(arb_midref(left), arb_midref(left), -1);
            if (arf_cmpabs_2exp_si(arb_midref(left), 0) > 0)
                arb_one(left);

            arb_set(right, left);
            arb_set(top, left);
            arb_set(bottom, left);
        }
        else
        {
            /* residue = (-log(z))^s / log(z) / z^a */
            acb_neg(residue, log_z);
            acb_pow(residue, residue, s, prec);
            acb_div(residue, residue, log_z, prec);
            acb_pow(t, z, a, prec);
            acb_div(residue, residue, t, prec);

            /* |im(log(z)) + 1 */
            arb_get_mag(im, acb_imagref(log_z));

            /* containing a unique pole is uncertain -- error out */
            if (mag_cmp_2exp_si(im, 1) > 0)
            {
                acb_indeterminate(res);
                goto cleanup1;
            }
            else
            {
                arf_set_mag(arb_midref(top), im);
                arb_add_ui(top, top, 1, prec);
                arb_set(bottom, top);

                /* max(0, -re(log(z))) + 1 */
                arb_get_lbound_arf(arb_midref(left), acb_realref(log_z), prec);
                arb_neg(left, left);
                if (arf_sgn(arb_midref(left)) < 0)
                    arb_one(left);
                else
                    arb_add_ui(left, left, 1, prec);

                /* max(0, re(log(z))) + 1 */
                arb_get_ubound_arf(arb_midref(right), acb_realref(log_z), prec);
                if (arf_sgn(arb_midref(right)) < 0)
                    arb_one(right);
                else
                    arb_add_ui(right, right, 1, prec);
            }
        }

        arb_neg(left, left);
        arb_neg(bottom, bottom);

        acb_zero(t);

        /* w = (-1)^(s-1) */
        acb_sub_ui(w, s, 1, prec);
        acb_exp_pi_i(w, w, prec);

        if (is_real)
        {
            /* right -> top right */
            arb_set(acb_realref(xa), right); arb_zero(acb_imagref(xa));
            arb_set(acb_realref(xb), right); arb_set(acb_imagref(xb), top);
            acb_calc_integrate(u, integrand, param, xa, xb, rel_goal, abs_tol, options, prec);
            acb_div(u, u, w, prec);
            acb_add(t, t, u, prec);

            /* top right -> top left */
            arb_set(acb_realref(xa), right); arb_set(acb_imagref(xa), top);
            arb_set(acb_realref(xb), left); arb_set(acb_imagref(xb), top);
            acb_calc_integrate(u, integrand, param, xa, xb, rel_goal, abs_tol, options, prec);
            acb_div(u, u, w, prec);
            acb_add(t, t, u, prec);

            /* top left -> left */
            arb_set(acb_realref(xa), left); arb_set(acb_imagref(xa), top);
            arb_set(acb_realref(xb), left); arb_zero(acb_imagref(xb));
            acb_calc_integrate(u, integrand2, param, xa, xb, rel_goal, abs_tol, options, prec);
            acb_add(t, t, u, prec);

            /* 2 * imaginary part */
            arb_zero(acb_realref(t));
            acb_mul_2exp_si(t, t, 1);
        }
        else
        {
            /* right -> top right */
            arb_set(acb_realref(xa), right); arb_zero(acb_imagref(xa));
            arb_set(acb_realref(xb), right); arb_set(acb_imagref(xb), top);
            acb_calc_integrate(u, integrand, param, xa, xb, rel_goal, abs_tol, options, prec);
            acb_div(u, u, w, prec);
            acb_add(t, t, u, prec);

            /* top right -> top left */
            arb_set(acb_realref(xa), right); arb_set(acb_imagref(xa), top);
            arb_set(acb_realref(xb), left); arb_set(acb_imagref(xb), top);
            acb_calc_integrate(u, integrand, param, xa, xb, rel_goal, abs_tol, options, prec);
            acb_div(u, u, w, prec);
            acb_add(t, t, u, prec);

            /* top left -> bottom left */
            arb_set(acb_realref(xa), left); arb_set(acb_imagref(xa), top);
            arb_set(acb_realref(xb), left); arb_set(acb_imagref(xb), bottom);
            acb_calc_integrate(u, integrand2, param, xa, xb, rel_goal, abs_tol, options, prec);
            acb_add(t, t, u, prec);

            /* bottom left -> bottom right */
            arb_set(acb_realref(xa), left); arb_set(acb_imagref(xa), bottom);
            arb_set(acb_realref(xb), right); arb_set(acb_imagref(xb), bottom);
            acb_calc_integrate(u, integrand, param, xa, xb, rel_goal, abs_tol, options, prec);
            acb_mul(u, u, w, prec);
            acb_add(t, t, u, prec);

            /* bottom right -> right */
            arb_set(acb_realref(xa), right); arb_set(acb_imagref(xa), bottom);
            arb_set(acb_realref(xb), right); arb_zero(acb_imagref(xb));
            acb_calc_integrate(u, integrand, param, xa, xb, rel_goal, abs_tol, options, prec);
            acb_mul(u, u, w, prec);
            acb_add(t, t, u, prec);
        }

        /* right -> infinity */
        arb_set(acb_realref(xa), right); arb_zero(acb_imagref(xa));
        arb_set(acb_realref(xb), acb_realref(N)); arb_zero(acb_imagref(xb));
        acb_calc_integrate(u, integrand, param, xa, xb, rel_goal, abs_tol, options, prec);

        acb_add_error_mag(u, tail_bound);

        /* (w - 1/w) */
        acb_inv(v, w, prec);
        acb_sub(v, w, v, prec);
        acb_addmul(t, u, v, prec);

        /* -gamma(1-s) * (t / (2 pi i) + residue) */
        acb_const_pi(u, prec);
        acb_mul_onei(u, u);
        acb_mul_2exp_si(u, u, 1);
        acb_div(t, t, u, prec);
        acb_add(t, t, residue, prec);

        acb_sub_ui(u, s, 1, prec);
        acb_neg(u, u);
        acb_gamma(u, u, prec);
        acb_neg(u, u);

        acb_mul(res, t, u, prec);

        if (is_real)
            arb_zero(acb_imagref(res));

cleanup1:
        arb_clear(left);
        arb_clear(right);
        arb_clear(bottom);
        arb_clear(top);
        acb_clear(residue);
        mag_clear(rm);
        mag_clear(im);
    }

    mag_clear(log_z_im_lower);
    acb_clear(xa);
    acb_clear(xb);

    acb_clear(log_z);
    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
    acb_clear(w);
    acb_clear(zero);
    acb_clear(N);
    mag_clear(abs_tol);
    mag_clear(tail_bound);
    _acb_vec_clear(param, 3);
}

void
acb_dirichlet_lerch_phi_integral(acb_t res, const acb_t z, const acb_t s, const acb_t a, slong prec)
{
    if (!acb_is_finite(z) || !acb_is_finite(s) || !acb_is_finite(a) ||
        acb_contains_zero(z) ||
        (arb_contains_si(acb_realref(z), 1) && arb_contains_zero(acb_imagref(z))))
    {
        acb_indeterminate(res);
        return;
    }

    if (acb_contains_int(a) && !arb_is_positive(acb_realref(a)))
    {
        if (!(acb_is_int(s) && arb_is_nonpositive(acb_realref(s))))
        {
            acb_indeterminate(res);
            return;
        }
    }

    /* wide intervals will just cause numerical integration
       to converge slowly, clogging unit tests */
    if (acb_rel_accuracy_bits(z) < 1 || acb_rel_accuracy_bits(s) < 1 || acb_rel_accuracy_bits(a) < 1)
    {
        acb_indeterminate(res);
        return;
    }

    if (arf_cmp_si(arb_midref(acb_realref(a)), -2 * prec) < 0)
    {
        acb_indeterminate(res);
        return;
    }

    if (arf_cmp_si(arb_midref(acb_realref(a)), 2) < 0)
    {
        slong n, N, wp;
        acb_t t, u, sum, negs;

        N = 2 - arf_get_si(arb_midref(acb_realref(a)), ARF_RND_FLOOR);
        wp = prec + 10;

        acb_init(t);
        acb_init(u);
        acb_init(sum);
        acb_init(negs);

        acb_one(t);
        acb_neg(negs, s);

        for (n = 0; n < N; n++)
        {
            if (n >= 1)
            {
                if (n % 8 == 0 && !acb_is_real(z))
                    acb_pow_ui(t, z, n, wp);
                else
                    acb_mul(t, t, z, wp);
            }

            acb_add_ui(u, a, n, wp);
            acb_pow(u, u, negs, wp);
            acb_mul(u, t, u, wp);
            acb_add(sum, sum, u, wp);
        }

        acb_add_ui(t, a, n, wp);
        _acb_dirichlet_lerch_phi_integral(u, z, s, t, prec);
        acb_pow_ui(t, z, n, prec);
        acb_mul(u, u, t, prec);

        acb_add(res, u, sum, prec);

        acb_clear(t);
        acb_clear(u);
        acb_clear(sum);
        acb_clear(negs);
    }
    else
    {
        _acb_dirichlet_lerch_phi_integral(res, z, s, a, prec);
    }
}
