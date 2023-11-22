/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

static int
acb_theta_ql_roots_1(acb_ptr rts, acb_srcptr z, arb_srcptr d,
    const arb_t f, const acb_mat_t tau, slong nb_steps, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_mat_t w;
    acb_ptr x;
    arb_t c, h;
    slong hprec, guard;
    slong k, a;
    int res = 1;

    acb_mat_init(w, g, g);
    x = _acb_vec_init(g);
    arb_init(c);
    arb_init(h);

    for (k = 0; (k < nb_steps) && res; k++)
    {
        acb_mat_scalar_mul_2exp_si(w, tau, k);
        _acb_vec_scalar_mul_2exp_si(x, z, g, k);
        arb_mul_2exp_si(c, f, k);
        arb_exp(c, c, prec);

        for (a = 0; (a < n) && res; a++)
        {
            arb_mul_2exp_si(h, &d[a], k);
            res = 0;
            for (guard = 16; (guard <= prec) && !res; guard += 16)
            {
                hprec = guard + acb_theta_dist_addprec(h);
                acb_theta_naive_fixed_ab(&rts[k * n + a], a << g, x, 1, w, hprec);
                if (acb_is_finite(&rts[k * n + a]) && !acb_contains_zero(&rts[k * n + a]))
                {
                    res = 1;
                }
            }
        }

        _acb_vec_scalar_mul_arb(rts + k * n, rts + k * n, n, c, prec);
    }

    acb_mat_clear(w);
    _acb_vec_clear(x, g);
    arb_clear(c);
    arb_clear(h);
    return res;
}

static int
acb_theta_ql_roots_3(acb_ptr rts, acb_srcptr t, acb_srcptr z, arb_srcptr d,
    const acb_mat_t tau, slong nb_steps, slong guard, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    int has_t = !_acb_vec_is_zero(t, g);
    arb_mat_t Yinv;
    acb_ptr x;
    arb_ptr y, w;
    arb_t f, pi;
    slong k;
    int res = 1;

    arb_mat_init(Yinv, g, g);
    x = _acb_vec_init(g);
    y = _arb_vec_init(g);
    w = _arb_vec_init(g);
    arb_init(f);
    arb_init(pi);

    acb_siegel_yinv(Yinv, tau, prec);
    _acb_vec_get_imag(y, z, g);
    arb_mat_vector_mul_col(w, Yinv, y, prec);
    arb_dot(f, NULL, 1, y, 1, w, 1, g, prec);
    arb_const_pi(pi, prec);
    arb_mul(f, f, pi, prec);

    if (!has_t)
    {
        res = acb_theta_ql_roots_1(rts, z, d, f, tau, nb_steps, guard);
    }
    else
    {
        for (k = 1; (k < 3) && res; k++)
        {
            _acb_vec_scalar_mul_ui(x, t, g, k, prec);
            _acb_vec_add(x, x, z, g, prec);
            res = acb_theta_ql_roots_1(rts + (k - 1) * nb_steps * n, x, d,
                f, tau, nb_steps, guard);
        }
    }

    arb_mat_clear(Yinv);
    _acb_vec_clear(x, g);
    _arb_vec_clear(y, g);
    _arb_vec_clear(w, g);
    arb_clear(f);
    arb_clear(pi);
    return res;
}

static int
acb_theta_ql_roots(acb_ptr rts, acb_srcptr t, acb_srcptr z, arb_srcptr d0,
    arb_srcptr d, const acb_mat_t tau, slong nb_steps, slong guard, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    int hasz = !_acb_vec_is_zero(z, g);
    int hast = !_acb_vec_is_zero(t, g);
    slong nbt = (hast ? 2 : 1);
    acb_ptr x;
    int res;

    x = _acb_vec_init(g);

    res = acb_theta_ql_roots_3(rts, t, x, d0, tau, nb_steps, guard, prec);
    if (res && hasz)
    {
        res = acb_theta_ql_roots_3(rts + nbt * n * nb_steps, t, z, d, tau,
            nb_steps, guard, prec);
    }

    _acb_vec_clear(x, g);
    return res;
}

static int
acb_theta_ql_a0_start(acb_ptr th, acb_srcptr t, acb_srcptr z, arb_srcptr d0,
    arb_srcptr d, const arb_t f, const acb_mat_t tau, slong nb_steps, slong s,
    slong guard, slong prec, acb_theta_ql_worker_t worker)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    int hast = !_acb_vec_is_zero(t, g);
    int hasz = !_acb_vec_is_zero(z, g);
    slong nbt = (hast ? 3 : 1);
    acb_mat_t w;
    acb_ptr x, u, zero;
    arb_ptr new_d0, new_d;
    arb_t c;
    int res;

    acb_mat_init(w, g, g);
    x = _acb_vec_init(g);
    u = _acb_vec_init(g);
    zero = _acb_vec_init(g);
    new_d0 = _arb_vec_init(n);
    new_d = _arb_vec_init(n);
    arb_init(c);

    acb_mat_scalar_mul_2exp_si(w, tau, nb_steps);
    _acb_vec_scalar_mul_2exp_si(u, t, g, nb_steps);
    _acb_vec_scalar_mul_2exp_si(x, z, g, nb_steps);
    _arb_vec_scalar_mul_2exp_si(new_d0, d0, n, nb_steps);
    _arb_vec_scalar_mul_2exp_si(new_d, d, n, nb_steps);
    arb_mul_2exp_si(c, f, nb_steps);
    arb_exp(c, c, prec);

    if (s > 0)
    {
        res = acb_theta_ql_a0_split(th, u, zero, new_d0, w, s, guard, prec, worker);
        if (res && hasz)
        {
            res = acb_theta_ql_a0_split(th + nbt * n, u, x, new_d, w, s, guard, prec, worker);
        }
    }
    else
    {
        res = acb_theta_ql_a0_naive(th, u, x, new_d0, new_d, w, guard, prec);
    }

    if (hasz)
    {
        _acb_vec_scalar_mul_arb(th + nbt * n, th + nbt * n, nbt * n, c, prec);
    }

    acb_mat_clear(w);
    _acb_vec_clear(x, g);
    _acb_vec_clear(u, g);
    _acb_vec_clear(zero, g);
    _arb_vec_clear(new_d0, n);
    _arb_vec_clear(new_d, n);
    arb_clear(c);
    return res;
}

static void
acb_theta_ql_step_1(acb_ptr res, acb_srcptr th0, acb_srcptr th, acb_srcptr rts,
    arb_srcptr d0, arb_srcptr d, slong g, slong prec)
{
    slong n = 1 << g;

    acb_theta_agm_mul_tight(res, th0, th, d0, d, g, prec);
    _acb_vec_scalar_mul_2exp_si(res, res, n, g);
    acb_theta_agm_sqrt(res, res, rts, n, prec);
}

static void
acb_theta_ql_step_2(acb_ptr res, acb_srcptr th0, acb_srcptr th, acb_srcptr rts,
    arb_srcptr d0, arb_srcptr d, slong g, slong prec)
{
    slong n = 1 << g;
    acb_ptr aux;

    aux = _acb_vec_init(3 * n);

    acb_theta_agm_mul_tight(aux + n, th0, th + n, d0, d, g, prec);
    acb_theta_agm_mul_tight(aux + 2 * n, th0, th + 2 * n, d0, d, g, prec);
    _acb_vec_scalar_mul_2exp_si(aux + n, aux + n, 2 * n, g);
    acb_theta_agm_sqrt(aux + n, aux + n, rts, 2 * n, prec);

    _acb_vec_set(res, aux, 3 * n);
    _acb_vec_clear(aux, 3 * n);
}


static void
acb_theta_ql_step_3(acb_ptr res, acb_srcptr th0, acb_srcptr th, acb_srcptr rts,
    arb_srcptr d0, arb_srcptr d, slong g, slong prec)
{
    slong n = 1 << g;
    acb_ptr aux;
    ulong a;

    aux = _acb_vec_init(3 * n);

    acb_theta_agm_mul_tight(aux + n, th0, th + n, d0, d, g, prec);
    acb_theta_agm_mul_tight(aux + 2 * n, th0, th + 2 * n, d0, d, g, prec);
    _acb_vec_scalar_mul_2exp_si(aux + n, aux + n, 2 * n, g);
    acb_theta_agm_sqrt(aux + n, aux + n, rts, 2 * n, prec);

    acb_theta_agm_mul_tight(aux, th0 + n, th + n, d0, d, g, prec);
    _acb_vec_scalar_mul_2exp_si(aux, aux, n, g);
    for (a = 0; a < n; a++)
    {
        acb_div(&aux[a], &aux[a], &aux[2 * n + a], prec);
    }

    _acb_vec_set(res, aux, 3 * n);
    _acb_vec_clear(aux, 3 * n);
}

static void
acb_theta_ql_a0_step(acb_ptr th, acb_srcptr all_rts, arb_srcptr d0, arb_srcptr d,
    slong k, slong nb_steps, int hast, int hasz, slong g, slong prec)
{
    slong n = 1 << g;
    arb_ptr new_d, new_d0;
    acb_ptr next;
    acb_ptr rts;
    slong nbt = (hast ? 3 : 1);
    slong nbr = (hast ? 2 : 1);
    slong nbz = (hasz ? 2 : 1);
    slong j;

    new_d = _arb_vec_init(n);
    new_d0 = _arb_vec_init(n);
    next = _acb_vec_init(nbz * nbt * n);
    rts = _acb_vec_init(nbr * nbz * n);

    _arb_vec_scalar_mul_2exp_si(new_d, d, n, k + 1);
    _arb_vec_scalar_mul_2exp_si(new_d0, d0, n, k + 1);
    for (j = 0; j < nbz * nbr; j++)
    {
        _acb_vec_set(rts + j * n, all_rts + j * nb_steps * n + k * n, n);
    }

    if (hast)
    {
        acb_theta_ql_step_3(next, th, th, rts, new_d0, new_d0, g, prec);
        if (hasz && (k == 0))
        {
            acb_theta_ql_step_3(next + nbt * n, th, th + nbt * n,
                rts + nbr * n, new_d0, new_d, g, prec);
        }
        else if (hasz)
        {
            acb_theta_ql_step_2(next + nbt * n, th, th + nbt * n,
                rts + nbr * n, new_d0, new_d, g, prec);
        }
    }
    else
    {
        acb_theta_ql_step_1(next, th, th, rts, new_d0, new_d0, g, prec);
        if (hasz)
        {
            acb_theta_ql_step_1(next + nbt * n, th, th + nbt * n,
                rts + nbr * n, new_d0, new_d, g, prec);
        }
    }
    _acb_vec_set(th, next, nbz * nbt * n);

    _arb_vec_clear(new_d, n);
    _arb_vec_clear(new_d0, n);
    _acb_vec_clear(next, nbz * nbt * n);
    _acb_vec_clear(rts, nbr * nbz * n);
}

int
acb_theta_ql_a0_steps(acb_ptr th, acb_srcptr t, acb_srcptr z, arb_srcptr d0,
    arb_srcptr d, const acb_mat_t tau, slong nb_steps, slong s,
    slong guard, slong prec, acb_theta_ql_worker_t worker)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    int hast = !_acb_vec_is_zero(t, g);
    int hasz = !_acb_vec_is_zero(z, g);
    slong nbt = (hast ? 3 : 1);
    slong nbr = (hast ? 2 : 1);
    slong nbz = (hasz ? 2 : 1);
    arb_mat_t Yinv;
    acb_ptr x, rts;
    arb_ptr y, w;
    arb_t f, c;
    slong k;
    int res = 1;

    arb_mat_init(Yinv, g, g);
    x = _acb_vec_init(g);
    y = _arb_vec_init(g);
    w = _arb_vec_init(g);
    rts = _acb_vec_init(nbz * nbr * n * nb_steps);
    arb_init(f);
    arb_init(c);

    acb_siegel_yinv(Yinv, tau, prec);
    _acb_vec_get_imag(y, z, g);
    arb_mat_vector_mul_col(w, Yinv, y, prec);
    arb_dot(f, NULL, 1, y, 1, w, 1, g, prec);
    arb_const_pi(c, prec);
    arb_mul(f, f, c, prec);

    res = acb_theta_ql_roots(rts, t, z, d0, d, tau, nb_steps, guard, prec);
    if (res)
    {
        res = acb_theta_ql_a0_start(th, t, z, d0, d,  f, tau, nb_steps, s,
            guard, prec, worker);
    }
    if (res)
    {
        for (k = nb_steps - 1; k >= 0; k--)
        {
            acb_theta_ql_a0_step(th, rts, d0, d, k, nb_steps, hast, hasz, g, prec);
        }
    }
    if (res && hasz)
    {
        arb_neg(f, f);
        arb_exp(c, f, prec);
        _acb_vec_scalar_mul_arb(th + nbt * n, th + nbt * n, n * nbt, c, prec);
    }

    arb_mat_clear(Yinv);
    _acb_vec_clear(x, g);
    _arb_vec_clear(y, g);
    _arb_vec_clear(w, g);
    _acb_vec_clear(rts, nbz * nbr * n * nb_steps);
    arb_clear(f);
    arb_clear(c);
    return res;
}
