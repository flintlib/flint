/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static int
acb_theta_ql_a0_start(acb_ptr th, acb_srcptr t, acb_srcptr z, arb_srcptr d0,
    arb_srcptr d, const acb_t f, const acb_mat_t tau, slong nb_steps, slong s,
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
    acb_t c;
    int res;

    acb_mat_init(w, g, g);
    x = _acb_vec_init(g);
    u = _acb_vec_init(g);
    zero = _acb_vec_init(g);
    new_d0 = _arb_vec_init(n);
    new_d = _arb_vec_init(n);
    acb_init(c);

    acb_mat_scalar_mul_2exp_si(w, tau, nb_steps);
    _acb_vec_scalar_mul_2exp_si(u, t, g, nb_steps);
    _acb_vec_scalar_mul_2exp_si(x, z, g, nb_steps);
    _arb_vec_scalar_mul_2exp_si(new_d0, d0, n, nb_steps);
    _arb_vec_scalar_mul_2exp_si(new_d, d, n, nb_steps);
    acb_mul_2exp_si(c, f, nb_steps);
    acb_exp_pi_i(c, c, prec);

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
        _acb_vec_scalar_mul(th + nbt * n, th + nbt * n, nbt * n, c, prec);
    }

    acb_mat_clear(w);
    _acb_vec_clear(x, g);
    _acb_vec_clear(u, g);
    _acb_vec_clear(zero, g);
    _arb_vec_clear(new_d0, n);
    _arb_vec_clear(new_d, n);
    acb_clear(c);
    return res;
}

static void
acb_theta_ql_a0_step(acb_ptr th, acb_srcptr rts, arb_srcptr d0, arb_srcptr d,
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
        _acb_vec_set(rts + j * n, rts + j * nb_steps * n + k * n, n);
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
    acb_ptr x, rts;
    acb_t f, c;
    slong k;
    int res = 1;

    x = _acb_vec_init(g);
    rts = _acb_vec_init(nbz * nbr * n * nb_steps);
    acb_init(f);
    acb_init(c);

    acb_theta_ql_log_rescale(f, z, tau, prec);

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
        acb_neg(c, f);
        acb_exp_pi_i(c, c, prec);
        _acb_vec_scalar_mul(th + nbt * n, th + nbt * n, n * nbt, c, prec);
    }

    _acb_vec_clear(x, g);
    _acb_vec_clear(rts, nbz * nbr * n * nb_steps);
    acb_clear(f);
    acb_clear(c);
    return res;
}
