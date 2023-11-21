/*
    Copyright (C) 2016 Fredrik Johansson
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

#ifdef __GNUC__
# define log __builtin_log
#else
# include <math.h>
#endif

void
acb_dirichlet_l_general(acb_t res, const acb_t s,
    const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
{
    /* this cutoff is probably too conservative when q is large */
    if (arf_cmp_d(arb_midref(acb_realref(s)), 8 + 0.5 * prec / log(prec)) >= 0)
    {
        acb_dirichlet_l_euler_product(res, s, G, chi, prec);
    }
    else
    {
        slong wp = prec + n_clog(G->phi_q, 2);
        acb_dirichlet_hurwitz_precomp_t pre;
        acb_dirichlet_hurwitz_precomp_init_num(pre, s, acb_is_one(s), G->phi_q, wp);
        acb_dirichlet_l_hurwitz(res, s, pre, G, chi, prec);
        acb_dirichlet_hurwitz_precomp_clear(pre);
    }
}

void
acb_dirichlet_l(acb_t res, const acb_t s,
    const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
{
    if (!acb_is_finite(s))
    {
        acb_indeterminate(res);
    }
    else if (G == NULL || G->q == 1)
    {
        acb_dirichlet_zeta(res, s, prec);
    }
    else if (dirichlet_char_is_primitive(G, chi) &&
        (arf_cmp_d(arb_midref(acb_realref(s)), -0.5) < 0 ||
            (G->q != 1 && dirichlet_parity_char(G, chi) == 0 &&
                arf_cmpabs_d(arb_midref(acb_imagref(s)), 0.125) < 0 &&
                arf_cmp_d(arb_midref(acb_realref(s)), 0.125) < 0)))
    {
        /* use functional equation */
        acb_t t, u, v;
        int parity;
        ulong q;

        parity = dirichlet_parity_char(G, chi);
        q = G->q;

        acb_init(t);
        acb_init(u);
        acb_init(v);

        /* gamma((1-s+p)/2) / gamma((s+p)/2) */
        acb_add_ui(t, s, parity, prec);
        acb_mul_2exp_si(t, t, -1);
        acb_rgamma(t, t, prec);

        if (!acb_is_zero(t))  /* assumes q != 1 when s = 0 */
        {
            acb_neg(u, s);
            acb_add_ui(u, u, 1 + parity, prec);
            acb_mul_2exp_si(u, u, -1);
            acb_gamma(u, u, prec);
            acb_mul(t, t, u, prec);

            /* epsilon */
            acb_dirichlet_root_number(u, G, chi, prec);
            acb_mul(t, t, u, prec);

            /* (pi/q)^(s-1/2) */
            acb_const_pi(u, prec);
            acb_div_ui(u, u, q, prec);
            acb_set_d(v, -0.5);
            acb_add(v, v, s, prec);
            acb_pow(u, u, v, prec);
            acb_mul(t, t, u, prec);

            acb_sub_ui(u, s, 1, prec);
            acb_neg(u, u);
            acb_conj(u, u);
            acb_dirichlet_l_general(u, u, G, chi, prec);
            acb_conj(u, u);
            acb_mul(t, t, u, prec);

            if (dirichlet_char_is_real(G, chi) && acb_is_real(s))
                arb_zero(acb_imagref(t));
        }

        acb_set(res, t);

        acb_clear(t);
        acb_clear(u);
        acb_clear(v);
    }
    else
    {
        acb_dirichlet_l_general(res, s, G, chi, prec);
    }
}

