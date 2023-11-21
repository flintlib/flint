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
#include "acb_poly.h"

void
acb_dirichlet_l_hurwitz(acb_t res, const acb_t s,
    const acb_dirichlet_hurwitz_precomp_t precomp,
    const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
{
    ulong order, chin, mult;
    acb_t t, u, a, w;
    dirichlet_char_t cn;
    acb_dirichlet_roots_t roots;
    int deflate;

    /* remove pole in Hurwitz zeta at s = 1 */
    deflate = 0;
    if (acb_is_one(s))
    {
        if (dirichlet_char_is_principal(G, chi))
        {
            acb_indeterminate(res);
            return;
        }
        deflate = 1;
    }

    dirichlet_char_init(cn, G);
    acb_init(t);
    acb_init(u);
    acb_init(a);
    acb_init(w);

    dirichlet_char_one(cn, G);
    acb_zero(t);

    prec += n_clog(G->phi_q, 2);

    order = dirichlet_order_char(G, chi);
    mult = G->expo / order;
    acb_dirichlet_roots_init(roots, order, dirichlet_group_size(G), prec);

    do {
        chin = dirichlet_pairing_char(G, chi, cn) / mult;

        if (precomp == NULL)
        {
            acb_set_ui(a, cn->n);
            acb_div_ui(a, a, G->q, prec);

            if (deflate == 0)
                acb_hurwitz_zeta(u, s, a, prec);
            else
                _acb_poly_zeta_cpx_series(u, s, a, 1, 1, prec);
        }
        else
        {
            acb_dirichlet_hurwitz_precomp_eval(u, precomp, cn->n, G->q, prec);
        }

        acb_dirichlet_root(w, roots, chin, prec);
        acb_addmul(t, u, w, prec);

    } while (dirichlet_char_next(cn, G) >= 0);

    acb_set_ui(u, G->q);
    acb_neg(a, s);
    acb_pow(u, u, a, prec);
    acb_mul(res, t, u, prec);

    dirichlet_char_clear(cn);

    acb_dirichlet_roots_clear(roots);
    acb_clear(t);
    acb_clear(u);
    acb_clear(a);
    acb_clear(w);
}

