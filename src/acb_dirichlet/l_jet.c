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

/* todo: move implemetation to the acb_dirichlet module */
void _acb_poly_zeta_cpx_reflect(acb_ptr t, const acb_t h,
    const acb_t a, int deflate, slong len, slong prec);

void
acb_dirichlet_l_jet(acb_ptr res, const acb_t s,
    const dirichlet_group_t G, const dirichlet_char_t chi,
    int deflate, slong len, slong prec)
{
    ulong order, chin, mult, phi;
    acb_t a, w;
    acb_ptr t, u;
    dirichlet_char_t cn;
    acb_dirichlet_roots_t roots;
    int deflate_hurwitz;

    if (len <= 0)
        return;

    /* special-case Riemann zeta */
    if (G == NULL || G->q == 1)
    {
        if (len == 1 && !deflate)
            acb_dirichlet_zeta(res, s, prec);
        else
            acb_dirichlet_zeta_jet(res, s, deflate, len, prec);
        return;
    }

    if (len == 1 && !(deflate && dirichlet_char_is_principal(G, chi)))
    {
        acb_dirichlet_l(res, s, G, chi, prec);
        return;
    }

    if (dirichlet_char_is_principal(G, chi))
        deflate_hurwitz = deflate;
    else
        deflate_hurwitz = acb_is_one(s);

    dirichlet_char_init(cn, G);
    t = _acb_vec_init(len);
    u = _acb_vec_init(len + 2);
    acb_init(a);
    acb_init(w);

    dirichlet_char_one(cn, G);

    prec += n_clog(G->phi_q, 2);

    order = dirichlet_order_char(G, chi);
    mult = G->expo / order;
    acb_dirichlet_roots_init(roots, order, dirichlet_group_size(G), prec);

    phi = 0;
    do
    {
        chin = dirichlet_pairing_char(G, chi, cn) / mult;
        acb_set_ui(a, cn->n);
        acb_div_ui(a, a, G->q, prec);
        _acb_poly_zeta_cpx_series(u, s, a, deflate_hurwitz, len, prec);
        acb_dirichlet_root(w, roots, chin, prec);
        _acb_vec_scalar_addmul(t, u, len, w, prec);
        phi++;
    }
    while (dirichlet_char_next(cn, G) >= 0);

    if (dirichlet_char_is_principal(G, chi) && deflate)
    {
        /* res = t * q^(-(s+x)) + [phi(q) * (q^(-(s+x)) - q^-1) / ((s+x)-1)] */
        if (acb_is_one(s))
        {
            acb_set_ui(a, G->q);
            _acb_poly_acb_invpow_cpx(u, a, s, len + 1, prec);
            _acb_poly_mullow(res, t, len, u, len, len, prec);
            acb_set_ui(u, phi);
            _acb_vec_scalar_addmul(res, u + 1, len, u, prec);
        }
        else
        {
            acb_sub_ui(u, s, 1, prec);
            acb_one(u + 1);

            acb_set_ui(a, G->q);
            _acb_poly_acb_invpow_cpx(u + 2, a, s, len, prec);
            _acb_poly_mullow(res, t, len, u + 2, len, len, prec);

            acb_inv(a, a, prec);
            acb_sub(u + 2, u + 2, a, prec);
            _acb_poly_div_series(t, u + 2, len, u, 2, len, prec);

            acb_set_ui(u, phi);
            _acb_vec_scalar_addmul(res, t, len, u, prec);
        }
    }
    else
    {
        /* res = t * q^(-(s+x)) */
        acb_set_ui(a, G->q);
        _acb_poly_acb_invpow_cpx(u, a, s, len, prec);
        _acb_poly_mullow(res, t, len, u, len, len, prec);
    }

    dirichlet_char_clear(cn);

    acb_dirichlet_roots_clear(roots);
    _acb_vec_clear(t, len);
    _acb_vec_clear(u, len + 2);
    acb_clear(a);
    acb_clear(w);
}

