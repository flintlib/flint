/*
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
acb_dirichlet_l_vec_hurwitz(acb_ptr res, const acb_t s,
    const acb_dirichlet_hurwitz_precomp_t precomp,
    const dirichlet_group_t G, slong prec)
{
    acb_t a, qs;
    acb_ptr zeta, z;
    dirichlet_char_t cn;
    int deflate;

    /* remove pole in Hurwitz zeta at s = 1 */
    deflate = acb_is_one(s);

    dirichlet_char_init(cn, G);
    acb_init(qs);
    acb_init(a);

    prec += n_clog(G->phi_q, 2);

    acb_set_ui(qs, G->q);
    acb_neg(a, s);
    acb_pow(qs, qs, a, prec);

    zeta = z = _acb_vec_init(G->phi_q);
    dirichlet_char_one(cn, G);
    do {
        if (precomp == NULL)
        {
            acb_set_ui(a, cn->n);
            acb_div_ui(a, a, G->q, prec);

            if (deflate == 0)
                acb_hurwitz_zeta(z, s, a, prec);
            else
                _acb_poly_zeta_cpx_series(z, s, a, 1, 1, prec);
        }
        else
        {
            acb_dirichlet_hurwitz_precomp_eval(z, precomp, cn->n, G->q, prec);
        }

        acb_mul(z, z, qs, prec);
        acb_conj(z, z);

        z++;
    } while (dirichlet_char_next(cn, G) >= 0);

    acb_dirichlet_dft_index(res, zeta, G, prec);

    {
        slong k;
        for (k = 0; k < G->phi_q; k++)
            acb_conj(res + k, res + k);
    }

    /* restore pole for the principal character */
    if (deflate)
        acb_indeterminate(res);

    dirichlet_char_clear(cn);
    _acb_vec_clear(zeta, G->phi_q);
    acb_clear(qs);
    acb_clear(a);
}

