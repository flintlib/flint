/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_dirichlet.h"

void
_acb_dirichlet_hardy_z_series(acb_ptr res, acb_srcptr s, slong slen,
    const dirichlet_group_t G, const dirichlet_char_t chi,
    slong len, slong prec)
{
    slen = FLINT_MIN(slen, len);

    if (len == 0)
        return;

    if (slen == 1)
    {
        acb_dirichlet_hardy_z(res, s, G, chi, 1, prec);
        _acb_vec_zero(res + 1, len - 1);
    }
    else
    {
        acb_ptr t, u;
        t = _acb_vec_init(len);
        u = _acb_vec_init(slen);

        acb_dirichlet_hardy_z(t, s, G, chi, len, prec);

        /* compose with nonconstant part */
        acb_zero(u);
        _acb_vec_set(u + 1, s + 1, slen - 1);
        _acb_poly_compose_series(res, t, len, u, slen, len, prec);

        _acb_vec_clear(t, len);
        _acb_vec_clear(u, slen);
    }
}

void
acb_dirichlet_hardy_z_series(acb_poly_t res, const acb_poly_t s,
    const dirichlet_group_t G, const dirichlet_char_t chi,
    slong len, slong prec)
{
    if (len == 0)
    {
        acb_poly_zero(res);
        return;
    }

    acb_poly_fit_length(res, len);

    if (s->length == 0)
    {
        acb_t t;
        acb_init(t);
        _acb_dirichlet_hardy_z_series(res->coeffs, t, 1, G, chi, len, prec);
        acb_clear(t);
    }
    else
    {
        _acb_dirichlet_hardy_z_series(res->coeffs, s->coeffs, s->length, G, chi, len, prec);
    }

    _acb_poly_set_length(res, len);
    _acb_poly_normalise(res);
}

