/*
    Copyright (C) 2019 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"

void
arb_hypgeom_coulomb(arb_t F, arb_t G, const arb_t l, const arb_t eta, const arb_t z, slong prec)
{
    acb_ptr tmp;

    tmp = _acb_vec_init(5);

    acb_set_arb(tmp + 2, l);
    acb_set_arb(tmp + 3, eta);
    acb_set_arb(tmp + 4, z);

    acb_hypgeom_coulomb(F ? tmp : NULL, G ? tmp + 1 : NULL,
        NULL, NULL, tmp + 2, tmp + 3, tmp + 4, prec);

    if (F != NULL)
    {
        if (acb_is_real(tmp))
            arb_set(F, acb_realref(tmp));
        else
            arb_indeterminate(F);
    }

    if (G != NULL)
    {
        if (acb_is_real(tmp + 1))
            arb_set(G, acb_realref(tmp + 1));
        else
            arb_indeterminate(G);
    }

    _acb_vec_clear(tmp, 5);
}

