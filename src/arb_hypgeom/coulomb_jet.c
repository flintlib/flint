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
arb_hypgeom_coulomb_jet(arb_ptr F, arb_ptr G, const arb_t l, const arb_t eta, const arb_t z, slong len, slong prec)
{
    acb_ptr tmp, tmpF, tmpG;
    slong k;

    if (len <= 0)
        return;

    if (len == 1)
    {
        arb_hypgeom_coulomb(F, G, l, eta, z, prec);
        return;
    }

    tmp = _acb_vec_init(3);
    tmpF = _acb_vec_init(len);
    tmpG = _acb_vec_init(len);

    acb_set_arb(tmp, l);
    acb_set_arb(tmp + 1, eta);
    acb_set_arb(tmp + 2, z);

    acb_hypgeom_coulomb_jet(F ? tmpF : NULL, G ? tmpG : NULL,
        NULL, NULL, tmp, tmp + 1, tmp + 2, len, prec);

    if (F != NULL)
    {
        if (acb_is_real(tmpF))
            for (k = 0; k < len; k++)
                arb_set(F + k, acb_realref(tmpF + k));
        else
            _arb_vec_indeterminate(F, len);
    }

    if (G != NULL)
    {
        if (acb_is_real(tmpG))
            for (k = 0; k < len; k++)
                arb_set(G + k, acb_realref(tmpG + k));
        else
            _arb_vec_indeterminate(G, len);
    }

    _acb_vec_clear(tmpF, len);
    _acb_vec_clear(tmpG, len);
    _acb_vec_clear(tmp, 3);
}

