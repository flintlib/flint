/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_ode.h"
#include "acb_poly.h"
#include "acb_types.h"

void
acb_ode_indicial_polynomial_from_exponents(acb_poly_t ind,
                                           const acb_ode_exponents_t expos,
                                           slong prec)
{
    slong len = acb_ode_exponents_length(expos);

    acb_struct * roots = _acb_vec_init(len);

    acb_struct * r = roots;
    for (slong g = 0; g < expos->ngroups; g++)
    {
        acb_ode_group_struct * group = expos->groups + g;

        for (slong s = 0; s < group->nshifts; s++)
        {
            acb_add_si(r, group->leader, group->shifts[s].n, prec);
            r++;

            for (slong m = 1; m < expos->groups[g].shifts[s].mult; m++)
            {
                acb_set(r, r-1);
                r++;
            }
        }
    }

    acb_poly_product_roots(ind, roots, len, prec);

    _acb_vec_clear(roots, len);
}
