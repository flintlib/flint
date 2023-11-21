/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dft.h"

acb_dft_step_ptr
_acb_dft_steps_prod(slong * cyc, slong num, slong prec)
{
    slong i, len;
    acb_dft_step_ptr s;
    s = flint_malloc(num * sizeof(acb_dft_step_struct));

    len = 1;
    for (i = 0; i < num; i++)
        len *= cyc[i];

    for (i = 0; i < num; i++)
    {
        slong m, M;
        m = cyc[i];
        M = (len /= m);
        s[i].m = m;
        s[i].M = M;
        s[i].dv = M;
        s[i].dz = 0;
        s[i].z = NULL;
        _acb_dft_precomp_init(s[i].pre, M, NULL, 0, m, prec);
    }

    return s;
}

void
acb_dft_prod_clear(acb_dft_prod_t t)
{
    slong i;
    for (i = 0; i < t->num; i++)
        acb_dft_precomp_clear(t->cyc[i].pre);
    flint_free(t->cyc);
}


void
acb_dft_prod_precomp(acb_ptr w, acb_srcptr v, const acb_dft_prod_t prod, slong prec)
{
    if (prod->num >= 1)
        acb_dft_step(w, v, prod->cyc, prod->num, prec);
}

void
acb_dft_prod(acb_ptr w, acb_srcptr v, slong * cyc, slong num, slong prec)
{
    acb_dft_prod_t t;
    acb_dft_prod_init(t, cyc, num, prec);
    acb_dft_prod_precomp(w, v, t, prec);
    acb_dft_prod_clear(t);
}
