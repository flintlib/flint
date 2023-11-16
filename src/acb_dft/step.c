/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dft.h"

#define REORDER 0

void
acb_dft_step(acb_ptr w, acb_srcptr v, acb_dft_step_ptr cyc, slong num, slong prec)
{
    acb_dft_step_struct c;
    if (num == 0)
        return;
    c = cyc[0];
    if (num == 1)
    {
        acb_dft_precomp(w, v, c.pre, prec);
        /*_acb_dft_base(w, v, c.dv, c.z, c.dz, c.m, prec);*/
    }
    else
    {
        slong i, j;
        slong m = c.m, M = c.M, dv = c.dv, dz = c.dz;
        acb_srcptr z = c.z;
        acb_ptr t;
#if REORDER
        acb_ptr w2;
#endif

        t = _acb_vec_init(m * M);

        if (w == v)
        {
            _acb_vec_set(t, v, m * M);
            v = t;
        }

        /* m DFT of size M */
        for (i = 0; i < m; i++)
            acb_dft_step(w + i * M, v + i * dv, cyc + 1, num - 1, prec);

        /* twiddle if non trivial product */
        if (c.z != NULL)
        {
            acb_ptr wi;
            for (wi = w + M, i = 1; i < m; i++, wi += M)
                for (j = 1; j < M; j++)
                {
                    if (DFT_VERB)
                        flint_printf("z[%wu*%wu]",dz,i*j);
                    acb_mul(wi + j, wi + j, z + dz * i * j, prec);
                }
            if (DFT_VERB)
                flint_printf("\n");
        }

#if REORDER
        /* reorder w to avoid dv shifts in next DFT */
        w2 = flint_malloc(m * M * sizeof(acb_struct));
        for (j = 0; j < M; j++)
            for (i = 0; i < m; i++)
                w2[j + M * i] = w[i + m * j];
#endif

        /* M DFT of size m */
        for (j = 0; j < M; j++)
            acb_dft_precomp(t + m * j, w + j, c.pre, prec);

        /* reorder */
        for (i = 0; i < m; i++)
            for (j = 0; j < M; j++)
                acb_set(w + j + M * i, t + i + m * j);

        _acb_vec_clear(t, m * M);
    }
}
