/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "dlog.h"

/* group components up to bound and return cofactor */
void
dlog_n_factor_group(n_factor_t * fac, ulong bound)
{
    int i, j, k;
    ulong m[FLINT_MAX_FACTORS_IN_LIMB];
    ulong c = 1;
    i = 0;

    for (k = fac->num - 1; k >= 0; k--)
    {
        ulong qe = n_pow(fac->p[k], fac->exp[k]);
        if (qe > bound)
            c *= qe;
        else
        {
            /* try to insert somewhere in m */
            for (j = 0; j < i && (m[j] * qe > bound); j++);
            if (j == i)
                m[i++] = qe;
            else
                m[j] *= qe;
        }
    }

    for (j = 0; j < i; j++)
    {
        fac->p[j] = m[j];
        fac->exp[j] = DLOG_G_SMALL;
    }

    if (c > 1)
    {
        fac->p[i] = c;
        fac->exp[i] = DLOG_G_BIG;
        i++;
    }

    fac->num = i;
}
