/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "fmpz_vec.h"

slong
_fmpz_vec_get_d_vec_2exp(double *appv, const fmpz * vec, slong len)
{
    slong *exp, i, maxexp = 0L;
    exp = (slong *) flint_malloc(len * sizeof(slong));

    for (i = 0; i < len; i++)
    {
        appv[i] = fmpz_get_d_2exp(&exp[i], vec + i);
        if (exp[i] > maxexp)
            maxexp = exp[i];
    }

    for (i = 0; i < len; i++)
        appv[i] = ldexp(appv[i], exp[i] - maxexp);

    flint_free(exp);
    return maxexp;
}
