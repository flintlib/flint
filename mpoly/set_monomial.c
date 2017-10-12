/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "mpoly.h"


void mpoly_set_monomial(ulong * poly_exps, const ulong * user_exps,
                                   slong bits, slong nfields, int deg, int rev)
{
    slong i = 0;
    ulong * tmp_exps, degree = 0;
    TMP_INIT;

    TMP_START;
    tmp_exps = (ulong *) TMP_ALLOC(nfields*sizeof(ulong));

    if (deg)
    {
        for (i = 0; i < nfields - 1; i++)
            degree += user_exps[i];
        tmp_exps[0] = degree;
    }

    if (rev)
    {
        for (i = deg; i < nfields; i++)
            tmp_exps[i] = user_exps[nfields - i - 1];
    } else
    {
        for (i = deg; i < nfields; i++)
            tmp_exps[i] = user_exps[i - deg];
    }

    mpoly_pack_vec(poly_exps, tmp_exps, bits, nfields, 1);

    TMP_END;
}
