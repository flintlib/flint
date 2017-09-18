/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "mpoly.h"


void mpoly_get_monomial(ulong * user_exps, const ulong * poly_exps,
                                   slong bits, slong nfields, int deg, int rev)
{
    slong i;
    ulong * tmp_exps;
    TMP_INIT;

    TMP_START;
    tmp_exps = (ulong *) TMP_ALLOC(nfields*sizeof(ulong));
    mpoly_unpack_vec(tmp_exps, poly_exps, bits, nfields, 1);

    if (rev)
    {
        for (i = nfields - 1; i >= deg; i--)
            user_exps[nfields - i - 1] = tmp_exps[i];
    } else
    {
        for (i = deg; i < nfields; i++)
            user_exps[i - deg] = tmp_exps[i];
    }

    TMP_END;
}
