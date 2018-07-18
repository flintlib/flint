/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

void mpoly_zipinfo_init(mpoly_zipinfo_t zinfo, slong nvars)
{
    zinfo->nvars = nvars;
    zinfo->Adegs = (slong *) flint_malloc(nvars*sizeof(slong));
    zinfo->Bdegs = (slong *) flint_malloc(nvars*sizeof(slong));
    zinfo->perm  = (slong *) flint_malloc(nvars*sizeof(slong));
}

void mpoly_zipinfo_clear(mpoly_zipinfo_t zinfo)
{
    flint_free(zinfo->Adegs);
    flint_free(zinfo->Bdegs);
    flint_free(zinfo->perm);
}
