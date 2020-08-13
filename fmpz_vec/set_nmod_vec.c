/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "nmod_poly.h"

void
_fmpz_vec_set_nmod_vec(fmpz * res, mp_srcptr poly, slong len, nmod_t mod)
{
    slong i;
        
    for (i = 0; i < len; i++)
       fmpz_set_ui_smod(res + i, poly[i], mod.n);
}
