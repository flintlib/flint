/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Fredrik Johansson

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

void
_fmpz_factor_fit_length(fmpz_factor_t factor, slong len)
{
    if (len > factor->alloc)
    {
        if (len < 2 * factor->alloc)
            len = 2 * factor->alloc;

        factor->p = (fmpz *) flint_realloc(factor->p, len * sizeof(fmpz));
        factor->exp = flint_realloc(factor->exp, len * sizeof(slong));

        if (len > factor->alloc)
        {
            flint_mpn_zero((mp_ptr)(factor->p + factor->alloc), len-factor->alloc);
            flint_mpn_zero((mp_ptr)(factor->exp + factor->alloc), len-factor->alloc);
        }

        factor->alloc = len;
    }
}
