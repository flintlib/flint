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
_fmpz_factor_set_length(fmpz_factor_t factor, slong newlen)
{
    if (factor->num > newlen)
    {
        slong i;
        for (i = newlen; i < factor->num; i++)
            _fmpz_demote(factor->p + i); 
    }
    factor->num = newlen;
}
