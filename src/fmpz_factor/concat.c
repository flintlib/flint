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
_fmpz_factor_concat(fmpz_factor_t factor1, fmpz_factor_t factor2, ulong exp)
{
    slong i;

    _fmpz_factor_fit_length(factor1, factor1->num + factor2->num);
    
    for (i = 0; i < factor2->num; i++)
    {
       fmpz_set(factor1->p + factor1->num + i, factor2->p + i);
       factor1->exp[factor1->num + i] = factor2->exp[i]*exp;
    }

    factor1->num += factor2->num;
}
