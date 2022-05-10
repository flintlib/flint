/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2009, 2010 Andy Novocin
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_lll.h"

int
fmpz_lll_is_reduced(const fmpz_mat_t B, const fmpz_lll_t fl, flint_bitcnt_t prec)
{
    if (fmpz_lll_is_reduced_d(B, fl))
        return 1;

    if (fmpz_lll_is_reduced_mpfr(B, fl, prec))
        return 1;

    return (fl->rt == Z_BASIS) ?
                fmpz_mat_is_reduced(B, fl->delta, fl->eta) :
                fmpz_mat_is_reduced_gram(B, fl->delta, fl->eta);
}
