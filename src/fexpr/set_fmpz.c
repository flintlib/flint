/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fexpr.h"

void
fexpr_set_fmpz(fexpr_t res, const fmpz_t c)
{
    if (!COEFF_IS_MPZ(*c))
    {
        fexpr_set_si(res, *c);
    }
    else
    {
        slong nlimbs;
        __mpz_struct * z = COEFF_TO_PTR(*c);

        nlimbs = FLINT_ABS(z->_mp_size);
        fexpr_fit_size(res, 1 + nlimbs);
        res->data[0] = ((z->_mp_size > 0) ? FEXPR_TYPE_BIG_INT_POS : FEXPR_TYPE_BIG_INT_NEG) | ((1 + nlimbs) << FEXPR_TYPE_BITS);

        flint_mpn_copyi(res->data + 1, z->_mp_d, nlimbs);
    }
}
