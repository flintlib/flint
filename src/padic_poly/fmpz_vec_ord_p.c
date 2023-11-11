/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "padic_poly.h"

/*
    Returns the minimum $p$-adic valuation of \code{(vec, len)}, assuming this fits into a \code{signed long}.

    If \code{len} is zero, returns $0$.
 */
slong _fmpz_vec_ord_p(const fmpz * vec, slong len, const fmpz_t p)
{
    if (len == 0)
    {
        return 0;
    }
    else
    {
        fmpz_t t;
        slong i, min = WORD_MAX, v;

        fmpz_init(t);
        for (i = 0; (min > 0) && (i < len); i++)
        {
            if (!fmpz_is_zero(vec + i))
            {
                v   = fmpz_remove(t, vec + i, p);
                min = FLINT_MIN(min, v);
            }
        }
        fmpz_clear(t);
        return (min < WORD_MAX) ? min : 0;
    }
}
