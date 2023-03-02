/*
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

int TEMPLATE(T, is_square)(const TEMPLATE(T, t) op,
                                      const TEMPLATE(T, ctx_t) ctx)
{
    fmpz_t ord;
    TEMPLATE(T, t) pow;
    int is_square = 0;

    if (TEMPLATE(T, is_zero)(op, ctx) || TEMPLATE(T, is_one)(op, ctx) ||
		              fmpz_cmp_ui(TEMPLATE(T, ctx_prime)(ctx), 2) == 0)
    {
        return 1;
    }

    fmpz_init(ord);
    TEMPLATE(T, init)(pow, ctx);

    TEMPLATE(T, ctx_order)(ord, ctx);
    fmpz_sub_ui(ord, ord, 1);
    fmpz_tdiv_q_2exp(ord, ord, 1);

    TEMPLATE(T, pow)(pow, op, ord, ctx);

    is_square = TEMPLATE(T, is_one)(pow, ctx);

    fmpz_clear(ord);
    TEMPLATE(T, clear)(pow, ctx);
    
    return is_square;
}

#endif
