/*
    Copyright (C) 2018 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

int TEMPLATE(T, multiplicative_order)(fmpz * ord, const TEMPLATE(T, t) op,
                                      const TEMPLATE(T, ctx_t) ctx)
{
    fmpz_t tmp;
    fmpz_factor_t ord_fact;
    TEMPLATE(T, t) one;
    slong i, j;
    int is_primitive = 1;

    if (ord == NULL)
        ord = tmp;

    fmpz_init(tmp);

    if (TEMPLATE(T, is_zero)(op, ctx))
    {
        fmpz_zero(ord);
        is_primitive = 0;
    }
    else
    {
        fmpz_factor_init(ord_fact);
        TEMPLATE(T, init)(one, ctx);

        TEMPLATE(T, ctx_order)(ord, ctx);
        fmpz_sub_ui(ord, ord, 1);
        fmpz_factor(ord_fact, ord);

        for (i = 0; i < ord_fact->num; i++)
        {
            for (j = ord_fact->exp[i]; j > 0; j--)
            {
                fmpz_cdiv_q(ord, ord, ord_fact->p + i);
                TEMPLATE(T, pow)(one, op, ord, ctx);
                if (!TEMPLATE(T, is_one)(one, ctx))
                    break;
                is_primitive = -1;
            }
            if (j > 0)
                fmpz_mul(ord, ord, ord_fact->p + i);
        }

        fmpz_factor_clear(ord_fact);
        TEMPLATE(T, clear)(one, ctx);
    }

    fmpz_clear(tmp);

    return is_primitive;
}

#endif
