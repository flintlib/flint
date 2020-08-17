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

int TEMPLATE(T, sqrt)(TEMPLATE(T, t) rop, const TEMPLATE(T, t) op,
                                                  const TEMPLATE(T, ctx_t) ctx)
{
    int res = 1;

    if (TEMPLATE(T, is_zero)(op, ctx) || TEMPLATE(T, is_one)(op, ctx))
        TEMPLATE(T, set)(rop, op, ctx);
    else if (fmpz_cmp_ui(TEMPLATE(T, ctx_prime)(ctx), 2) == 0)
        TEMPLATE(T, pth_root)(rop, op, ctx);
    else
    {
        TEMPLATE(T, t) z, c, t, b, temp;
        fmpz_t ord, Q, Q2;
        flint_bitcnt_t S;
        slong M, i, j;

        TEMPLATE(T, init)(z, ctx);
        TEMPLATE(T, init)(c, ctx);
        TEMPLATE(T, init)(t, ctx);
        TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(temp, ctx);

        fmpz_init(ord);
        fmpz_init(Q);
        fmpz_init(Q2);

        if (ctx->is_conway)
            TEMPLATE(T, gen)(z, ctx);
        else
        {
            flint_rand_t state;

            flint_randinit(state);

            while (TEMPLATE(T, is_square)(z, ctx))
                TEMPLATE(T, rand)(z, state, ctx);

            flint_randclear(state);
        }	

        TEMPLATE(T, ctx_order)(ord, ctx);
        fmpz_sub_ui(ord, ord, 1);
        S = fmpz_val2(ord);
        fmpz_tdiv_q_2exp(Q, ord, S);
        fmpz_add_ui(Q2, Q, 1);
        fmpz_tdiv_q_2exp(Q2, Q2, 1);

        M = S;
        TEMPLATE(T, pow)(c, z, Q, ctx);
        TEMPLATE(T, pow)(t, op, Q, ctx);
        TEMPLATE(T, pow)(rop, op, Q2, ctx);

        while (1)
        {
            if (TEMPLATE(T, is_zero)(t, ctx))
            {
                TEMPLATE(T, zero)(rop, ctx);
                break;
            }

	    if (TEMPLATE(T, is_one)(t, ctx))
                break;

            i = 1;
	    TEMPLATE(T, sqr)(temp, t, ctx);

            while (i < M && !TEMPLATE(T, is_one)(temp, ctx))
	    {
                TEMPLATE(T, sqr)(temp, temp, ctx);
                i++;
            }

            if (i == M)
            {
                res = 0;
		break;
            }

            TEMPLATE(T, set)(b, c, ctx);

            for (j = 0; j < M - i - 1; j++)
                TEMPLATE(T, sqr)(b, b, ctx);

            M = i;
            TEMPLATE(T, sqr)(c, b, ctx);
            TEMPLATE(T, mul)(t, t, c, ctx);
            TEMPLATE(T, mul)(rop, rop, b, ctx);
        }

        fmpz_clear(Q2);
        fmpz_clear(Q);
        fmpz_clear(ord);

        TEMPLATE(T, clear)(temp, ctx);
        TEMPLATE(T, clear)(b, ctx);
        TEMPLATE(T, clear)(t, ctx);
        TEMPLATE(T, clear)(c, ctx);
        TEMPLATE(T, clear)(z, ctx);
    }

    return res;
}

#endif
