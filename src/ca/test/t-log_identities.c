/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ca.h"

/* atan(x) = -i/2 log((1+ix)/(1-ix)) */
/* valid for -inf < x < inf */
void
simple_ca_atan(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    ca_t t, u, v;

    ca_init(t, ctx);
    ca_init(u, ctx);
    ca_init(v, ctx);

    ca_i(t, ctx);
    ca_mul(u, x, t, ctx);
    /* v = 1 + i x */
    ca_add_ui(v, u, 1, ctx);
    /* res = 1 - i x */
    ca_sub_ui(res, u, 1, ctx);
    ca_neg(res, res, ctx);

    ca_div(res, v, res, ctx);
    ca_log(res, res, ctx);

    ca_mul(res, res, t, ctx);

    ca_div_ui(res, res, 2, ctx);
    ca_neg(res, res, ctx);

    ca_clear(t, ctx);
    ca_clear(u, ctx);
    ca_clear(v, ctx);
}

void
simple_ca_atan_p_q(ca_t res, ulong p, ulong q, ca_ctx_t ctx)
{
    ca_set_ui(res, p, ctx);
    ca_div_ui(res, res, q, ctx);
    simple_ca_atan(res, res, ctx);
}

/* valid for -1 < x < 1 */
void
simple_ca_atanh(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    ca_t t, u;

    ca_init(t, ctx);
    ca_init(u, ctx);

    ca_add_ui(t, x, 1, ctx);
    ca_sub_ui(u, x, 1, ctx);
    ca_neg(u, u, ctx);
    ca_div(res, t, u, ctx);
    ca_log(res, res, ctx);
    ca_div_ui(res, res, 2, ctx);

    ca_clear(t, ctx);
    ca_clear(u, ctx);
}

void
simple_ca_atanh_p_q(ca_t res, ulong p, ulong q, ca_ctx_t ctx)
{
    ca_set_ui(res, p, ctx);
    ca_div_ui(res, res, q, ctx);
    simple_ca_atanh(res, res, ctx);
}

#define NUM_FORMULAS 8

slong machin_formulas[NUM_FORMULAS][4][2] = {
    {{1, 1}, {0, 0}, {0, 0}, {0, 0}},

    {{1, 2}, {1, 3}, {0, 0}, {0, 0}},
    {{2, 2}, {-1, 7}, {0, 0}, {0, 0}},
    {{2, 3}, {1, 7}, {0, 0}, {0, 0}},
    {{4, 5}, {-1, 239}, {0, 0}, {0, 0}},

    {{1, 2}, {1, 5}, {1, 8}, {0, 0}},
    {{1, 3}, {1, 4}, {1, 7}, {1, 13}},
    {{12, 49}, {32, 57}, {-5, 239}, {12, 110443}},
};

#define NUM_FORMULAS2 7

slong hyperbolic_logs[NUM_FORMULAS2] = {2, 3, 5, 2, 3, 5, 7};

slong hyperbolic_machin_formulas[NUM_FORMULAS2][4][2] = {
    {{14, 31}, {10, 49}, {6, 161}, {0, 0}},
    {{22, 31}, {16, 49}, {10, 161}, {0, 0}},
    {{32, 31}, {24, 49}, {14, 161}, {0, 0}},

    {{144, 251}, {54, 449}, {-38, 4801}, {62, 8749}},
    {{228, 251}, {86, 449}, {-60, 4801}, {98, 8749}},
    {{334, 251}, {126, 449}, {-88, 4801}, {144, 8749}},
    {{404, 251}, {152, 449}, {-106, 4801}, {174, 8749}},
};

TEST_FUNCTION_START(ca_log_identities, state)
{
    slong iter;

    for (iter = 0; iter < 1; iter++)
    {
        ca_ctx_t ctx;
        ca_t x, y, pi4;
        slong i, j, c, q;

        ca_ctx_init(ctx);
        ca_init(x, ctx);
        ca_init(y, ctx);
        ca_init(pi4, ctx);

        ca_pi(pi4, ctx);
        ca_div_ui(pi4, pi4, 4, ctx);

        for (i = 0; i < NUM_FORMULAS; i++)
        {
            ca_zero(x, ctx);
            for (j = 0; j < 4; j++)
            {
                c = machin_formulas[i][j][0];
                q = machin_formulas[i][j][1];

                if (c != 0)
                {
                    simple_ca_atan_p_q(y, 1, q, ctx);
                    ca_mul_si(y, y, c, ctx);
                    ca_add(x, x, y, ctx);
                }
            }

            ca_sub(x, x, pi4, ctx);

            if (!(CA_IS_QQ(x, ctx) && fmpq_is_zero(CA_FMPQ(x))))
            {
                flint_printf("FAIL\n\n");
                for (j = 0; j < 4; j++)
                {
                    c = machin_formulas[i][j][0];
                    q = machin_formulas[i][j][1];
                    flint_printf("%wd %wd   ", c, q);
                }
                flint_printf("\n\n");
                flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
                flint_abort();
            }
        }

        for (i = 0; i < NUM_FORMULAS2; i++)
        {
            ca_zero(x, ctx);
            for (j = 0; j < 4; j++)
            {
                c = hyperbolic_machin_formulas[i][j][0];
                q = hyperbolic_machin_formulas[i][j][1];

                if (c != 0)
                {
                    simple_ca_atanh_p_q(y, 1, q, ctx);
                    ca_mul_si(y, y, c, ctx);
                    ca_add(x, x, y, ctx);
                }
            }

            ca_set_ui(y, hyperbolic_logs[i], ctx);
            ca_log(y, y, ctx);
            ca_sub(x, x, y, ctx);

            if (!(CA_IS_QQ(x, ctx) && fmpq_is_zero(CA_FMPQ(x))))
            {
                flint_printf("FAIL\n\n");
                for (j = 0; j < 4; j++)
                {
                    c = hyperbolic_machin_formulas[i][j][0];
                    q = hyperbolic_machin_formulas[i][j][1];
                    flint_printf("%wd %wd   ", c, q);
                }
                flint_printf("\n\n");
                flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
                flint_abort();
            }
        }

        ca_clear(x, ctx);
        ca_clear(y, ctx);
        ca_clear(pi4, ctx);
        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
