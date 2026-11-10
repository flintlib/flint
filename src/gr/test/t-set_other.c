/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* generated using Claude Opus 5.5 */

#include "test_helpers.h"
#include "gr.h"

TEST_GR_FUNCTION_START(gr_set_other, state, count_success, count_domain, count_unable)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_struct ctx[3];
        gr_ptr x, y, z, w;
        int st01, st12, st02;
        slong j;
        int check_domain = 1;

        for (j = 0; j < 3; j++)
        {
            switch (n_randint(state, 16))
            {
                /* most of the time, test conversions between subrings of the
                   complex numbers */
                case 0: gr_ctx_init_fmpz(ctx + j); break;
                case 1: gr_ctx_init_fmpq(ctx + j); break;
                case 2: gr_ctx_init_fmpzi(ctx + j); break;
                case 3: gr_ctx_init_real_arb(ctx + j, 2 + n_randint(state, 200)); break;
                case 4: gr_ctx_init_complex_acb(ctx + j, 2 + n_randint(state, 200)); break;
                case 5: gr_ctx_init_real_ca(ctx + j); break;
                case 6: gr_ctx_init_complex_ca(ctx + j); break;
                case 7: gr_ctx_init_real_algebraic_ca(ctx + j); break;
                case 8: gr_ctx_init_complex_algebraic_ca(ctx + j); break;
                case 9: gr_ctx_init_complex_extended_ca(ctx + j); break;
                case 10: gr_ctx_init_real_qqbar(ctx + j); break;
                case 11: gr_ctx_init_complex_qqbar(ctx + j); break;

                default:
                         gr_ctx_init_random(ctx + j, state);
                         check_domain = 0;
            }
        }

        x = gr_heap_init(ctx + 0);
        y = gr_heap_init(ctx + 1);
        z = gr_heap_init(ctx + 2);
        w = gr_heap_init(ctx + 2);

        GR_MUST_SUCCEED(gr_randtest(x, state, ctx + 0));

        st01 = gr_set_other(y, x, ctx + 0, ctx + 1);
        st12 = gr_set_other(z, y, ctx + 1, ctx + 2);
        st02 = gr_set_other(w, x, ctx + 0, ctx + 2);

        if (((st01 | st12 | st02) == GR_SUCCESS && gr_equal(z, w, ctx + 2) == T_FALSE)
            || (check_domain && (st01 | st12) == GR_SUCCESS && st02 == GR_DOMAIN)
            || (check_domain && (st01 | st02) == GR_SUCCESS && st12 == GR_DOMAIN))
        {
            flint_printf("FAIL\n\n");
            flint_printf("%{gr_ctx}\n --> %{gr_ctx}\n --> %{gr_ctx}\n",
                            ctx, ctx + 1, ctx + 2);
            flint_printf("st12 = %d\n", st12);
            flint_printf("x = "); gr_println(x, ctx + 0);
            flint_printf("y = "); gr_println(y, ctx + 1);
            flint_printf("z = "); gr_println(z, ctx + 2);
            flint_printf("w = "); gr_println(w, ctx + 2);
            flint_abort();
        }

        count_success += (st12 == GR_SUCCESS);
        count_domain += ((st12 & GR_DOMAIN) != 0);
        count_unable += ((st12 & GR_UNABLE) != 0);

        gr_heap_clear(x, ctx + 0);
        gr_heap_clear(y, ctx + 1);
        gr_heap_clear(z, ctx + 2);
        gr_heap_clear(w, ctx + 2);

        for (j = 0; j < 3; j++)
            gr_ctx_clear(ctx + j);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
