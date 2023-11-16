/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"
#include "dlog.h"

TEST_FUNCTION_START(dlog, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        dlog_table_t table;
        dlog_bsgs_t bsgs;
        dlog_crt_t crt;
        dlog_precomp_t pre1, pre100;
        ulong p, a, k;
        nmod_t modp;

        if (iter < 10)
            p = n_nth_prime(iter + 2);
        else
            p = n_randprime(state, 15, 0);
        nmod_init(&modp, p);

        a = n_primitive_root_prime(p);

        dlog_table_init(table, a, p);
        dlog_bsgs_init(bsgs, a, p, p-1, dlog_bsgs_size(p, 1));
        dlog_crt_init(crt, a, p, p-1, 10);
        dlog_precomp_n_init(pre1, a, p, p-1, 1);
        dlog_precomp_n_init(pre100, a, p, p-1, 100);

        for (k = 1; k < 100 && k < p; k++)
        {
            ulong l0, l1, l2, l3, l4, l5;

            l1 = dlog_table(table, k);
            l2 = dlog_bsgs(bsgs, k);
            l3 = dlog_crt(crt, k);
            l4 = dlog_precomp(pre1, k);
            l5 = dlog_precomp(pre100, k);

            if (iter < 50 && k <= 7)
                l0 = dlog_once(k, a, modp, p-1);
            else
                l0 = l1;

            if (l0 != l1 || l1 != l2 || l1 != l3 || l1 != l4 || l1 != l5)
            {
                flint_printf("\n\nFAIL: log(%wu,%wu) mod %wu\n\n",k,a,p);
                flint_printf("once: %wu\ntable: %wu\nbsgs: %wu\ncrt: %wu\nprecomp1: %wu\nprecomp100: %wu\n\n",
                        l0, l1, l2, l3, l4, l5);
                flint_abort();
            }
        }
        dlog_table_clear(table);
        dlog_bsgs_clear(bsgs);
        dlog_crt_clear(crt);
        dlog_precomp_clear(pre1);
        dlog_precomp_clear(pre100);
    }

    TEST_FUNCTION_END(state);
}
