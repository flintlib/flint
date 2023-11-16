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

#if FLINT_BITS == 64
#define LIM UWORD(1000000000000)
#else
#define LIM UWORD(1000000000)
#endif

TEST_FUNCTION_START(dlog_modpe, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        ulong p, e, pe, a;

        if (iter == 0)
        {
            p = 2;
            pe = 8;
            a = 5;
            e = 3;
        }
        else
        {
            p = pe = n_randprime(state, 10, 0);
            a = (p == 40487) ? 10 : n_primitive_root_prime(p);
            e = 1;
        }

        for (; pe < LIM; pe *= p, e++)
        {
            ulong k, phi;
            nmod_t mod;

            dlog_modpe_t modpe;

            nmod_init(&mod, pe);
            phi = (p == 2) ? pe / 4 : pe - pe / p;

            dlog_modpe_init(modpe, a, p, e, pe, 10);

            for (k = 0; k < 100 && k < p; k++)
            {
                ulong l, b, x;

                l = n_randint(state, phi);
                b = nmod_pow_ui(a, l, mod);

                if ((x = dlog_modpe(modpe, b)) != l)
                {
                    flint_printf("FAIL modpe: %wu^%wu = %wu [%wu^%wu]\n\n",
                            a, l, b, p, e);
                    flint_printf("modpe returned %wu\n\n", x);
                    flint_abort();
                }
            }

            dlog_modpe_clear(modpe);

            /* multiplication can overflow on 32-bit */
            if ((double) pe * p > LIM)
                break;
        }
    }

    TEST_FUNCTION_END(state);
}
