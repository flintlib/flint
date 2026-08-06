/*
    Copyright (C) 2026 Vincent Neiger, Kevin Tran

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "test_helpers.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_extrapolate_geometric, state)
{
    int i, result = 1;

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t P;
        nn_ptr ival, oval, val;
        nn_ptr ipts, opts;
        ulong modn, r, q;
        nmod_t mod;
        slong ilen, olen, istart, offset;

        ilen = (i < 10) ? i : n_randint(state, 500);
        olen = (i < 10) ? i : n_randint(state, 500);

        ival = FLINT_ARRAY_ALLOC(ilen, ulong);
        oval = FLINT_ARRAY_ALLOC(olen, ulong);
        val = FLINT_ARRAY_ALLOC(olen, ulong);
        ipts = FLINT_ARRAY_ALLOC(ilen, ulong);
        opts = FLINT_ARRAY_ALLOC(olen, ulong);

        /* forward: offset >= ilen */
        {
            istart = n_randint(state, 50);
            offset = ilen + n_randint(state, 10);

            /* choose modn large enough for r to exist:  */
            do 
            { 
                modn = n_randtest_prime(state, 1); 
            }
            while (modn <= 2*((ulong)istart+offset+olen) + 1);

            nmod_init(&mod, modn);
            nmod_poly_init(P, modn);
            r = n_primitive_root_prime(modn);
            q = nmod_mul(r, r, mod);

            _nmod_vec_randtest(ival, state, ilen, mod);

            ulong c = 1 + n_randint(state, modn-1);

            /* input points are c * q**k, k = istart...istart+ilen-1 */
            if (ilen > 0)
            {
                ipts[0] = nmod_pow_ui(q, istart, mod);
                ipts[0] = nmod_mul(ipts[0], c, mod);
                for (slong k = 1; k < ilen; k++)
                    ipts[k] = nmod_mul(ipts[k-1], q, mod);
            }

            /* output points are c * q**k, k = istart+offset...istart+offset+olen-1 */
            if (olen > 0)
            {
                opts[0] = nmod_pow_ui(q, istart+offset, mod);
                opts[0] = nmod_mul(opts[0], c, mod);
                for (slong k = 1; k < olen; k++)
                    opts[k] = nmod_mul(opts[k-1], q, mod);
            }

            nmod_poly_interpolate_nmod_vec(P, ipts, ival, ilen);
            nmod_poly_evaluate_nmod_vec(val, P, opts, olen);

            nmod_poly_extrapolate_geometric(oval, olen, ival, ilen, offset, r, mod);

            result = _nmod_vec_equal(oval, val, olen);

            if (!result)
            {
                flint_printf("FAIL (forward):\n");
                flint_printf("mod=%wu, ilen=%wd, olen=%wd, offset=%wd\n\n", modn, ilen, olen, offset);
                flint_printf("ival: "); _nmod_vec_print_pretty(ival, ilen, mod); flint_printf("\n\n");
                flint_printf("oval: "); _nmod_vec_print_pretty(oval, olen, mod); flint_printf("\n\n");
                flint_printf("correct: "); _nmod_vec_print_pretty(val, olen, mod); flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }

            nmod_poly_clear(P);
        }

        /* backward: offset + olen <= 0 */
        {
            /* here we build the points explicitly, so we also need istart+offset >= 0 */
            istart = olen + 10 + n_randint(state, 50);
            offset = - (olen + (slong)n_randint(state, 10));

            /* choose modn large enough for r to exist */
            do 
            { 
                modn = n_randtest_prime(state, 1); 
            }
            while (modn <= 2*((ulong)istart+ilen) + 1);

            nmod_init(&mod, modn);
            nmod_poly_init(P, modn);
            r = n_primitive_root_prime(modn);
            q = nmod_mul(r, r, mod);

            _nmod_vec_randtest(ival, state, ilen, mod);

            ulong c = 1 + n_randint(state, modn-1);

            /* input points are c * q**k, k = istart...istart+ilen-1 */
            if (ilen > 0)
            {
                ipts[0] = nmod_pow_ui(q, istart, mod);
                ipts[0] = nmod_mul(ipts[0], c, mod);
                for (slong k = 1; k < ilen; k++)
                    ipts[k] = nmod_mul(ipts[k-1], q, mod);
            }

            /* output points are c * q**k, k = istart+offset...istart+offset+olen-1 */
            if (olen > 0)
            {
                opts[0] = nmod_pow_ui(q, istart+offset, mod);
                opts[0] = nmod_mul(opts[0], c, mod);
                for (slong k = 1; k < olen; k++)
                    opts[k] = nmod_mul(opts[k-1], q, mod);
            }

            nmod_poly_interpolate_nmod_vec(P, ipts, ival, ilen);
            nmod_poly_evaluate_nmod_vec(val, P, opts, olen);

            nmod_poly_extrapolate_geometric(oval, olen, ival, ilen, offset, r, mod);

            result = _nmod_vec_equal(oval, val, olen);

            if (!result)
            {
                flint_printf("FAIL (backward):\n");
                flint_printf("mod=%wu, ilen=%wd, olen=%wd, ", modn, ilen, olen);
                flint_printf("istart=%wd, offset=%wd, q=%wu, c=%wu\n\n", istart, offset, q, c);
                if (ilen < 20 && olen < 20)
                {
                    flint_printf("ival: "); _nmod_vec_print_pretty(ival, ilen, mod); flint_printf("\n\n");
                    flint_printf("oval: "); _nmod_vec_print_pretty(oval, olen, mod); flint_printf("\n\n");
                    flint_printf("correct: "); _nmod_vec_print_pretty(val, olen, mod); flint_printf("\n\n");
                }
                fflush(stdout);
                flint_abort();
            }
            nmod_poly_clear(P);
        }

        flint_free(ival);
        flint_free(oval);
        flint_free(val);
        flint_free(ipts);
        flint_free(opts);
    }

    TEST_FUNCTION_END(state);
}
