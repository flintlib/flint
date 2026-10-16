/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mp_real.h"

/* Tiny arguments, down to a quarter of the safe exponent range: the
   remainder bounds of the first-order values (sin x = x, cos x = 1,
   tan x = x, atan x = x, exp x = 1, and the same for the functions of
   pi x and for balls around zero) must not pad the outputs down to the
   scale of the remainder, which would take up to B^(2^55) limbs. */

TEST_FUNCTION_START(mp_real_tiny_arg, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        mp_real_t x, y, z;
        slong prec = 2 + n_randint(state, 500), e, maxsize;
        ulong d = n_randtest_not_zero(state);
        int f;

        mp_real_init(x);
        mp_real_init(y);
        mp_real_init(z);

        e = -(slong) n_randint(state, FLINT_BITS * (MP_REAL_EXP_MAX / 4));
        _mp_real_set_mpn_2exp(x, &d, 1, e);
        if (n_randint(state, 2))
            mp_real_neg(x, x);
        if (n_randint(state, 2))
            mp_real_add_error_2exp_si(x, mp_real_abs_bound_lt_2exp_si(x) - prec - n_randint(state, 100));

        maxsize = mp_real_prec_bits(prec) + 8;

        for (f = 0; f < 8; f++)
        {
            switch (f)
            {
                case 0: mp_real_sin_cos_bits(y, z, x, prec); break;
                case 1: mp_real_tan_bits(y, x, prec); mp_real_zero(z); break;
                case 2: mp_real_atan_bits(y, x, prec); mp_real_zero(z); break;
                case 3: mp_real_exp_bits(y, x, prec); mp_real_zero(z); break;
                case 4: mp_real_sin_cos_pi_bits(y, z, x, prec); break;
                case 5: mp_real_tan_pi_bits(y, x, prec); mp_real_zero(z); break;
                case 6:
                    mp_real_zero(z);
                    mp_real_add_error_2exp_si(z, e);
                    mp_real_exp_bits(y, z, prec);
                    break;
                default:
                    mp_real_zero(z);
                    mp_real_add_error_2exp_si(z, e);
                    mp_real_sin_cos_bits(z, y, z, prec);
                    break;
            }

            if (y->size > maxsize || z->size > maxsize)
            {
                flint_printf("FAIL: f = %d, e = %wd, prec = %wd, sizes %wd %wd\n",
                    f, e, prec, y->size, z->size);
                flint_abort();
            }
        }

        mp_real_clear(x);
        mp_real_clear(y);
        mp_real_clear(z);
    }

    TEST_FUNCTION_END(state);
}
