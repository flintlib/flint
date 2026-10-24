/*
    Copyright (C) 2026 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "nmod_mpoly_factor.h"

TEST_FUNCTION_START(nmod_mpoly_factor_irred_smprime, state)
{
    const char * vars[] = {"x", "y", "z"};
    const char * inputs[] = {
        "x^3*y^3+x^3*y*z+6*x^2*y^2+3*x^2*y+11*x*y+2*x+6",
        "x^4*y^3+x^4*y*z+x^3*y^2+x^3*y+x^3+2*x^2*y^2+x^2*y+2*x^2+3*x+2"
    };
    const ulong primes[] = {7, 11};
    const ulong seeds[] = {74, 146};
    slong i;
    int algo, success;
    flint_rand_t retry_state;

    /*
        Along y = z = t the inputs become
          (t*x+1)*(t*x+2)*((t+1)*x+3) and
          (t^2*x^2+x+1)*((t+1)*x^2+x+2).
        Both inputs are primitive and linear in z, hence irreducible.
        Their leading coefficient is y*(y^2+z), but the specialized
        factors can assign y too many times. The fixed seeds exercise
        the subsequent and initial failed divisions, respectively.
    */
    flint_rand_init(retry_state);
    for (i = 0; i < 2; i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, lc;
        nmod_mpoly_factor_t lcf;
        nmod_mpolyv_t fac;
        nmod_mpoly_ctx_init(ctx, 3, ORD_LEX, primes[i]);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(lc, ctx);
        nmod_mpoly_factor_init(lcf, ctx);
        nmod_mpolyv_init(fac, ctx);
        FLINT_TEST(nmod_mpoly_set_str_pretty(
                                          a, inputs[i], vars, ctx) == 0);

        /* Keep y first in the leading coefficient factorization. */
        nmod_mpoly_factor_fit_length(lcf, 2, ctx);
        nmod_mpoly_gen(lcf->poly + 0, 1, ctx);
        nmod_mpoly_mul(lc, lcf->poly + 0, lcf->poly + 0, ctx);
        nmod_mpoly_gen(lcf->poly + 1, 2, ctx);
        nmod_mpoly_add(lcf->poly + 1, lcf->poly + 1, lc, ctx);
        nmod_mpoly_mul(lc, lcf->poly + 0, lcf->poly + 1, ctx);
        fmpz_one(lcf->exp + 0);
        fmpz_one(lcf->exp + 1);
        lcf->num = 2;

        for (algo = 0; algo < 2; algo++)
        {
            flint_rand_set_seed(retry_state, seeds[i], seeds[i] + 1);
            success = algo == 0 ?
                nmod_mpoly_factor_irred_smprime_wang(
                                        fac, a, lcf, lc, ctx, retry_state) :
                nmod_mpoly_factor_irred_smprime_zippel(
                                        fac, a, lcf, lc, ctx, retry_state);
            /* The evaluation attempts may be exhausted. */
            FLINT_TEST(success >= 0);
            if (success)
            {
                FLINT_TEST(fac->length == 1);
                FLINT_TEST(nmod_mpoly_equal(fac->coeffs + 0, a, ctx));
            }
        }

        nmod_mpolyv_clear(fac, ctx);
        nmod_mpoly_factor_clear(lcf, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(lc, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }
    flint_rand_clear(retry_state);

    TEST_FUNCTION_END(state);
}
