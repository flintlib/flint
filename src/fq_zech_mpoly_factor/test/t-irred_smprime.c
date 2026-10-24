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
#include "nmod_mpoly.h"
#include "fq_zech_mpoly_factor.h"

TEST_FUNCTION_START(fq_zech_mpoly_factor_irred_smprime, state)
{
    const char * vars[] = {"x", "y", "z"};
    const char * inputs[] = {
        "x^3*y^3+x^3*y*z+6*x^2*y^2+3*x^2*y+11*x*y+2*x+6",
        "x^4*y^3+x^4*y*z+x^3*y^2+x^3*y+x^3+2*x^2*y^2+x^2*y+2*x^2+3*x+2"
    };
    const ulong primes[] = {7, 3};
    const slong degrees[] = {2, 3};
    const ulong seeds[] = {655, 110};
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
        fq_zech_mpoly_ctx_t ctx;
        fq_zech_mpoly_t a, lc;
        fq_zech_mpoly_factor_t lcf;
        fq_zech_mpolyv_t fac;
        nmod_mpoly_ctx_t nctx;
        nmod_mpoly_t t;

        fq_zech_mpoly_ctx_init_deg(ctx, 3, ORD_LEX, primes[i], degrees[i]);
        nmod_mpoly_ctx_init(nctx, 3, ORD_LEX, primes[i]);
        nmod_mpoly_init(t, nctx);
        fq_zech_mpoly_init(a, ctx);
        fq_zech_mpoly_init(lc, ctx);
        fq_zech_mpoly_factor_init(lcf, ctx);
        fq_zech_mpolyv_init(fac, ctx);
        FLINT_TEST(nmod_mpoly_set_str_pretty(t, inputs[i], vars, nctx) == 0);
        _fq_zech_mpoly_set_nmod_mpoly(a, ctx, t, nctx);

        /* Keep y first in the leading coefficient factorization. */
        fq_zech_mpoly_factor_fit_length(lcf, 2, ctx);
        fq_zech_mpoly_gen(lcf->poly + 0, 1, ctx);
        fq_zech_mpoly_mul(lc, lcf->poly + 0, lcf->poly + 0, ctx);
        fq_zech_mpoly_gen(lcf->poly + 1, 2, ctx);
        fq_zech_mpoly_add(lcf->poly + 1, lcf->poly + 1, lc, ctx);
        fq_zech_mpoly_mul(lc, lcf->poly + 0, lcf->poly + 1, ctx);
        fmpz_one(lcf->exp + 0);
        fmpz_one(lcf->exp + 1);
        lcf->num = 2;

        for (algo = 0; algo < 2; algo++)
        {
            flint_rand_set_seed(retry_state, seeds[i], seeds[i] + 1);
            success = algo == 0 ?
                fq_zech_mpoly_factor_irred_smprime_wang(
                                        fac, a, lcf, lc, ctx, retry_state) :
                fq_zech_mpoly_factor_irred_smprime_zippel(
                                        fac, a, lcf, lc, ctx, retry_state);
            /* The evaluation attempts may be exhausted. */
            FLINT_TEST(success >= 0);
            if (success)
            {
                FLINT_TEST(fac->length == 1);
                FLINT_TEST(fq_zech_mpoly_equal(fac->coeffs + 0, a, ctx));
            }
        }

        fq_zech_mpolyv_clear(fac, ctx);
        fq_zech_mpoly_factor_clear(lcf, ctx);
        fq_zech_mpoly_clear(a, ctx);
        fq_zech_mpoly_clear(lc, ctx);
        fq_zech_mpoly_ctx_clear(ctx);
        nmod_mpoly_clear(t, nctx);
        nmod_mpoly_ctx_clear(nctx);
    }
    flint_rand_clear(retry_state);

    TEST_FUNCTION_END(state);
}
