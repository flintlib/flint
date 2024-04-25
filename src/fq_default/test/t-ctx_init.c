/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* NOTE: Here we check that we get the contexts we expect. */

#include "test_helpers.h"
#include "fq_default.h"

static const char str_fq_zech[] = "fq_zech";
static const char str_fq_nmod[] = "fq_nmod";
static const char str_fq[] = "fq";
static const char str_nmod[] = "nmod";
static const char str_fmpz_mod[] = "fmpz_mod";

static const char * get_str(int type)
{
    switch (type)
    {
        case FQ_DEFAULT_FQ_ZECH: return str_fq_zech;
        case FQ_DEFAULT_FQ_NMOD: return str_fq_nmod;
        case FQ_DEFAULT_FQ: return str_fq;
        case FQ_DEFAULT_NMOD: return str_nmod;
        case FQ_DEFAULT_FMPZ_MOD: return str_fmpz_mod;
        default: FLINT_UNREACHABLE;
    }
}

static slong rand_deg_p_fq_zech(flint_rand_t state, fmpz_t p)
{
    /* d > 1 and d * bits(p) <= 16 */
    /* 2 <= bits(p) <= 16 / d */
    slong d = 2 + n_randint(state, 7); /* in {2, ..., 8} */
    flint_bitcnt_t bits = 2 + n_randint(state, 16 / d - 1);
    fmpz_randprime(p, state, bits, 0);
    return d;
}

static slong rand_deg_p_fq_nmod(flint_rand_t state, fmpz_t p)
{
    /* d > 1 and fits_ui(p) and d * bits(p) > 16 */
    slong d = 2 + n_randint(state, 20);
    flint_bitcnt_t bits = FLINT_MAX((16 + d - 1) / d + 1, 2) + n_randint(state, FLINT_BITS - FLINT_MAX((16 + d - 1) / d + 1, 2) + 1);
    fmpz_randprime(p, state, bits, 0);
    return d;
}

static slong rand_deg_p_nmod(flint_rand_t state, fmpz_t p)
{
    /* d == 1 and fits_ui(p) */
    slong d = 1;
    flint_bitcnt_t bits = 2 + n_randint(state, FLINT_BITS - 1);
    fmpz_randprime(p, state, bits, 0);
    return d;
}

static slong rand_deg_p_fmpz_mod(flint_rand_t state, fmpz_t p)
{
    /* d == 1 and not fits_ui(p) */
    slong d = 1;
    flint_bitcnt_t bits = FLINT_BITS + 1 + n_randint(state, 64);
    fmpz_randprime(p, state, bits, 0);
    return d;
}

static slong rand_deg_p_fq(flint_rand_t state, fmpz_t p)
{
    /* d > 1 and not fits_ui(p) */
    slong d = 2 + n_randint(state, 20);
    flint_bitcnt_t bits = FLINT_BITS + 1 + n_randint(state, 64);
    fmpz_randprime(p, state, bits, 0);
    return d;
}

#define TYPE_MIN FQ_DEFAULT_FQ_ZECH
#define TYPE_MAX FQ_DEFAULT_FMPZ_MOD

TEST_FUNCTION_START(fq_default_ctx_init, state)
{
    slong ix;
    int result;
    fq_default_ctx_t ctx;

    /* fq_default_ctx_init_type */
    for (ix = 0; ix < 100 * flint_test_multiplier(); ix++)
    {
        int type;
        fmpz_t p;
        slong d;
        char var[] = "x";

        type = TYPE_MIN + n_randint(state, TYPE_MAX + TYPE_MIN - 1);
        fmpz_init(p);

        switch (type)
        {
            case FQ_DEFAULT_FQ_ZECH: d = rand_deg_p_fq_zech(state, p); break;
            case FQ_DEFAULT_FQ_NMOD: d = rand_deg_p_fq_nmod(state, p); break;
            case FQ_DEFAULT_FQ: d = rand_deg_p_fq(state, p); break;
            case FQ_DEFAULT_NMOD: d = rand_deg_p_nmod(state, p); break;
            case FQ_DEFAULT_FMPZ_MOD: d = rand_deg_p_fmpz_mod(state, p); break;
            default: FLINT_UNREACHABLE;
        }

        fq_default_ctx_init_type(ctx, p, d, var, 0);

        result = (fq_default_ctx_type(ctx) == type);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "fq_default_ctx_init_type:\n"
                    "ix = %wd\n"
                    "p = %{fmpz}\n"
                    "d = %{slong}\n"
                    "Expected type: %s\n"
                    "Got type:      %s\n",
                    ix, p, d, get_str(type), get_str(fq_default_ctx_type(ctx)));

        fq_default_ctx_clear(ctx);
        fmpz_clear(p);
    }

    /* fq_default_ctx_init_modulus */
    for (ix = 0; ix < 100 * flint_test_multiplier(); ix++)
    {
        int type;
        fmpz_mod_poly_t modulus;
        fmpz_mod_ctx_t mod_ctx;
        fmpz_t p;
        slong d;
        char var[] = "x";

        type = TYPE_MIN + n_randint(state, TYPE_MAX + TYPE_MIN - 1);
        fmpz_init(p);

        /* Get degree and prime */
        switch (type)
        {
            case FQ_DEFAULT_FQ_ZECH: d = rand_deg_p_fq_zech(state, p); break;
            case FQ_DEFAULT_FQ_NMOD: d = rand_deg_p_fq_nmod(state, p); break;
            case FQ_DEFAULT_FQ: d = rand_deg_p_fq(state, p); break;
            case FQ_DEFAULT_NMOD: d = rand_deg_p_nmod(state, p); break;
            case FQ_DEFAULT_FMPZ_MOD: d = rand_deg_p_fmpz_mod(state, p); break;
            default: FLINT_UNREACHABLE;
        }

        fmpz_mod_ctx_init(mod_ctx, p);
        fmpz_mod_poly_init(modulus, mod_ctx);

        if (type != FQ_DEFAULT_FQ_ZECH)
            fmpz_mod_poly_randtest_monic_irreducible(modulus, state, d + 1, mod_ctx);
        else
            fmpz_mod_poly_randtest_monic_primitive(modulus, state, d + 1, mod_ctx);

        fq_default_ctx_init_modulus(ctx, modulus, mod_ctx, var);

        result = (fq_default_ctx_type(ctx) == type);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "fq_default_ctx_init_modulus:\n"
                    "ix = %wd\n"
                    "p = %{fmpz}\n"
                    "d = %{slong}\n"
                    "Expected type: %s\n"
                    "Got type:      %s\n",
                    ix, p, d, get_str(type), get_str(fq_default_ctx_type(ctx)));

        fq_default_ctx_clear(ctx);
        fmpz_mod_poly_clear(modulus, mod_ctx);
        fmpz_mod_ctx_clear(mod_ctx);
        fmpz_clear(p);
    }

    /* fq_default_ctx_init_modulus_nmod */
    for (ix = 0; ix < 100 * flint_test_multiplier(); ix++)
    {
        int type;
        nmod_poly_t modulus;
        nmod_t mod_ctx;
        fmpz_t p;
        ulong px;
        slong d;
        char var[] = "x";

        type = TYPE_MIN + n_randint(state, TYPE_MAX + TYPE_MIN - 1);
        fmpz_init(p);

        /* Get degree and prime */
        switch (type)
        {
            case FQ_DEFAULT_FQ_ZECH: d = rand_deg_p_fq_zech(state, p); break;
            case FQ_DEFAULT_FQ_NMOD: d = rand_deg_p_fq_nmod(state, p); break;
            case FQ_DEFAULT_FQ: fmpz_clear(p); continue;
            case FQ_DEFAULT_NMOD: d = rand_deg_p_nmod(state, p); break;
            case FQ_DEFAULT_FMPZ_MOD: fmpz_clear(p); continue;
            default: FLINT_UNREACHABLE;
        }

        px = fmpz_get_ui(p);
        fmpz_clear(p);

        nmod_init(&mod_ctx, px);
        nmod_poly_init_mod(modulus, mod_ctx);

        if (type != FQ_DEFAULT_FQ_ZECH)
            nmod_poly_randtest_monic_irreducible(modulus, state, d + 1);
        else
            nmod_poly_randtest_monic_primitive(modulus, state, d + 1);

        fq_default_ctx_init_modulus_nmod(ctx, modulus, var);

        result = (fq_default_ctx_type(ctx) == type);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "fq_default_ctx_init_modulus_nmod:\n"
                    "ix = %wd\n"
                    "p = %{fmpz}\n"
                    "d = %{slong}\n"
                    "Expected type: %s\n"
                    "Got type:      %s\n",
                    ix, p, d, get_str(type), get_str(fq_default_ctx_type(ctx)));

        fq_default_ctx_clear(ctx);
        nmod_poly_clear(modulus);
    }

    TEST_FUNCTION_END(state);
}

#undef TYPE_MIN
#undef TYPE_MAX
