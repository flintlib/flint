/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz_mod.h"
#include "NTL-interface.h"

NTL_CLIENT

TEST_FUNCTION_START(ZZ_to_fmpz, state)
{
    int i, result;
    mp_bitcnt_t bits, randbits;
    fmpz_t int1, int2;

    ZZ z;

    /* Check conversion */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        bits = n_randint(state, 1000) + 1;
        randbits = n_randint(state, bits);

        fmpz_init(int1);
        fmpz_init(int2);

        fmpz_randbits(int1, state, randbits);

        fmpz_get_ZZ(z, int1);
        fmpz_set_ZZ(int2, z);

        result = fmpz_equal(int1, int2);
        if (!result)
        {
           flint_printf("FAIL:\n");
           flint_printf("int1 = %wd  ", *int1); fmpz_print(int1); flint_printf("\n");
           flint_printf("int2 = %wd  ", *int2); fmpz_print(int2); flint_printf("\n");
           flint_abort();
        }

        fmpz_clear(int1);
        fmpz_clear(int2);
    }

    TEST_FUNCTION_END(state);
}

TEST_FUNCTION_START(ZZX_to_fmpz_poly, state)
{
    fmpz_poly_t f_poly1, f_poly2;
    slong length;
    mp_bitcnt_t bits;
    int i, result;

    /* Check aliasing of a and c */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        bits = n_randint(state, 1000) + 1;
        length = n_randint(state, 1000);

        fmpz_poly_init(f_poly1);
        fmpz_poly_init(f_poly2);

        ZZX ZZX_poly;

        fmpz_poly_randtest(f_poly1, state, length, bits);

        fmpz_poly_get_ZZX(ZZX_poly, f_poly1);
        fmpz_poly_set_ZZX(f_poly2, ZZX_poly);

        result = fmpz_poly_equal(f_poly1, f_poly2);
        if (!result)
        {
           flint_printf("FAIL:\n");
           flint_printf("f_poly1 = "); fmpz_poly_print(f_poly1); flint_printf("\n");
           flint_printf("f_poly2 = "); fmpz_poly_print(f_poly2); flint_printf("\n");
           flint_abort();
        }

        fmpz_poly_clear(f_poly1);
        fmpz_poly_clear(f_poly2);
    }

    TEST_FUNCTION_END(state);
}

TEST_FUNCTION_START(ZZ_pX_to_fmpz_mod_poly, state)
{
    fmpz_t p;
    fmpz_mod_poly_t f_poly1, f_poly2;
    fmpz_mod_ctx_t ctx;
    slong length;
    int i, result;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Check aliasing of a and c */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        ZZ_pX ZZ_pX_poly;
        ZZ mod;

        length = n_randint(state, 1000);

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_get_ZZ(mod, p);
        ZZ_p::init(mod);

        fmpz_mod_poly_init(f_poly1, ctx);
        fmpz_mod_poly_init(f_poly2, ctx);

        fmpz_mod_poly_randtest(f_poly1, state, length, ctx);

        fmpz_mod_poly_get_ZZ_pX(ZZ_pX_poly, f_poly1, ctx);
        fmpz_mod_poly_set_ZZ_pX(f_poly2, ZZ_pX_poly, ctx);

        result = fmpz_mod_poly_equal(f_poly1, f_poly2, ctx);
        if (!result)
        {
           flint_printf("FAIL:\n");
           flint_printf("f_poly1 = "); fmpz_mod_poly_print(f_poly1, ctx); flint_printf("\n");
           flint_printf("f_poly2 = "); fmpz_mod_poly_print(f_poly2, ctx); flint_printf("\n");
           flint_abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(f_poly1, ctx);
        fmpz_mod_poly_clear(f_poly2, ctx);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}

TEST_FUNCTION_START(zz_pX_to_fmpz_mod_poly, state)
{
    fmpz_mod_poly_t f_poly1, f_poly2;
    fmpz_mod_ctx_t ctx;
    slong length;
    int i, result;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Check aliasing of a and c */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        zz_pX zz_pX_poly;

        length = n_randint(state, 1000);

        fmpz_mod_ctx_set_modulus_ui(ctx, n_randprime(state, 16, 1));

        zz_p::init(fmpz_get_si(fmpz_mod_ctx_modulus(ctx)));

        fmpz_mod_poly_init(f_poly1, ctx);
        fmpz_mod_poly_init(f_poly2, ctx);

        fmpz_mod_poly_randtest(f_poly1, state, length, ctx);

        fmpz_mod_poly_get_zz_pX(zz_pX_poly, f_poly1, ctx);
        fmpz_mod_poly_set_zz_pX(f_poly2, zz_pX_poly, ctx);

        result = fmpz_mod_poly_equal(f_poly1, f_poly2, ctx);
        if (!result)
        {
           flint_printf("FAIL:\n");
           flint_printf("f_poly1 = "); fmpz_mod_poly_print(f_poly1, ctx); flint_printf("\n");
           flint_printf("f_poly2 = "); fmpz_mod_poly_print(f_poly2, ctx); flint_printf("\n");
           flint_abort();
        }

        fmpz_mod_poly_clear(f_poly1, ctx);
        fmpz_mod_poly_clear(f_poly2, ctx);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}

TEST_FUNCTION_START(ZZ_pE_to_fq, state)
{
    fq_t f1, f2;
    int i, result;
    fmpz_mod_ctx_t ctxp;
    fq_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctxp, 2);

    /* Check aliasing of a and c */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t fmod;
        slong d;
        ZZ prime;
        ZZ_pX mod;

        fmpz_mod_ctx_set_modulus_ui(ctxp, n_randprime(state, 2 + n_randint(state, 6), 1));

        d = n_randint(state, 10) + 1;

        fmpz_get_ZZ(prime, fmpz_mod_ctx_modulus(ctxp));
        ZZ_p::init(prime);

        BuildIrred(mod, d);
        ZZ_pE::init(mod);

        fmpz_mod_poly_init(fmod, ctxp);
        fmpz_mod_poly_set_ZZ_pX(fmod, mod, ctxp);

        fq_ctx_init_modulus(ctx, fmod, ctxp, "a");

        ZZ_pE zzpe;

        fq_init(f1, ctx);
        fq_init(f2, ctx);

        fq_randtest(f1, state, ctx);

        fq_get_ZZ_pE(zzpe, f1, ctx);
        fq_set_ZZ_pE(f2, zzpe, ctx);

        result = fq_equal(f1, f2, ctx);
        if (!result)
        {
           flint_printf("FAIL:\n");
           flint_printf("p = "); fmpz_print(fmpz_mod_ctx_modulus(ctxp)); flint_printf("\n");
           flint_printf("mod = "); fmpz_mod_poly_print_pretty(fmod, "x", ctxp); flint_printf("\n");
           flint_printf("f1 = "); fq_print_pretty(f1, ctx); flint_printf(" - %wd", f1->length); flint_printf("\n");
           flint_printf("zzpe:"); cout << zzpe; flint_printf("\n");
           flint_printf("f2 = "); fq_print_pretty(f2, ctx); flint_printf(" - %wd", f2->length); flint_printf("\n");
           flint_abort();
        }

        fq_clear(f1, ctx);
        fq_clear(f2, ctx);
        fq_ctx_clear(ctx);

        fmpz_mod_poly_clear(fmod, ctxp);
    }

    fmpz_mod_ctx_clear(ctxp);

    TEST_FUNCTION_END(state);
}

TEST_FUNCTION_START(ZZ_pEX_to_fq_poly, state)
{
    fq_poly_t f1, f2;
    slong length;
    int i, result;
    fmpz_mod_ctx_t ctxp;
    fq_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctxp, 2);

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t fmod;
        slong d;
        ZZ prime;
        ZZ_pX mod;

        fmpz_mod_ctx_set_modulus_ui(ctxp, n_randprime(state, 2 + n_randint(state, 6), 1));

        d = n_randint(state, 10) + 1;

        fmpz_get_ZZ(prime, fmpz_mod_ctx_modulus(ctxp));
        ZZ_p::init(prime);

        BuildIrred(mod, d);
        ZZ_pE::init(mod);

        fmpz_mod_poly_init(fmod, ctxp);
        fmpz_mod_poly_set_ZZ_pX(fmod, mod, ctxp);

        fq_ctx_init_modulus(ctx, fmod, ctxp, "a");

        ZZ_pEX zzpex;

        fq_poly_init(f1, ctx);
        fq_poly_init(f2, ctx);

        length = n_randint(state, 1000);

        fq_poly_randtest(f1, state, length, ctx);

        fq_poly_get_ZZ_pEX(zzpex, f1, ctx);
        fq_poly_set_ZZ_pEX(f2, zzpex, ctx);

        result = fq_poly_equal(f1, f2, ctx);
        if (!result)
        {
           flint_printf("FAIL:\n");
           flint_printf("p = "); fmpz_print(fmpz_mod_ctx_modulus(ctxp)); flint_printf("\n");
           flint_printf("mod = "); fmpz_mod_poly_print_pretty(fmod, "x", ctxp); flint_printf("\n");
           flint_printf("f1 = "); fq_poly_print_pretty(f1, "x", ctx); flint_printf("\n");
           flint_printf("zzpex:"); cout << zzpex; flint_printf("\n");
           flint_printf("f2 = "); fq_poly_print_pretty(f2, "x", ctx); flint_printf("\n");
           flint_abort();
        }

        fq_poly_clear(f1, ctx);
        fq_poly_clear(f2, ctx);
        fq_ctx_clear(ctx);

        fmpz_mod_poly_clear(fmod, ctxp);
    }

    fmpz_mod_ctx_clear(ctxp);

    TEST_FUNCTION_END(state);
}

TEST_FUNCTION_START(zz_pE_to_fq, state)
{
    fq_t f1, f2;
    int i, result;
    fmpz_mod_ctx_t ctxp;
    fq_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctxp, 2);

    /* Check aliasing of a and c */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t fmod;
        slong d;
        zz_pX mod;

        fmpz_mod_ctx_set_modulus_ui(ctxp, n_randprime(state, 2 + n_randint(state, 6), 1));

        d = n_randint(state, 10) + 1;

        zz_p::init(fmpz_get_si(fmpz_mod_ctx_modulus(ctxp)));

        BuildIrred(mod, d);
        zz_pE::init(mod);

        fmpz_mod_poly_init(fmod, ctxp);
        fmpz_mod_poly_set_zz_pX(fmod, mod, ctxp);

        fq_ctx_init_modulus(ctx, fmod, ctxp, "a");

        zz_pE zzpe;

        fq_init(f1, ctx);
        fq_init(f2, ctx);

        fq_randtest(f1, state, ctx);

        fq_get_zz_pE(zzpe, f1, ctx);
        fq_set_zz_pE(f2, zzpe, ctx);

        result = fq_equal(f1, f2, ctx);
        if (!result)
        {
           flint_printf("FAIL:\n");
           flint_printf("p = "); fmpz_print(fmpz_mod_ctx_modulus(ctxp)); flint_printf("\n");
           flint_printf("mod = "); fmpz_mod_poly_print_pretty(fmod, "x", ctxp); flint_printf("\n");
           flint_printf("f1 = "); fq_print_pretty(f1, ctx); flint_printf(" - %wd", f1->length); flint_printf("\n");
           flint_printf("zzpe:"); cout << zzpe; flint_printf("\n");
           flint_printf("f2 = "); fq_print_pretty(f2, ctx); flint_printf(" - %wd", f2->length); flint_printf("\n");
           flint_abort();
        }

        fq_clear(f1, ctx);
        fq_clear(f2, ctx);
        fq_ctx_clear(ctx);

        fmpz_mod_poly_clear(fmod, ctxp);
    }

    fmpz_mod_ctx_clear(ctxp);

    TEST_FUNCTION_END(state);
}

TEST_FUNCTION_START(zz_pEX_to_fq_poly, state)
{
    fq_poly_t f1, f2;
    slong length;
    int i, result;
    fmpz_mod_ctx_t ctxp;
    fq_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctxp, 2);

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t fmod;
        slong d;
        zz_pX mod;

        fmpz_mod_ctx_set_modulus_ui(ctxp, n_randprime(state, 2 + n_randint(state, 6), 1));

        d = n_randint(state, 10) + 1;

        zz_p::init(fmpz_get_si(fmpz_mod_ctx_modulus(ctxp)));

        BuildIrred(mod, d);
        zz_pE::init(mod);

        fmpz_mod_poly_init(fmod, ctxp);
        fmpz_mod_poly_set_zz_pX(fmod, mod, ctxp);

        fq_ctx_init_modulus(ctx, fmod, ctxp, "a");

        zz_pEX zzpex;

        fq_poly_init(f1, ctx);
        fq_poly_init(f2, ctx);

        length = n_randint(state, 1000);

        fq_poly_randtest(f1, state, length, ctx);

        fq_poly_get_zz_pEX(zzpex, f1, ctx);
        fq_poly_set_zz_pEX(f2, zzpex, ctx);

        result = fq_poly_equal(f1, f2, ctx);
        if (!result)
        {
           flint_printf("FAIL:\n");
           flint_printf("p = "); fmpz_print(fmpz_mod_ctx_modulus(ctxp)); flint_printf("\n");
           flint_printf("mod = "); fmpz_mod_poly_print_pretty(fmod, "x", ctxp); flint_printf("\n");
           flint_printf("f1 = "); fq_poly_print_pretty(f1, "x", ctx); flint_printf("\n");
           flint_printf("zzpex:"); cout << zzpex; flint_printf("\n");
           flint_printf("f2 = "); fq_poly_print_pretty(f2, "x", ctx); flint_printf("\n");
           flint_abort();
        }

        fq_poly_clear(f1, ctx);
        fq_poly_clear(f2, ctx);
        fq_ctx_clear(ctx);

        fmpz_mod_poly_clear(fmod, ctxp);
    }

    fmpz_mod_ctx_clear(ctxp);

    TEST_FUNCTION_END(state);
}

/* ISO C++ forbids converting a string constant to char* */
char str1[] = "ZZ_to_fmpz";
char str2[] = "ZZX_to_fmpz_poly";
char str3[] = "ZZ_pX_to_fmpz_mod_poly";
char str4[] = "zz_pX_to_fmpz_mod_poly";
char str5[] = "ZZ_pE_to_fq";
char str6[] = "ZZ_pEX_to_fq_poly";
char str7[] = "zz_pE_to_fq";
char str8[] = "zz_pEX_to_fq_poly";

test_struct tests[] =
{
    { CAT(test, ZZ_to_fmpz), str1 },
    { CAT(test, ZZX_to_fmpz_poly), str2 },
    { CAT(test, ZZ_pX_to_fmpz_mod_poly), str3 },
    { CAT(test, zz_pX_to_fmpz_mod_poly), str4 },
    { CAT(test, ZZ_pE_to_fq), str5 },
    { CAT(test, ZZ_pEX_to_fq_poly), str6 },
    { CAT(test, zz_pE_to_fq), str7 },
    { CAT(test, zz_pEX_to_fq_poly), str8 }
};

TEST_MAIN(tests)
