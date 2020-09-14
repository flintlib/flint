/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"
#include "NTL-interface.h"

NTL_CLIENT

int test_ZZ_to_fmpz()
{
    int i, result;
    mp_bitcnt_t bits, randbits;
    fmpz_t int1, int2;
   
    ZZ z;
    FLINT_TEST_INIT(state);
    
    flint_printf("ZZ_to_fmpz....");
    fflush(stdout);

    

    /* Check conversion */
    for (i = 0; i < 10000; i++)
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
           return 0;
        }

        fmpz_clear(int1);
        fmpz_clear(int2);
    }    

    FLINT_TEST_CLEANUP(state);

    return 1;
}

int test_ZZX_to_fmpz_poly()
{
    fmpz_poly_t f_poly1, f_poly2;
    slong length;
    mp_bitcnt_t bits;
    int i, result;
   
    FLINT_TEST_INIT(state);
    
    flint_printf("ZZX_to_fmpz_poly....");
    fflush(stdout);

    

    /* Check aliasing of a and c */
    for (i = 0; i < 10000; i++)
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
           return 0;
        }
          
        fmpz_poly_clear(f_poly1);
        fmpz_poly_clear(f_poly2);
    }
      
    FLINT_TEST_CLEANUP(state);

    return 1;
}

int test_ZZ_pX_to_fmpz_mod_poly()
{
    fmpz_t p;
    fmpz_mod_poly_t f_poly1, f_poly2;
    fmpz_mod_ctx_t ctx;
    slong length;
    int i, result;

    FLINT_TEST_INIT(state);
    
    flint_printf("ZZ_pX_to_fmpz_mod_poly....");
    fflush(stdout);

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Check aliasing of a and c */
    for (i = 0; i < 10000; i++)
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
           return 0;
        }
        
        fmpz_clear(p);
        fmpz_mod_poly_clear(f_poly1, ctx);
        fmpz_mod_poly_clear(f_poly2, ctx);
    }

    fmpz_mod_ctx_clear(ctx);
      
    FLINT_TEST_CLEANUP(state);

    return 1;
}

int test_zz_pX_to_fmpz_mod_poly()
{
    fmpz_mod_poly_t f_poly1, f_poly2;
    fmpz_mod_ctx_t ctx;
    slong length;
    int i, result;
   
    FLINT_TEST_INIT(state);
    
    flint_printf("zz_pX_to_fmpz_mod_poly....");
    fflush(stdout);

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Check aliasing of a and c */
    for (i = 0; i < 10000; i++)
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
           return 0;
        }

        fmpz_mod_poly_clear(f_poly1, ctx);
        fmpz_mod_poly_clear(f_poly2, ctx);
    }
      
    fmpz_mod_ctx_clear(ctx);

    FLINT_TEST_CLEANUP(state);

    return 1;
}

int test_ZZ_pE_to_fq()
{
    fq_t f1, f2;
    int i, result;
    fmpz_mod_ctx_t ctxp;
    fq_ctx_t ctx;
   
    FLINT_TEST_INIT(state);
    
    flint_printf("ZZ_pE_to_fq....");
    fflush(stdout);

    fmpz_mod_ctx_init_ui(ctxp, 2);

    /* Check aliasing of a and c */
    for (i = 0; i < 10000; i++)
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
           return 0;
        }

        fq_clear(f1, ctx);
        fq_clear(f2, ctx);
        fq_ctx_clear(ctx);

        fmpz_mod_poly_clear(fmod, ctxp);
    }

    fmpz_mod_ctx_clear(ctxp);

    FLINT_TEST_CLEANUP(state);

    return 1;
}

int test_ZZ_pEX_to_fq_poly()
{
    fq_poly_t f1, f2;
    slong length;
    int i, result;
    fmpz_mod_ctx_t ctxp;
    fq_ctx_t ctx;
   
    FLINT_TEST_INIT(state);
    
    flint_printf("ZZ_pEX_to_fq_poly....");
    fflush(stdout);

    fmpz_mod_ctx_init_ui(ctxp, 2);

    for (i = 0; i < 10000; i++)
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
           return 0;
        }

        fq_poly_clear(f1, ctx);
        fq_poly_clear(f2, ctx);
        fq_ctx_clear(ctx);

        fmpz_mod_poly_clear(fmod, ctxp);
    }

    fmpz_mod_ctx_clear(ctxp);

    FLINT_TEST_CLEANUP(state);

    return 1;
}

int test_zz_pE_to_fq()
{
    fq_t f1, f2;
    int i, result;
    fmpz_mod_ctx_t ctxp;
    fq_ctx_t ctx;
   
    FLINT_TEST_INIT(state);
    
    flint_printf("zz_pE_to_fq....");
    fflush(stdout);

    fmpz_mod_ctx_init_ui(ctxp, 2);

    /* Check aliasing of a and c */
    for (i = 0; i < 10000; i++)
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
           return 0;
        }
        
        fq_clear(f1, ctx);
        fq_clear(f2, ctx);
        fq_ctx_clear(ctx);

        fmpz_mod_poly_clear(fmod, ctxp);
    }

    fmpz_mod_ctx_clear(ctxp);

    FLINT_TEST_CLEANUP(state);

    return 1;
}

int test_zz_pEX_to_fq_poly()
{
    fq_poly_t f1, f2;
    slong length;
    int i, result;
    fmpz_mod_ctx_t ctxp;
    fq_ctx_t ctx;
   
    FLINT_TEST_INIT(state);
    
    flint_printf("zz_pEX_to_fq_poly....");
    fflush(stdout);

    fmpz_mod_ctx_init_ui(ctxp, 2);

    for (i = 0; i < 10000; i++)
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
           return 0;
        }
        
        fq_poly_clear(f1, ctx);
        fq_poly_clear(f2, ctx);
        fq_ctx_clear(ctx);

        fmpz_mod_poly_clear(fmod, ctxp);
    }
      
    fmpz_mod_ctx_clear(ctxp);

    FLINT_TEST_CLEANUP(state);

    return 1;
}


int
main(void)
{
    int r = 1;
    
    if ((r &= test_ZZ_to_fmpz())) flint_printf("PASS\n");
    if ((r &= test_ZZX_to_fmpz_poly())) flint_printf("PASS\n");
    if ((r &= test_ZZ_pX_to_fmpz_mod_poly())) flint_printf("PASS\n");
    if ((r &= test_ZZ_pE_to_fq())) flint_printf("PASS\n");
    if ((r &= test_ZZ_pEX_to_fq_poly())) flint_printf("PASS\n");
    if ((r &= test_zz_pX_to_fmpz_mod_poly())) flint_printf("PASS\n");
    if ((r &= test_zz_pE_to_fq())) flint_printf("PASS\n");
    if ((r &= test_zz_pEX_to_fq_poly())) flint_printf("PASS\n");

    if (!r) abort();

    return 0;
}
