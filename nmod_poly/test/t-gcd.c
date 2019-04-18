/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);
    

    flint_printf("gcd....");
    fflush(stdout);

    if (FLINT_BITS == 64)
    {
        nmod_poly_t g, a, b;
        mp_limb_t p = 4611686018427388039;
        nmod_poly_init(g, p);
        nmod_poly_init(a, p);
        nmod_poly_init(b, p);
nmod_poly_set_coeff_ui(a, 2100, UWORD(2006661583596935265));
nmod_poly_set_coeff_ui(a, 2078, UWORD(1089917328005271071));
nmod_poly_set_coeff_ui(a, 2003, UWORD(941394622217952937));
nmod_poly_set_coeff_ui(a, 2000, UWORD(2006661583596935265));
nmod_poly_set_coeff_ui(a, 1102, UWORD(1968267277683696712));
nmod_poly_set_coeff_ui(a, 1101, UWORD(3546399977828707305));
nmod_poly_set_coeff_ui(a, 1100, UWORD(541301798270834258));
nmod_poly_set_coeff_ui(a, 1080, UWORD(3098613960755694602));
nmod_poly_set_coeff_ui(a, 1079, UWORD(1502374134347507339));
nmod_poly_set_coeff_ui(a, 1078, UWORD(2667994756610227708));
nmod_poly_set_coeff_ui(a, 1005, UWORD(2063265285853271828));
nmod_poly_set_coeff_ui(a, 1004, UWORD(3984226001478627868));
nmod_poly_set_coeff_ui(a, 1003, UWORD(650685849759317055));
nmod_poly_set_coeff_ui(a, 1002, UWORD(1968267277683696712));
nmod_poly_set_coeff_ui(a, 1001, UWORD(3546399977828707305));
nmod_poly_set_coeff_ui(a, 1000, UWORD(541301798270834258));
nmod_poly_set_coeff_ui(a, 104, UWORD(2267212073552292930));
nmod_poly_set_coeff_ui(a, 103, UWORD(1699524650453150818));
nmod_poly_set_coeff_ui(a, 102, UWORD(419342695940634527));
nmod_poly_set_coeff_ui(a, 101, UWORD(3942558637419700805));
nmod_poly_set_coeff_ui(a, 100, UWORD(1521));
nmod_poly_set_coeff_ui(a, 82, UWORD(1037763813349290299));
nmod_poly_set_coeff_ui(a, 81, UWORD(145226791149400957));
nmod_poly_set_coeff_ui(a, 80, UWORD(1148529091590782859));
nmod_poly_set_coeff_ui(a, 79, UWORD(2981475204432541474));
nmod_poly_set_coeff_ui(a, 78, UWORD(4035866613882571826));
nmod_poly_set_coeff_ui(a, 7, UWORD(3856387922096858087));
nmod_poly_set_coeff_ui(a, 6, UWORD(128272702618979181));
nmod_poly_set_coeff_ui(a, 5, UWORD(3482400894885422989));
nmod_poly_set_coeff_ui(a, 4, UWORD(865794471887245894));
nmod_poly_set_coeff_ui(a, 3, UWORD(3257437664047763809));
nmod_poly_set_coeff_ui(a, 2, UWORD(419342695940634527));
nmod_poly_set_coeff_ui(a, 1, UWORD(3942558637419700805));
nmod_poly_set_coeff_ui(a, 0, UWORD(1521));

nmod_poly_set_coeff_ui(b, 1200, UWORD(1662416798900483325));
nmod_poly_set_coeff_ui(b, 1178, UWORD(2788122389242055195));
nmod_poly_set_coeff_ui(b, 1156, UWORD(3704061168200567502));
nmod_poly_set_coeff_ui(b, 1103, UWORD(1790409644260772760));
nmod_poly_set_coeff_ui(b, 1100, UWORD(3324833597800966650));
nmod_poly_set_coeff_ui(b, 1081, UWORD(2307667962473032650));
nmod_poly_set_coeff_ui(b, 1078, UWORD(2788122389242055195));
nmod_poly_set_coeff_ui(b, 1006, UWORD(2052993745029114342));
nmod_poly_set_coeff_ui(b, 1003, UWORD(1790409644260772760));
nmod_poly_set_coeff_ui(b, 1000, UWORD(1662416798900483325));
nmod_poly_set_coeff_ui(b, 202, UWORD(672878238860598452));
nmod_poly_set_coeff_ui(b, 201, UWORD(819159959961227504));
nmod_poly_set_coeff_ui(b, 200, UWORD(39));
nmod_poly_set_coeff_ui(b, 180, UWORD(4089995938177947931));
nmod_poly_set_coeff_ui(b, 179, UWORD(312944801058649168));
nmod_poly_set_coeff_ui(b, 178, UWORD(2926679785938335091));
nmod_poly_set_coeff_ui(b, 158, UWORD(4401633899136726058));
nmod_poly_set_coeff_ui(b, 157, UWORD(4359923666897560538));
nmod_poly_set_coeff_ui(b, 156, UWORD(2466580787167147313));
nmod_poly_set_coeff_ui(b, 105, UWORD(3032517747292651873));
nmod_poly_set_coeff_ui(b, 104, UWORD(200562934235634078));
nmod_poly_set_coeff_ui(b, 103, UWORD(1853618366758975553));
nmod_poly_set_coeff_ui(b, 102, UWORD(1345756477721196904));
nmod_poly_set_coeff_ui(b, 101, UWORD(1638319919922455008));
nmod_poly_set_coeff_ui(b, 100, UWORD(78));
nmod_poly_set_coeff_ui(b, 83, UWORD(1046315994897039504));
nmod_poly_set_coeff_ui(b, 82, UWORD(950061260731447402));
nmod_poly_set_coeff_ui(b, 81, UWORD(3204665336321833643));
nmod_poly_set_coeff_ui(b, 80, UWORD(4089995938177947931));
nmod_poly_set_coeff_ui(b, 79, UWORD(312944801058649168));
nmod_poly_set_coeff_ui(b, 78, UWORD(2926679785938335091));
nmod_poly_set_coeff_ui(b, 8, UWORD(3913060101590872477));
nmod_poly_set_coeff_ui(b, 7, UWORD(3229367217583684204));
nmod_poly_set_coeff_ui(b, 6, UWORD(164279788043871567));
nmod_poly_set_coeff_ui(b, 5, UWORD(3032517747292651873));
nmod_poly_set_coeff_ui(b, 4, UWORD(200562934235634078));
nmod_poly_set_coeff_ui(b, 3, UWORD(1853618366758975553));
nmod_poly_set_coeff_ui(b, 2, UWORD(672878238860598452));
nmod_poly_set_coeff_ui(b, 1, UWORD(819159959961227504));
nmod_poly_set_coeff_ui(b, 0, UWORD(39));

printf("a: "); nmod_poly_print_pretty(a, "v"); printf("\n");
printf("b: "); nmod_poly_print_pretty(b, "v"); printf("\n");
printf("calling gcd\n");
        nmod_poly_gcd(g, a, b);
printf("gcd returned\n");
printf("g: "); nmod_poly_print_pretty(g, "v"); printf("\n");

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
    }

    /* 
       Find coprime polys, multiply by another poly 
       and check the GCD is that poly 
    */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, g;

        mp_limb_t n;
        do n = n_randtest_not_zero(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_init(g, n);
        
        do {
            nmod_poly_randtest(a, state, n_randint(state, 1000));
            nmod_poly_randtest(b, state, n_randint(state, 1000));
            nmod_poly_gcd_euclidean(g, a, b);
        } while (g->length != 1);

        do {
            nmod_poly_randtest(c, state, n_randint(state, 1000));
        } while (c->length < 2);
        nmod_poly_make_monic(c, c);
        
        nmod_poly_mul(a, a, c);
        nmod_poly_mul(b, b, c);

        nmod_poly_gcd_euclidean(g, a, b);

        result = (nmod_poly_equal(g, c));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(c), flint_printf("\n\n");
            nmod_poly_print(g), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            abort();
        }
        
        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(g);
    }

    /* Check aliasing of a and g */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, g;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(g, n);
        nmod_poly_randtest(a, state, n_randint(state, 1000));
        nmod_poly_randtest(b, state, n_randint(state, 1000));
        
        nmod_poly_gcd_euclidean(g, a, b);
        nmod_poly_gcd_euclidean(a, a, b);

        result = (nmod_poly_equal(a, g));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(g), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
    }

    /* Check aliasing of b and g */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, g;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(g, n);
        nmod_poly_randtest(a, state, n_randint(state, 1000));
        nmod_poly_randtest(b, state, n_randint(state, 1000));
       
        nmod_poly_gcd_euclidean(g, a, b);
        nmod_poly_gcd_euclidean(b, a, b);

        result = (nmod_poly_equal(b, g));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(g), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
