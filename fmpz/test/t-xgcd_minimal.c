/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, j, result;
    fmpz_t maxval;
    fmpz_t d, x, y, a, b;
    fmpz_t tx, ty;

    FLINT_TEST_INIT(state);
    fmpz_set_str(maxval, "999999999999999999999999999999999999999999", 10);

    flint_printf("xgcd_minimal....");
    fflush(stdout);


    fmpz_init(d);
    fmpz_init(x);
    fmpz_init(y);
    fmpz_init(a);
    fmpz_init(b);

    fmpz_xgcd_minimal(d, x, y, a, b); /* check small (0, 0)
                                       * gives (0, 0, 0) */
    result = (fmpz_is_zero(d)
           && fmpz_is_zero(x)
           && fmpz_is_zero(y)
           && fmpz_is_zero(a)
           && fmpz_is_zero(b));

    if (!result)
    {
        flint_printf("FAIL:\n\n");
        flint_printf("d = "), fmpz_print(d), flint_printf("\n");
        flint_printf("x = "), fmpz_print(x), flint_printf("\n");
        flint_printf("y = "), fmpz_print(y), flint_printf("\n");
        flint_printf("a = "), fmpz_print(a), flint_printf("\n");
        flint_printf("b = "), fmpz_print(b), flint_printf("\n");
        abort();
    }


    _fmpz_promote_val(a);
    _fmpz_promote_val(b);
    fmpz_zero(a);
    fmpz_zero(b);

    fmpz_xgcd_minimal(d, x, y, a, b); /* check promoted (0, 0)
                                       * gives (0, 0, 0) */
    result = (fmpz_is_zero(d)
           && fmpz_is_zero(x)
           && fmpz_is_zero(y));

    if (!result)
    {
        flint_printf("FAIL:\n\n");
        flint_printf("d %d\n",fmpz_is_zero(d));
        flint_printf("x %d\n",fmpz_is_zero(x));
        flint_printf("y %d\n",fmpz_is_zero(y));
        flint_printf("a %d\n",fmpz_is_zero(a));
        flint_printf("b %d\n",fmpz_is_zero(b));
        flint_printf("d = "), fmpz_print(d), flint_printf("\n");
        flint_printf("x = "), fmpz_print(x), flint_printf("\n");
        flint_printf("y = "), fmpz_print(y), flint_printf("\n");
        flint_printf("a = "), fmpz_print(a), flint_printf("\n");
        flint_printf("b = "), fmpz_print(b), flint_printf("\n");
        abort();
    }


    fmpz_one(a);
    fmpz_set_str(b, "999", 10);
    fmpz_xgcd_minimal(d, x, y, a, b); /* check (1, b) for small b gives
                                       * (1, **long expressions**) */

    /* i = |sgn((b - 1) * (b + 1))| */
    fmpz_init_set(tx, b);
    fmpz_init_set(ty, b);
    fmpz_sub_si(tx, tx, 1);
    fmpz_add_ui(ty, ty, 1);
    fmpz_mul(tx, tx, ty);
    i = abs(fmpz_sgn(tx));

    /* j = sgn(b) * (sgn(b + 1) - sgn(b - 1)) */
    fmpz_set(tx, b);
    fmpz_set(ty, b);
    fmpz_sub_si(tx, tx, 1);
    fmpz_add_ui(ty, ty, 1);
    j = fmpz_sgn(b) * (fmpz_sgn(ty) - fmpz_sgn(tx));

    result = (fmpz_is_one(d)
           && fmpz_equal_si(x, i)
           && fmpz_equal_si(y, j));

    if (!result)
    {
        flint_printf("FAIL:\n\n");
        flint_printf("d = "), fmpz_print(d), flint_printf("\n");
        flint_printf("x = "), fmpz_print(x), flint_printf("\n");
        flint_printf("y = "), fmpz_print(y), flint_printf("\n");
        flint_printf("a = "), fmpz_print(a), flint_printf("\n");
        flint_printf("b = "), fmpz_print(b), flint_printf("\n");
        abort();
    }


    fmpz_one(a);
    fmpz_set_str(b, "99999999999999999999999999999999999999999999", 10);
    fmpz_xgcd_minimal(d, x, y, a, b); /* check (1, b) for big b gives
                                       * (1, **long expressions**) */

    /* i = |sgn((b - 1) * (b + 1))| */
    fmpz_init_set(tx, b);
    fmpz_init_set(ty, b);
    fmpz_sub_si(tx, tx, 1);
    fmpz_add_ui(ty, ty, 1);
    fmpz_mul(tx, tx, ty);
    i = abs(fmpz_sgn(tx));

    /* j = sgn(b) * (sgn(b + 1) - sgn(b - 1)) */
    fmpz_set(tx, b);
    fmpz_set(ty, b);
    fmpz_sub_si(tx, tx, 1);
    fmpz_add_ui(ty, ty, 1);
    j = fmpz_sgn(b) * (fmpz_sgn(ty) - fmpz_sgn(tx));

    result = (fmpz_is_one(d)
           && fmpz_equal_si(x, i)
           && fmpz_equal_si(y, j));

    if (!result)
    {
        flint_printf("FAIL:\n\n");
        flint_printf("d = "), fmpz_print(d), flint_printf("\n");
        flint_printf("x = "), fmpz_print(x), flint_printf("\n");
        flint_printf("y = "), fmpz_print(y), flint_printf("\n");
        flint_printf("a = "), fmpz_print(a), flint_printf("\n");
        flint_printf("b = "), fmpz_print(b), flint_printf("\n");
        abort();
    }


    fmpz_set_str(a, "999", 10);
    fmpz_one(b);
    fmpz_xgcd_minimal(d, x, y, a, b); /* check (a, 1) for small a gives
                                       * (1, 0, 1) */

    result = (fmpz_is_one(d)
           && fmpz_is_zero(x)
           && fmpz_is_one(y));

    if (!result)
    {
        flint_printf("FAIL:\n\n");
        flint_printf("d = "), fmpz_print(d), flint_printf("\n");
        flint_printf("x = "), fmpz_print(x), flint_printf("\n");
        flint_printf("y = "), fmpz_print(y), flint_printf("\n");
        flint_printf("a = "), fmpz_print(a), flint_printf("\n");
        flint_printf("b = "), fmpz_print(b), flint_printf("\n");
        abort();
    }


    fmpz_set_str(a, "9771293786198273698712639876198276398716298736999", 10);
    fmpz_one(b);
    fmpz_xgcd_minimal(d, x, y, a, b); /* check (a, 1) for big a gives
                                       * (1, 0, 1) */

    result = (fmpz_is_one(d)
           && fmpz_is_zero(x)
           && fmpz_is_one(y));

    if (!result)
    {
        flint_printf("FAIL:\n\n");
        flint_printf("d = "), fmpz_print(d), flint_printf("\n");
        flint_printf("x = "), fmpz_print(x), flint_printf("\n");
        flint_printf("y = "), fmpz_print(y), flint_printf("\n");
        flint_printf("a = "), fmpz_print(a), flint_printf("\n");
        flint_printf("b = "), fmpz_print(b), flint_printf("\n");
        abort();
    }


    fmpz_set_str(a, "977", 10);
    fmpz_set_str(b, "977", 10);
    fmpz_xgcd_minimal(d, x, y, a, b); /* check (a, a) for small a gives
                                       * (|a|, 0, sgn(a)) */
    
    fmpz_abs(tx, a);
    result = (fmpz_equal(d, tx)
           && fmpz_is_zero(x)
           && fmpz_equal_si(y, fmpz_sgn(a)));

    if (!result)
    {
        flint_printf("FAIL:\n\n");
        flint_printf("d = "), fmpz_print(d), flint_printf("\n");
        flint_printf("x = "), fmpz_print(x), flint_printf("\n");
        flint_printf("y = "), fmpz_print(y), flint_printf("\n");
        flint_printf("a = "), fmpz_print(a), flint_printf("\n");
        flint_printf("b = "), fmpz_print(b), flint_printf("\n");
        abort();
    }


    fmpz_set_str(a, "977817236128736817263876121723681726387618273681", 10);
    fmpz_set_str(b, "977817236128736817263876121723681726387618273681", 10);
    fmpz_xgcd_minimal(d, x, y, a, b); /* check (a, a) for big a gives
                                       * (|a|, 0, sgn(a)) */
    
    fmpz_abs(tx, a);
    result = (fmpz_equal(d, tx)
           && fmpz_is_zero(x)
           && fmpz_equal_si(y, fmpz_sgn(a)));

    if (!result)
    {
        flint_printf("FAIL:\n\n");
        flint_printf("d = "), fmpz_print(d), flint_printf("\n");
        flint_printf("x = "), fmpz_print(x), flint_printf("\n");
        flint_printf("y = "), fmpz_print(y), flint_printf("\n");
        flint_printf("a = "), fmpz_print(a), flint_printf("\n");
        flint_printf("b = "), fmpz_print(b), flint_printf("\n");
        abort();
    }


    fmpz_set_str(a, "97177", 10);
    fmpz_set_str(b, "-97177", 10);
    fmpz_xgcd_minimal(d, x, y, a, b); /* check (a, -a) for small a gives
                                       * (|a|, 0, -sgn(a)) */
    
    fmpz_abs(tx, a);
    result = (fmpz_equal(d, tx)
           && fmpz_is_zero(x)
           && fmpz_equal_si(y, -fmpz_sgn(a)));

    if (!result)
    {
        flint_printf("FAIL:\n\n");
        flint_printf("d = "), fmpz_print(d), flint_printf("\n");
        flint_printf("x = "), fmpz_print(x), flint_printf("\n");
        flint_printf("y = "), fmpz_print(y), flint_printf("\n");
        flint_printf("a = "), fmpz_print(a), flint_printf("\n");
        flint_printf("b = "), fmpz_print(b), flint_printf("\n");
        abort();
    }


    fmpz_set_str(a, "9778172361736817263876121723681726387618273681", 10);
    fmpz_set_str(b, "-9778172361736817263876121723681726387618273681", 10);
    fmpz_xgcd_minimal(d, x, y, a, b); /* check (a, -a) for big a gives
                                       * (|a|, 0, -sgn(a)) */
    
    fmpz_abs(tx, a);
    result = (fmpz_equal(d, tx)
           && fmpz_is_zero(x)
           && fmpz_equal_si(y, -fmpz_sgn(a)));

    if (!result)
    {
        flint_printf("FAIL:\n\n");
        flint_printf("d =akjsdn "), fmpz_print(d), flint_printf("\n");
        flint_printf("x = "), fmpz_print(x), flint_printf("\n");
        flint_printf("y = "), fmpz_print(y), flint_printf("\n");
        flint_printf("a = "), fmpz_print(a), flint_printf("\n");
        flint_printf("b = "), fmpz_print(b), flint_printf("\n");
        abort();
    }



    /* Test the switch (f, g) -> (g, f) is ok */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t d1, d2, a1, a2, b1, b2, f, g, tmp1, tmp2, tmp3, tmp4;

        fmpz_init(d1);
        fmpz_init(a1);
        fmpz_init(b1);
        fmpz_init(d2);
        fmpz_init(a2);
        fmpz_init(b2);
        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(tmp1);
        fmpz_init(tmp2);
        fmpz_init(tmp3);
        fmpz_init(tmp4);

        fmpz_randm(f, state, maxval);
        fmpz_randm(g, state, maxval);
        if (fmpz_is_zero(g)) /* ensure not dividing by 0 */
            fmpz_one(g);
        if (n_randint(state, 2)) fmpz_neg(g, g);
        if (n_randint(state, 2)) fmpz_neg(f, f);

        fmpz_xgcd_minimal(d1, a1, b1, f, g);
        fmpz_xgcd_minimal(d2, b2, a2, g, f);

        fmpz_divexact(tmp1, g, d1);
        fmpz_divexact(tmp2, f, d1);
        fmpz_set(tmp3, a1);
        fmpz_mul(tmp3, tmp3, f);
        fmpz_addmul(tmp3, b1, g);
        fmpz_set(tmp4, a2);
        fmpz_mul(tmp4, tmp4, f);
        fmpz_addmul(tmp4, b2, g);
        result = (fmpz_equal(d1, d2)
               && (fmpz_is_zero(g) || fmpz_is_zero(f)
                 || (fmpz_cmpabs(a1, tmp1) <= 0 && fmpz_cmpabs(b1, tmp2) <= 0)
                 || (fmpz_cmpabs(a2, tmp1) <= 0 && fmpz_cmpabs(b2, tmp2) <= 0))
               && fmpz_equal(tmp3, d1)
               && fmpz_equal(tmp4, d2));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("d1 = "), fmpz_print(d1), flint_printf("\n");
            flint_printf("a1 = "), fmpz_print(a1), flint_printf("\n");
            flint_printf("b1 = "), fmpz_print(b1), flint_printf("\n");
            flint_printf("d2 = "), fmpz_print(d2), flint_printf("\n");
            flint_printf("a2 = "), fmpz_print(a2), flint_printf("\n");
            flint_printf("b2 = "), fmpz_print(b2), flint_printf("\n");
            flint_printf("f = "), fmpz_print(f), flint_printf("\n");
            flint_printf("g = "), fmpz_print(g), flint_printf("\n");
            abort();
        }

        fmpz_clear(d1);
        fmpz_clear(a1);
        fmpz_clear(b1);
        fmpz_clear(d2);
        fmpz_clear(a2);
        fmpz_clear(b2);
        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_init(tmp1);
        fmpz_init(tmp2);
        fmpz_init(tmp3);
        fmpz_init(tmp4);
    }


    /* Test aliasing of d and f, a and g */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t d, a, b, c, f, g, F, G;

        fmpz_init(d);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(F);
        fmpz_init(G);

        fmpz_randtest_unsigned(G, state, 200);
        fmpz_add_ui(G, G, 1);
        fmpz_randm(F, state, G);
        if (n_randint(state, 2)) fmpz_neg(G, G);
        if (n_randint(state, 2)) fmpz_neg(F, F);
        fmpz_set(f, F);
        fmpz_set(g, G);

        fmpz_xgcd_minimal(d, a, b, f, g);
        fmpz_xgcd_minimal(f, g, c, f, g);

        result = (fmpz_equal(d, f)
               && fmpz_equal(b, c)
               && fmpz_equal(a, g));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("d = "), fmpz_print(d), flint_printf("\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("c = "), fmpz_print(c), flint_printf("\n");
            flint_printf("f = "), fmpz_print(f), flint_printf("\n");
            flint_printf("g = "), fmpz_print(g), flint_printf("\n");
            flint_printf("F = "), fmpz_print(F), flint_printf("\n");
            flint_printf("G = "), fmpz_print(G), flint_printf("\n");
            abort();
        }

        fmpz_clear(d);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(F);
        fmpz_clear(G);
    }


    /* Test aliasing of a and f, d and g */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t d, a, b, c, f, g, F, G;

        fmpz_init(d);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(F);
        fmpz_init(G);

        fmpz_randtest_unsigned(G, state, 200);
        fmpz_add_ui(G, G, 1);
        fmpz_randm(F, state, G);
        if (n_randint(state, 2)) fmpz_neg(G, G);
        if (n_randint(state, 2)) fmpz_neg(F, F);
        fmpz_set(f, F);
        fmpz_set(g, G);

        fmpz_xgcd_minimal(d, a, b, f, g);
        fmpz_xgcd_minimal(g, f, c, f, g);

        result = (fmpz_equal(d, g)
               && fmpz_equal(b, c)
               && fmpz_equal(a, f));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("d = "), fmpz_print(d), flint_printf("\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("c = "), fmpz_print(c), flint_printf("\n");
            flint_printf("f = "), fmpz_print(f), flint_printf("\n");
            flint_printf("g = "), fmpz_print(g), flint_printf("\n");
            flint_printf("F = "), fmpz_print(F), flint_printf("\n");
            flint_printf("G = "), fmpz_print(G), flint_printf("\n");
            abort();
        }

        fmpz_clear(d);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(F);
        fmpz_clear(G);
    }


    /* Test aliasing of d and f, b and g */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t d, a, b, c, f, g, F, G;

        fmpz_init(d);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(F);
        fmpz_init(G);

        fmpz_randtest_unsigned(G, state, 200);
        fmpz_add_ui(G, G, 1);
        fmpz_randm(F, state, G);
        if (n_randint(state, 2)) fmpz_neg(G, G);
        if (n_randint(state, 2)) fmpz_neg(F, F);
        fmpz_set(f, F);
        fmpz_set(g, G);

        fmpz_xgcd_minimal(d, a, b, f, g);
        fmpz_xgcd_minimal(f, c, g, f, g);

        result = (fmpz_equal(d, f) 
               && fmpz_equal(a, c)
               && fmpz_equal(b, g));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("d = "), fmpz_print(d), flint_printf("\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("c = "), fmpz_print(c), flint_printf("\n");
            flint_printf("f = "), fmpz_print(f), flint_printf("\n");
            flint_printf("g = "), fmpz_print(g), flint_printf("\n");
            flint_printf("F = "), fmpz_print(F), flint_printf("\n");
            flint_printf("G = "), fmpz_print(G), flint_printf("\n");
            abort();
        }

        fmpz_clear(d);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(F);
        fmpz_clear(G);
    }


    /* Test aliasing of b and f, d and g */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t d, a, b, c, f, g, F, G;

        fmpz_init(d);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(F);
        fmpz_init(G);

        fmpz_randtest_unsigned(G, state, 200);
        fmpz_add_ui(G, G, 1);
        fmpz_randm(F, state, G);
        if (n_randint(state, 2)) fmpz_neg(G, G);
        if (n_randint(state, 2)) fmpz_neg(F, F);
        fmpz_set(f, F);
        fmpz_set(g, G);

        fmpz_xgcd_minimal(d, a, b, f, g);
        fmpz_xgcd_minimal(g, c, f, f, g);

        result = (fmpz_equal(d, g)
               && fmpz_equal(a, c)
               && fmpz_equal(b, f));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("d = "), fmpz_print(d), flint_printf("\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("c = "), fmpz_print(c), flint_printf("\n");
            flint_printf("f = "), fmpz_print(f), flint_printf("\n");
            flint_printf("g = "), fmpz_print(g), flint_printf("\n");
            flint_printf("F = "), fmpz_print(F), flint_printf("\n");
            flint_printf("G = "), fmpz_print(G), flint_printf("\n");
            abort();
        }

        fmpz_clear(d);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(F);
        fmpz_clear(G);
    }


    /* Test a f  + b g == d and d >= 0 */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t d, a, b, f, g, t1, t2;

        fmpz_init(d);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(t1);
        fmpz_init(t2);

        fmpz_randtest_unsigned(g, state, 200);
        fmpz_add_ui(g, g, 1);
        fmpz_randm(f, state, g);
        if (n_randint(state, 2)) fmpz_neg(g, g);
        if (n_randint(state, 2)) fmpz_neg(f, f);
        
        fmpz_xgcd_minimal(d, a, b, f, g);

        fmpz_mul(t1, a, f);
        fmpz_mul(t2, b, g);
        fmpz_add(t1, t1, t2);

        result = fmpz_equal(t1, d) && fmpz_sgn(d) >= 0;
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("d = "), fmpz_print(d), flint_printf("\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("f = "), fmpz_print(f), flint_printf("\n");
            flint_printf("g = "), fmpz_print(g), flint_printf("\n");
            flint_printf("t1 = "), fmpz_print(t1), flint_printf("\n");
            flint_printf("t2 = "), fmpz_print(t2), flint_printf("\n");
            abort();
        }

        fmpz_clear(d);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(t1);
        fmpz_clear(t2);
    }


    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
