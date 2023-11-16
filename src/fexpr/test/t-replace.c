/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "calcium.h"
#include "fexpr.h"
#include "fexpr_builtin.h"

void
fexpr_randtest_atom1(fexpr_t expr, flint_rand_t state)
{
    if (n_randint(state, 2))
    {
        if (n_randint(state, 10) == 5)
            fexpr_set_si(expr, (slong) n_randtest(state));
        else
            fexpr_set_si(expr, ((slong) n_randint(state, 3)) - 1);
    }
    else
    {
        if (n_randint(state, 2))
            fexpr_set_symbol_builtin(expr, FEXPR_alpha);
        else
            fexpr_set_symbol_builtin(expr, FEXPR_beta);
    }
}

void
fexpr_randtest1(fexpr_t expr, flint_rand_t state, slong max_leaves, slong max_depth)
{
    if (max_leaves <= 1 || max_depth <= 1 || n_randint(state, 10) == 0)
    {
        fexpr_randtest_atom1(expr, state);
    }
    else
    {
        slong i, nargs;
        fexpr_t f;
        fexpr_ptr args;

        fexpr_init(f);
        fexpr_randtest1(f, state, max_leaves / 4, max_depth - 1);
        max_leaves -= fexpr_num_leaves(f);

        nargs = n_randint(state, 1 + FLINT_MIN(5, max_leaves / max_depth));

        args = _fexpr_vec_init(nargs);

        for (i = 0; i < nargs; i++)
        {
            fexpr_randtest1(args + i, state, max_leaves, max_depth - 1);
            max_leaves -= fexpr_num_leaves(args + i);
        }

        fexpr_call_vec(expr, f, args, nargs);

        _fexpr_vec_clear(args, nargs);
        fexpr_clear(f);
    }
}

int
fexpr_replace_vec_naive(fexpr_t res, const fexpr_t expr, const fexpr_vec_t xs, const fexpr_vec_t ys)
{
    slong i, num_rules;

    num_rules = xs->length;

    if (num_rules != ys->length)
    {
        flint_printf("fexpr_replace_vec: vectors don't match\n");
        flint_abort();
    }

    for (i = 0; i < num_rules; i++)
    {
        if (fexpr_equal(expr, xs->entries + i))
        {
            fexpr_set(res, ys->entries + i);
            return 1;
        }
    }

    if (fexpr_is_atom(expr))
    {
        fexpr_set(res, expr);
        return 0;
    }
    else
    {
        fexpr_ptr args;
        fexpr_t func, view;
        slong nargs;
        int changed;

        nargs = fexpr_nargs(expr);
        args = NULL;

        fexpr_init(func);

        fexpr_view_func(view, expr);
        changed = fexpr_replace_vec_naive(func, view, xs, ys);

        if (nargs >= 1)
        {
            args = _fexpr_vec_init(nargs);
            fexpr_view_arg(view, expr, 0);
            for (i = 0; i < nargs; i++)
            {
                changed |= fexpr_replace_vec_naive(args + i, view, xs, ys);
                if (i < nargs - 1)
                    fexpr_view_next(view);
            }

            fexpr_call_vec(res, func, args, nargs);
            _fexpr_vec_clear(args, nargs);
        }
        else
        {
            fexpr_call0(res, func);
        }

        fexpr_clear(func);

        return changed;
    }
}

TEST_FUNCTION_START(fexpr_replace, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fexpr_t expr, res1, res2, res3;
        fexpr_vec_t xs, ys;
        slong i, len;
        int changed1, changed2, changed3;

        fexpr_init(expr);
        fexpr_init(res1);
        fexpr_init(res2);
        fexpr_init(res3);

        fexpr_randtest1(expr, state, 1 + n_randint(state, 100), 2 + n_randint(state, 10));
        fexpr_randtest1(res1, state, 5, 2 + n_randint(state, 5));
        fexpr_randtest1(res2, state, 5, 2 + n_randint(state, 5));

        len = n_randint(state, 4);
        fexpr_vec_init(xs, len);
        fexpr_vec_init(ys, len);
        for (i = 0; i < len; i++)
        {
            fexpr_randtest1(xs->entries + i, state, 1 + n_randint(state, 3), 1 + n_randint(state, 2));
            fexpr_randtest1(ys->entries + i, state, 1 + n_randint(state, 10), 2 + n_randint(state, 5));
        }

        if (xs->length == 1 && n_randint(state, 2))
        {
            changed1 = fexpr_replace(res1, expr, xs->entries, ys->entries);

            /* also test contains() */
            if (n_randint(state, 2))
                changed1 = fexpr_contains(expr, xs->entries);
        }
        else
            changed1 = fexpr_replace_vec(res1, expr, xs, ys);

        fexpr_set(res2, expr);

        if (xs->length == 1 && n_randint(state, 2))
            changed2 = fexpr_replace(res2, res2, xs->entries, ys->entries);
        else
            changed2 = fexpr_replace_vec(res2, res2, xs, ys);

        changed3 = fexpr_replace_vec_naive(res3, expr, xs, ys);

        if ((changed1 != changed3) || !fexpr_equal(res1, res3))
        {
            flint_printf("FAIL\n");
            flint_printf("expr = "); fexpr_print(expr); flint_printf("\n\n");
            flint_printf("xs = "); fexpr_vec_print(xs); flint_printf("\n\n");
            flint_printf("ys = "); fexpr_vec_print(ys); flint_printf("\n\n");
            flint_printf("res1 = "); fexpr_print(res1); flint_printf("\n\n");
            flint_printf("res3 = "); fexpr_print(res3); flint_printf("\n\n");
            flint_printf("changed1 = %d, changed3 = %d\n\n", changed1, changed3);
            flint_abort();
        }

        if ((changed1 != changed2) || !fexpr_equal(res1, res2))
        {
            flint_printf("FAIL (aliasing)\n");
            flint_printf("expr = "); fexpr_print(expr); flint_printf("\n\n");
            flint_printf("xs = "); fexpr_vec_print(xs); flint_printf("\n\n");
            flint_printf("ys = "); fexpr_vec_print(ys); flint_printf("\n\n");
            flint_printf("res1 = "); fexpr_print(res1); flint_printf("\n\n");
            flint_printf("res2 = "); fexpr_print(res2); flint_printf("\n\n");
            flint_printf("changed1 = %d, changed2 = %d\n\n", changed1, changed2);
            flint_abort();
        }

        fexpr_vec_clear(xs);
        fexpr_vec_clear(ys);

        fexpr_clear(expr);
        fexpr_clear(res1);
        fexpr_clear(res2);
        fexpr_clear(res3);
    }

    TEST_FUNCTION_END(state);
}
