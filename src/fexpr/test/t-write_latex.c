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

#ifdef __GNUC__
# define strlen __builtin_strlen
#else
# include <string.h>
#endif

void
fexpr_randtest_atom(fexpr_t expr, flint_rand_t state)
{
    if (n_randint(state, 2))
    {
        fexpr_set_si(expr, (slong) n_randtest(state));
    }
    else if (n_randint(state, 2))
    {
        switch (n_randint(state, 10))
        {
            case 0: fexpr_set_symbol_builtin(expr, FEXPR_For); break;
            case 1: fexpr_set_symbol_builtin(expr, FEXPR_Where); break;
            case 2: fexpr_set_symbol_builtin(expr, FEXPR_Def); break;
            case 3: fexpr_set_symbol_builtin(expr, FEXPR_Neg); break;
            case 4: fexpr_set_symbol_builtin(expr, FEXPR_Pos); break;
            case 5: fexpr_set_symbol_builtin(expr, FEXPR_Add); break;
            case 6: fexpr_set_symbol_builtin(expr, FEXPR_Sub); break;
            case 7: fexpr_set_symbol_builtin(expr, FEXPR_Mul); break;
            case 8: fexpr_set_symbol_builtin(expr, FEXPR_Div); break;
            default: fexpr_set_symbol_builtin(expr, FEXPR_Pow); break;
        }
    }
    else
    {
        fexpr_set_symbol_builtin(expr, n_randint(state, FEXPR_BUILTIN_LENGTH));
    }
}

/* todo: doesn't actually respect max_leaves */
void
fexpr_randtest_gibberish(fexpr_t expr, flint_rand_t state, slong max_leaves, slong max_depth)
{
    if (max_leaves <= 1 || max_depth <= 1 || n_randint(state, 10) == 0)
    {
        fexpr_randtest_atom(expr, state);
    }
    else
    {
        slong i, nargs;
        fexpr_t f;
        fexpr_ptr args;

        fexpr_init(f);
        if (n_randint(state, 10) != 0)
        {
            fexpr_set_symbol_builtin(f, n_randint(state, FEXPR_BUILTIN_LENGTH));
        }
        else
        {
            fexpr_randtest_gibberish(f, state, max_leaves / 4, max_depth - 1);
            max_leaves -= fexpr_num_leaves(f);
        }

        nargs = n_randint(state, 1 + FLINT_MIN(5, max_leaves / max_depth));

        args = _fexpr_vec_init(nargs);

        for (i = 0; i < nargs; i++)
        {
            fexpr_randtest_gibberish(args + i, state, max_leaves, max_depth - 1);
            max_leaves -= fexpr_num_leaves(args + i);
        }

        fexpr_call_vec(expr, f, args, nargs);

        _fexpr_vec_clear(args, nargs);
        fexpr_clear(f);
    }
}

TEST_FUNCTION_START(fexpr_write_latex, state)
{
    slong iter;

    /* Generate gibberish and just check that we get valid strings */
    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fexpr_t expr;
        char * s;
        ulong flags;
        slong len;

        fexpr_init(expr);
        fexpr_randtest_gibberish(expr, state, 2 + n_randint(state, 30), 2 + n_randint(state, 5));

        flags = 0;
        if (n_randint(state, 2))
            flags |= FEXPR_LATEX_SMALL;
        if (n_randint(state, 2))
            flags |= FEXPR_LATEX_LOGIC;

        /* fexpr_print(expr); printf("\n\n"); */

        s = fexpr_get_str_latex(expr, flags);
        len = strlen(s);
        (void)len;
        flint_free(s);

        fexpr_clear(expr);
    }

    TEST_FUNCTION_END(state);
}
