/*
    Copyright (C) 2026 Ricardo Buring

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "perm.h"
#include "fmpz.h"

TEST_FUNCTION_START(perm_next, state)
{
    /* check |S_n| = n! and |A_n| = n! / 2 */

    slong n;
    slong lex_count, lex_count_even;
    slong heap_count, heap_count_even, *perm, *stack;
    fmpz_t fact, half_fact;
    fmpz_init(fact);
    fmpz_init(half_fact);

    for (n = 2; n < 11; n++)
    {
        fmpz_fac_ui(fact, (ulong) n);
        fmpz_divexact_ui(half_fact, fact, (ulong) 2);

        lex_count = lex_count_even = 0;
        perm = _perm_init(n);
        do {
            //flint_printf("%{slong*}\n", perm, n);
            ++lex_count;
            if (!_perm_parity(perm, n)) ++lex_count_even;
        } while (_perm_next_lex(perm, n));

        heap_count = heap_count_even = 0;
        _perm_one(perm, n);
        stack = (slong*) flint_calloc(n, sizeof(slong));
        do {
            //flint_printf("%{slong*}\n", perm, n);
            int odd = _perm_parity(perm, n);
            if (odd != (heap_count % 2)) {
                TEST_FUNCTION_FAIL("unexpected sign in Heap's algorithm\n"
                        "parity = %d\n"
                        "index = %wd\n",
                        odd, heap_count);
            }
            ++heap_count;
            if (!odd) ++heap_count_even;
        } while (_perm_next_heap(perm, n, stack));

        _perm_clear(perm);
        flint_free(stack);

        if (!fmpz_equal_si(fact, lex_count) ||
            !fmpz_equal_si(half_fact, lex_count_even) ||
            !fmpz_equal_si(fact, heap_count) ||
            !fmpz_equal_si(half_fact, heap_count_even))
        {
            TEST_FUNCTION_FAIL(
                    "n = %wd\n"
                    "n! = %{fmpz}\n"
                    "n!/2 = %{fmpz}\n",
                    "lex_count = %wd\n"
                    "lex_count_even = %wd\n"
                    "heap_count = %wd\n"
                    "heap_count_even = %wd\n",
                    n, fact, half_fact,
                    lex_count, lex_count_even,
                    heap_count, heap_count_even);
        }
    }

    fmpz_clear(fact);
    fmpz_clear(half_fact);

    TEST_FUNCTION_END(state);
}
