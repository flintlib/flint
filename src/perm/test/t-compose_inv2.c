/*
    Copyright (C) 2026 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "perm.h"

TEST_FUNCTION_START(perm_compose_inv2, state)
{
    int i;

    /* check perm_compose_inv2(a, b) = perm_compose(a, inv(b)) without aliasing */
    for (i = 0; i < 10000; i++)
    {
        slong n, *a, *b, *binv, *c, *d;

        n = n_randint(state, 100);

        a = _perm_init(n);
        b = _perm_init(n);
        binv = _perm_init(n);
        c = _perm_init(n);
        d = _perm_init(n);

        _perm_randtest(a, n, state);
        _perm_randtest(b, n, state);
        _perm_inv(binv, b, n);

        _perm_compose_inv2(c, a, b, n);
        _perm_compose(d, a, binv, n);

        if (!_perm_equal(c, d, n))
            TEST_FUNCTION_FAIL("FAIL (1):\n"
                    "n = %wd\n"
                    "a = %{slong*}\n\n"
                    "b = %{slong*}\n\n"
                    "binv = %{slong*}\n\n"
                    "c = %{slong*}\n\n"
                    "d = %{slong*}\n\n",
                    n,
                    a, n,
                    b, n,
                    binv, n,
                    c, n,
                    d, n);

        _perm_clear(a);
        _perm_clear(b);
        _perm_clear(binv);
        _perm_clear(c);
        _perm_clear(d);
    }

    /* check aliasing the first argument */
    for (i = 0; i < 10000; i++)
    {
        slong n, *a, *b, *c;

        n = n_randint(state, 100);

        a = _perm_init(n);
        b = _perm_init(n);
        c = _perm_init(n);

        _perm_randtest(a, n, state);
        _perm_randtest(b, n, state);

        _perm_compose_inv2(c, a, b, n);
        _perm_compose_inv2(a, a, b, n);

        if (!_perm_equal(c, a, n))
            TEST_FUNCTION_FAIL("FAIL (2):\n"
                    "n = %wd\n"
                    "a = %{slong*}\n\n"
                    "b = %{slong*}\n\n"
                    "c = %{slong*}\n\n",
                    n,
                    a, n,
                    b, n,
                    c, n);

        _perm_clear(a);
        _perm_clear(b);
        _perm_clear(c);
    }

    /* check aliasing the second argument */
    for (i = 0; i < 10000; i++)
    {
        slong n, *a, *b, *c;

        n = n_randint(state, 100);

        a = _perm_init(n);
        b = _perm_init(n);
        c = _perm_init(n);

        _perm_randtest(a, n, state);
        _perm_randtest(b, n, state);

        _perm_compose_inv2(c, a, b, n);
        _perm_compose_inv2(b, a, b, n);

        if (!_perm_equal(c, b, n))
            TEST_FUNCTION_FAIL("FAIL (3):\n"
                    "n = %wd\n"
                    "a = %{slong*}\n\n"
                    "b = %{slong*}\n\n"
                    "c = %{slong*}\n\n",
                    n,
                    a, n,
                    b, n,
                    c, n);

        _perm_clear(a);
        _perm_clear(b);
        _perm_clear(c);
    }

    /* check aliasing both arguments */
    for (i = 0; i < 10000; i++)
    {
        slong n, *a, *c;

        n = n_randint(state, 100);

        a = _perm_init(n);
        c = _perm_init(n);

        _perm_randtest(a, n, state);

        _perm_compose_inv2(c, a, a, n);
        _perm_compose_inv2(a, a, a, n);

        if (!_perm_equal(c, a, n))
            TEST_FUNCTION_FAIL("FAIL (4):\n"
                    "n = %wd\n"
                    "a = %{slong*}\n\n"
                    "c = %{slong*}\n\n",
                    n,
                    a, n,
                    c, n);

        _perm_clear(a);
        _perm_clear(c);
    }

    TEST_FUNCTION_END(state);
}
