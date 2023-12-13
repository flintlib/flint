/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "perm.h"

TEST_FUNCTION_START(perm_parity, state)
{
    int i;

    /* check inv(inv(a)) == a */
    for (i = 0; i < 10000; i++)
    {
        slong n, *a, *b, *c;
        int ap, bp, cp, ap2, bp2, cp2;

        n = n_randint(state, 100);

        a = _perm_init(n);
        b = _perm_init(n);
        c = _perm_init(n);

        ap = _perm_randtest(a, n, state);
        bp = _perm_randtest(b, n, state);

        _perm_compose(c, a, b, n);
        cp = ap ^ bp;

        ap2 = _perm_parity(a, n);
        bp2 = _perm_parity(b, n);
        cp2 = _perm_parity(c, n);

        if (ap != ap2 || bp != bp2 || cp != cp2)
        {
            flint_throw(FLINT_TEST_FAIL,
                    "n = %wd\n"
                    "a = %{slong*}\n\n"
                    "b = %{slong*}\n\n"
                    "c = %{slong*}\n\n"
                    "ap = %d\n"
                    "bp = %d\n"
                    "cp = %d\n"
                    "ap2 = %d\n"
                    "bp2 = %d\n"
                    "cp2 = %d\n",
                    n,
                    a, n,
                    b, n,
                    c, n,
                    ap,
                    bp,
                    cp,
                    ap2,
                    bp2,
                    cp2);
        }

        _perm_clear(a);
        _perm_clear(b);
        _perm_clear(c);
    }

    TEST_FUNCTION_END(state);
}
