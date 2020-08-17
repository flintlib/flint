/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "perm.h"
#include "ulong_extras.h"

int main(void)
{
    int i;
    FLINT_TEST_INIT(state);
    

    flint_printf("compose....");
    fflush(stdout);

    /* check (b^(-1))(b(a)) = a */
    for (i = 0; i < 10000; i++)
    {
        slong n, *a, *b, *binv, *c;

        n = n_randint(state, 100);

        a = _perm_init(n);
        b = _perm_init(n);
        binv = _perm_init(n);
        c = _perm_init(n);

        _perm_randtest(a, n, state);
        _perm_randtest(b, n, state);
        _perm_inv(binv, b, n);

        _perm_compose(c, b, a, n);
        _perm_compose(c, binv, c, n);

        if (!_perm_equal(a, c, n))
        {
            flint_printf("FAIL:\n");
            flint_printf("a: "); _perm_print(a, n); flint_printf("\n\n");
            flint_printf("b: "); _perm_print(b, n); flint_printf("\n\n");
            flint_printf("binv: "); _perm_print(binv, n); flint_printf("\n\n");
            flint_printf("c: "); _perm_print(c, n); flint_printf("\n\n");
            abort();
        }

        _perm_clear(a);
        _perm_clear(b);
        _perm_clear(binv);
        _perm_clear(c);
    }

    /* check aliasing with first argument */
    for (i = 0; i < 10000; i++)
    {
        slong n, *a, *b, *c;

        n = n_randint(state, 100);

        a = _perm_init(n);
        b = _perm_init(n);
        c = _perm_init(n);

        _perm_randtest(a, n, state);
        _perm_randtest(b, n, state);

        _perm_compose(c, b, a, n);
        _perm_compose(b, b, a, n);

        if (!_perm_equal(b, c, n))
        {
            flint_printf("FAIL:\n");
            flint_printf("a: "); _perm_print(a, n); flint_printf("\n\n");
            flint_printf("b: "); _perm_print(b, n); flint_printf("\n\n");
            flint_printf("c: "); _perm_print(c, n); flint_printf("\n\n");
            abort();
        }

        _perm_clear(a);
        _perm_clear(b);
        _perm_clear(c);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
