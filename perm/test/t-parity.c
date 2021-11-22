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
    

    flint_printf("parity....");
    fflush(stdout);

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
            flint_printf("FAIL:\n");
            flint_printf("a: "); _perm_print(a, n); flint_printf("\n\n");
            flint_printf("b: "); _perm_print(b, n); flint_printf("\n\n");
            flint_printf("c: "); _perm_print(c, n); flint_printf("\n\n");
            flint_printf("ap = %d\n", ap);
            flint_printf("bp = %d\n", bp);
            flint_printf("cp = %d\n", cp);
            flint_printf("ap2 = %d\n", ap2);
            flint_printf("bp2 = %d\n", bp2);
            flint_printf("cp2 = %d\n", cp2);
            fflush(stdout);
            flint_abort();
        }

        _perm_clear(a);
        _perm_clear(b);
        _perm_clear(c);
    }
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
