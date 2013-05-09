/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "perm.h"
#include "ulong_extras.h"

int main(void)
{
    int i;
    flint_rand_t state;
    flint_randinit(state);

    printf("inv....");
    fflush(stdout);

    /* check inv(inv(a)) == a */
    for (i = 0; i < 10000; i++)
    {
        len_t n, *a, *b, *c;

        n = n_randint(state, 100);

        a = _perm_init(n);
        b = _perm_init(n);
        c = _perm_init(n);

        _perm_randtest(a, n, state);

        _perm_inv(b, a, n);
        _perm_inv(c, b, n);

        if (!_perm_equal(a, c, n))
        {
            printf("FAIL:\n");
            printf("a: "); _perm_print(a, n); printf("\n\n");
            printf("b: "); _perm_print(b, n); printf("\n\n");
            printf("c: "); _perm_print(c, n); printf("\n\n");
            abort();
        }

        _perm_clear(a);
        _perm_clear(b);
        _perm_clear(c);
    }

    /* check aliasing */
    for (i = 0; i < 10000; i++)
    {
        len_t n, *a, *b;

        n = n_randint(state, 100);

        a = _perm_init(n);
        b = _perm_init(n);

        _perm_randtest(a, n, state);

        _perm_inv(b, a, n);
        _perm_inv(a, a, n);

        if (!_perm_equal(a, b, n))
        {
            printf("FAIL:\n");
            printf("a: "); _perm_print(a, n); printf("\n\n");
            printf("b: "); _perm_print(b, n); printf("\n\n");
            abort();
        }

        _perm_clear(a);
        _perm_clear(b);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
