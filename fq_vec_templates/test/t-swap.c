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

    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

#include <stdio.h>
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    printf("swap....");
    fflush(stdout);

    for (i = 0; i < 100; i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, struct) * a, *b, *c;
        slong len = n_randint(state, 100);

        TEMPLATE(T, ctx_randtest) (ctx, state);

        a = _TEMPLATE(T, vec_init) (len, ctx);
        b = _TEMPLATE(T, vec_init) (len, ctx);
        c = _TEMPLATE(T, vec_init) (len, ctx);

        _TEMPLATE(T, vec_randtest) (a, state, len, ctx);
        _TEMPLATE(T, vec_randtest) (b, state, len, ctx);

        _TEMPLATE(T, vec_set) (c, b, len, ctx);
        _TEMPLATE(T, vec_swap) (a, b, len, ctx);

        result = (_TEMPLATE(T, vec_equal) (a, c, len, ctx));
        if (!result)
        {
            printf("FAIL:\n");
            _TEMPLATE(T, vec_print) (a, len, ctx), printf("\n\n");
            _TEMPLATE(T, vec_print) (b, len, ctx), printf("\n\n");
            _TEMPLATE(T, vec_print) (c, len, ctx), printf("\n\n");
            abort();
        }

        _TEMPLATE(T, vec_clear) (a, len, ctx);
        _TEMPLATE(T, vec_clear) (b, len, ctx);
        _TEMPLATE(T, vec_clear) (c, len, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);
    printf("PASS\n");
    return 0;
}


#endif
