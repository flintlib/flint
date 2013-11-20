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
    flint_rand_t state;

    printf("zero....");
    fflush(stdout);

    flint_randinit(state);

    /* Check it's zero */
    for (i = 0; i < 100; i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, struct) * a;
        slong len = n_randint(state, 100);

        TEMPLATE(T, ctx_randtest) (ctx, state);

        a = _TEMPLATE(T, vec_init) (len, ctx);
        _TEMPLATE(T, vec_randtest) (a, state, len, ctx);

        _TEMPLATE(T, vec_zero) (a, len, ctx);

        result = (_TEMPLATE(T, vec_is_zero) (a, len, ctx));
        if (!result)
        {
            printf("FAIL:\n");
            _TEMPLATE(T, vec_print) (a, len, ctx), printf("\n\n");
            abort();
        }

        _TEMPLATE(T, vec_clear) (a, len, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return 0;
}


#endif
