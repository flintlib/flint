/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 William Hart
                  2020 Julian Rüth

******************************************************************************/

#include <stdio.h>
#include "nf.h"
#include "nf_elem.h"

int
main(void)
{
    int i;
    flint_rand_t state;

    flint_printf("init/clear....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100 * antic_test_multiplier(); i++)
    {
        nf_t nf;
        nf_elem_t a;

        nf_init_randtest(nf, state, 40, 200);

        nf_elem_init(a, nf);

        if (!nf_elem_is_zero(a, nf))
        {
            flint_printf("FAIL\n");
            abort();
        }

        nf_elem_randtest(a, state, 200, nf);
        nf_elem_clear(a, nf);
        
        nf_clear(nf);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
