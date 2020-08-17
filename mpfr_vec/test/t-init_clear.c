/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "mpfr_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);
    

    flint_printf("init/clear....");
    fflush(stdout);

    for (i = 0; i < 10000; i++)
    {
        __mpfr_struct *a;
        slong j, length = n_randint(state, 100);
        mp_prec_t prec = n_randint(state, 200) + MPFR_PREC_MIN;

        a = _mpfr_vec_init(length, prec);

        for (j = 0; j < length; j++)
            mpfr_set_ui(a + j, 0, GMP_RNDN);

        _mpfr_vec_clear(a, length);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
