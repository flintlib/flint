/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

int main()
{
    flint_rand_t state;

    flint_printf("field_init_clear....");
    fflush(stdout);

    flint_randinit(state);

    {
        ca_field_t K;
        ca_extension_t I, Pi;
        qqbar_t t;

        qqbar_init(t);
        qqbar_i(t);

        ca_extension_init_qqbar(I, t);
        ca_extension_init_const(Pi, CA_Pi);

        ca_field_init_mpoly_q(K, 2);
        ca_field_set_ext(K, 0, I);
        ca_field_set_ext(K, 1, Pi);

        flint_printf("\n");
        ca_field_print(K);
        flint_printf("\n");

        ca_extension_clear(I);
        ca_extension_clear(Pi);
        ca_field_clear(K);
        qqbar_clear(t);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

