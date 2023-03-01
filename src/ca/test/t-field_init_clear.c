/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"
#include "ca_ext.h"
#include "ca_field.h"

int main()
{
    flint_rand_t state;

    flint_printf("field_init_clear....");
    fflush(stdout);

    flint_randinit(state);

    {
        ca_ctx_t ctx;
        ca_field_t K;
        ca_field_t I, Pi;
        qqbar_t t;

        qqbar_init(t);
        qqbar_i(t);

        ca_ctx_init(ctx);
        ca_field_init_nf(I, t, ctx);
        ca_field_init_const(Pi, CA_Pi, ctx);
        ca_field_init_multi(K, 2, ctx);

/*
        ca_field_set_ext(K, 0, I);
        ca_field_set_ext(K, 1, Pi);
        flint_printf("\n");
        ca_field_print(K);
        flint_printf("\n");
*/

        ca_field_clear(I, ctx);
        ca_field_clear(Pi, ctx);
        ca_field_clear(K, ctx);

        ca_ctx_clear(ctx);
        qqbar_clear(t);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

