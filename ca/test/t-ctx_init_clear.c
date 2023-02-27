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
    ca_ctx_t ctx;

    flint_printf("ctx_init_clear....");
    fflush(stdout);

    flint_randinit(state);

    ca_ctx_init(ctx);
    ca_ctx_clear(ctx);

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

