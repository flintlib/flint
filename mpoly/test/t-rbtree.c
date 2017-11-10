/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "mpoly.h"
#include "ulong_extras.h"

int
main(void)
{
    int new;
    slong i;
    void ** data;
    slong * keys;
    mpoly_rbtree_t t;
    mpoly_rbnode_struct * x;
    FLINT_TEST_INIT(state);

    flint_printf("rbtree....");
    fflush(stdout);

    mpoly_rbtree_init(t);

    for (i = 0; i < 1000; i++)
    {
        x = mpoly_rbtree_get(&new, t, (n_randint(state, 1000)));
        if (new)
            x->data = (void *) flint_malloc(8);
    }

    if (t->size == 0) {
        printf("FAIL\n");
        flint_printf("empty tree\n");
        flint_abort();
    }

    data = (void **) flint_malloc(t->size * sizeof(void *));
    keys = (slong *) flint_malloc(t->size * sizeof(slong));

    mpoly_rbtree_clear(t, data, keys);

    for (i = 0; i < t->size; i++)
    {
        if (i > 0 && keys[i-1] <= keys[i]) {
            printf("FAIL\n");
            flint_printf("keys out of order\n");
            flint_abort();
        }
        flint_free(data[i]);
    }

    flint_free(data);    
    flint_free(keys);

    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}




