/*
    Copyright (C) 2009 William Hart

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
#include "hashmap.h"
#include "ulong_extras.h"

int main(void)
{
    slong rep, num, *keys, *vals, nreps = 1000, i, *ptr;
    hashmap_t h;
    FLINT_TEST_INIT(state);

    flint_printf("hashmap....");
    fflush(stdout);
    for (rep = 0; rep < nreps; rep++)
    {
        num = n_randint(state, 1000);
        keys = flint_malloc(num*sizeof(*keys));
        vals = flint_malloc(num*sizeof(*vals));

        hashmap_init(h, num);
        for (i = 0; i < num; ++i) 
        {
            do keys[i] = n_randtest(state);
            while (hashmap_get(h, keys[i]) != NULL);
            vals[i] = n_randtest(state);
            hashmap_put(h, keys[i], &vals[i]);
        }
        if (h->num != num)
        {
            flint_printf("FAIL: Added %wd keys, hashmap only has %wd\n", num, h->num);
            abort();
        }

        for (i = 0; i < num; ++i)
        {
            if (h->keys[i] != keys[i])
            {
                flint_printf("FAIL: Mismatch key at %wd:%wd\n", i, keys[i]);
                abort();
            }
            ptr = hashmap_get(h, keys[i]);
            if (ptr != &vals[i])
            {
                flint_printf("FAIL: Mismatch value at %wd:%wd\n", i, keys[i]);
                abort();
            }
        }
        for (i = 0; i < num; ++i)
        {
            hashmap_rem(h, keys[i]);
            ptr = hashmap_get(h, keys[i]);
            if (ptr != NULL)
            {
                flint_printf("FAIL: removed key %wd, still exists in hashtable\n");
                abort();
            }
        }
        flint_free(keys);
        flint_free(vals);
        hashmap_clear(h);
    }
    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
