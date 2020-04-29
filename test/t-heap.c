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
#include "heap.h"
#include "ulong_extras.h"

void check_heap(heap_t h, slong *vals)
{
    slong i, j;
    for (i = h->num - 1; i >= 0; --i)
    {
        if (h->pos[h->idx[i]] != i)
        {
            flint_printf("Mismatch between idx and pos\n");
            abort();
        }
        if (h->val[h->idx[i]] != vals[h->idx[i]])
        {
            flint_printf("Mismatch between inside and outside scores\n");
            abort();
        }
        if (i==0) break;
        j = (i - 1)/2;
        if (h->val[h->idx[i]] < h->val[h->idx[j]])
        {
            flint_printf("Not a heap: %2wd:%2wd:%2wd < %2wd:%2wd:%2wd\n",
            i, h->idx[i], h->val[h->idx[i]], j, h->idx[j], h->val[h->idx[j]]);
            abort();
        }
    }
}

int main(void)
{
    slong rep, num, *vals, nreps = 1000, i, j, val, oval;
    heap_t h;
    FLINT_TEST_INIT(state);

    flint_printf("heap....");
    fflush(stdout);
    for (rep = 0; rep < nreps; rep++)
    {
        num = n_randint(state, 1000);
        vals = flint_malloc(num*sizeof(*vals));

        heap_init(h, num);
        for (i = 0; i < num; ++i) 
        {
            vals[i] = n_randtest(state);
            heap_push(h, vals[i]);
        }
        if (h->num != num)
        {
            flint_printf("FAIL: Added %wd values, heap only has %wd\n", num, h->num);
            abort();
        }
        check_heap(h, vals);

        /* Modify each value and check if heap maintained */
        for (i = 0; i < num; ++i)
        {
            vals[i] = n_randtest(state);
            heap_adjust(h, i, vals[i]);
            check_heap(h, vals);
        }

        /* Pop items from heap and check that in sorted order */
        for (i = 0; i < num; ++i)
        {
            heap_pop(h, &val);
            if (i > 0 && val < oval)
            {
                flint_printf("FAIL: popping from heap out of order\n");
                abort();
            }
            check_heap(h, vals);
            oval = val;
        }
        flint_free(vals);
        heap_clear(h);
    }
    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
