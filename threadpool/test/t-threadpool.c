/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "threadpool.h"
#include "fmpz.h"

threadpool_t global_thread_pool;

/******************************************************************************
    test1:
        calculate x = n! using a simple non-recursive residue scheme
        the master thread also participates in the calculation
*******************************************************************************/

typedef struct
{
    ulong modulus;
    ulong residue;
    ulong n;
    fmpz_t ans;
}
worker1_arg_struct;

void worker1(void * varg)
{
    worker1_arg_struct * arg = (worker1_arg_struct *) varg;
    ulong i;

    fmpz_one(arg->ans);
    for (i = arg->residue; i <= arg->n; i += arg->modulus)
    {
        fmpz_mul_ui(arg->ans, arg->ans, i);
    }
}

void test1(fmpz_t x, ulong n)
{
    ulong i, modulus;
    slong k, req, num_workers;
    worker1_arg_struct * args;
    threadpool_threadhandle * handles;

    /* request as many workers as possible */
    req = threadpool_get_size(global_thread_pool);
    handles = (threadpool_threadhandle *) flint_malloc(sizeof(threadpool_threadhandle) * req);
    num_workers = threadpool_request(global_thread_pool, handles, req);

    args = (worker1_arg_struct *) flint_malloc(sizeof(worker1_arg_struct) * FLINT_MAX(num_workers, 1));

    modulus = num_workers + 1;

    /* start up each of the available workers */
    for (k = 0; k < num_workers; k++)
    {
        args[k].residue = k + 1;
        args[k].modulus = modulus;
        args[k].n = n;
        fmpz_init(args[k].ans); /* worker expects an inited ans */
        threadpool_wake(global_thread_pool, handles[k], worker1, &args[k]);
    }

    /* do some work ourselves */
    fmpz_one(x);
    for (i = modulus; i <= n; i += modulus)
    {
        fmpz_mul_ui(x, x, i);        
    }

    /* wait for each of the workers and combine their answers */
    for (k = 0; k < num_workers; k++)
    {
        threadpool_wait(global_thread_pool, handles[k]);
        fmpz_mul(x, x, args[k].ans);
        fmpz_clear(args[k].ans); /* worker of course did not clean up its ans */
        threadpool_giveback(global_thread_pool, handles[k]);
    }

    flint_free(args);
    flint_free(handles);
}


/******************************************************************************
    test2 - calculate x = n! by recursively splitting the work when necessary
*******************************************************************************/

typedef struct
{
    ulong min;
    ulong max;
    fmpz_t ans;
}
worker2_arg_struct;

void test2_helper(fmpz_t x, ulong min, ulong max);

void worker2(void * varg)
{
    worker2_arg_struct * arg = (worker2_arg_struct *) varg;

    test2_helper(arg->ans, arg->min, arg->max);
}


/* set x = product of numbers in (min, max] */
void test2_helper(fmpz_t x, ulong min, ulong max)
{
    ulong i, mid;
    slong num_workers;
    threadpool_threadhandle handles[1];
    worker2_arg_struct args[1];

    FLINT_ASSERT(max >= min);

    if (max - min > UWORD(20)
        && ((num_workers = threadpool_request(global_thread_pool, handles, 1)) != 0))
    {
        FLINT_ASSERT(num_workers == 1);

        mid = min + ((max - min)/UWORD(2));

        /* start up the worker */
        args[0].min = min;
        args[0].max = mid;
        fmpz_init(args[0].ans);
        threadpool_wake(global_thread_pool, handles[0], worker2, &args[0]);

        /* do some work ourselves */
        test2_helper(x, mid, max);

        /* wait for the worker and combine its answer */
        threadpool_wait(global_thread_pool, handles[0]);
        fmpz_mul(x, x, args[0].ans);
        fmpz_clear(args[0].ans);
        threadpool_giveback(global_thread_pool, handles[0]);
    }
    else
    {
        fmpz_one(x);
        for (i = max; i > min; i--)
        {
            fmpz_mul_ui(x, x, i);
        }
    }
}

void test2(fmpz_t x, ulong n)
{
    test2_helper(x, 0, n);
}


int
main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("threadpool....");
    fflush(stdout);

    for (i = 0; i < 10*flint_test_multiplier(); i++)
    {
        fmpz_t x, y;

        fmpz_init(x);
        fmpz_init(y);
        threadpool_init(global_thread_pool, n_randint(state, 10));

        for (j = 0; j < 10; j++)
        {
            ulong n = n_randint(state, 1000);

            fmpz_fac_ui(y, n);

            test1(x, n);
            if (!fmpz_equal(x, y))
            {
                flint_printf("n: %wu\n", n);
                printf("x: "); fmpz_print(x); printf("\n");
                printf("y: "); fmpz_print(y); printf("\n");
                printf("test1 failed\n");
                flint_abort();
            }

            test2(x, n);
            if (!fmpz_equal(x, y))
            {
                flint_printf("n: %wu\n", n);
                printf("x: "); fmpz_print(x); printf("\n");
                printf("y: "); fmpz_print(y); printf("\n");
                printf("test2 failed\n");
                flint_abort();
            }
        }

        threadpool_clear(global_thread_pool);
        fmpz_clear(y);
        fmpz_clear(x);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

