/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "thread_support.h"

typedef struct
{
    slong idx;
    slong num;
    slong len;
    fmpz * vec;
}
worker_arg_struct;

void worker1(void * varg)
{
    worker_arg_struct * arg = (worker_arg_struct *) varg;
    slong idx = arg->idx;
    slong num = arg->num;
    slong len = arg->len;
    fmpz * vec = arg->vec;
    slong i;
    slong start = (idx + 0)*len/num;
    slong stop = (idx + 1)*len/num;

    for (i = start; i < stop; i++)
    {
        fmpz_set_ui(vec + i, i);
        fmpz_mul_2exp(vec + i, vec + i, 100);
    }
}

void worker2(void * varg)
{
    worker_arg_struct * arg = (worker_arg_struct *) varg;
    slong idx = arg->idx;
    slong num = arg->num;
    slong len = arg->len;
    fmpz * vec = arg->vec;
    slong i;

    for (i = idx + num; i < len; i += num)
    {
        fmpz_add(vec + idx, vec + idx, vec + i);
        fmpz_zero(vec + i);
    }
}

void worker3(void * varg)
{
    worker_arg_struct * arg = (worker_arg_struct *) varg;
    slong idx = arg->idx;
    slong num = arg->num;
    slong len = arg->len;
    fmpz * vec = arg->vec;
    slong i;

    for (i = idx + num; i < len; i += num)
    {
        fmpz_set_ui(vec + i, i);
        fmpz_mul_2exp(vec + i, vec + i, 110);
    }
}

TEST_FUNCTION_START(fmpz_stress, state)
{
    slong i, j, k;
    slong max_num_threads = 5;
    thread_pool_handle * handles;
    slong num_handles;
    worker_arg_struct * wargs;

    wargs = (worker_arg_struct *) flint_malloc(max_num_threads*
                                                    sizeof(worker_arg_struct));

    /* checking parallel clearing of fmpz's allocated from the same thread */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        slong n;
        fmpz * v;
        fmpz_t total, check;

        fmpz_init(total);
        fmpz_init(check);

        n = n_randint(state, 10000) + 10000;
        v = _fmpz_vec_init(n);

        for (j = 0; j < 4; j++)
        {
            flint_set_num_threads(n_randint(state, max_num_threads) + 1);
            num_handles = flint_request_threads(&handles, max_num_threads);

            for (k = 0; k <= num_handles; k++)
            {
                wargs[k].idx = k;
                wargs[k].num = num_handles + 1;
                wargs[k].len = n;
                wargs[k].vec = v;
            }

            /* fill v with blocks of large fmpz's from the same thread */
            for (k = 0; k < num_handles; k++)
                thread_pool_wake(global_thread_pool, handles[k], 0, worker1, &wargs[k]);
            worker1(&wargs[num_handles]);
            for (k = 0; k < num_handles; k++)
                thread_pool_wait(global_thread_pool, handles[k]);

            /* rip up the blocks in succession, zeroing as we go */
            for (k = 0; k < num_handles; k++)
                thread_pool_wake(global_thread_pool, handles[k], 0, worker2, &wargs[k]);
            worker2(&wargs[num_handles]);
            for (k = 0; k < num_handles; k++)
                thread_pool_wait(global_thread_pool, handles[k]);

            fmpz_zero(total);
            for (k = 0; k <= num_handles; k++)
                fmpz_add(total, total, v + k);

            /* put large fmpz's back in succession */
            for (k = 0; k < num_handles; k++)
                thread_pool_wake(global_thread_pool, handles[k], 0, worker3, &wargs[k]);
            worker3(&wargs[num_handles]);
            for (k = 0; k < num_handles; k++)
                thread_pool_wait(global_thread_pool, handles[k]);

            flint_give_back_threads(handles, num_handles);

            fmpz_set_ui(check, n);
            fmpz_mul_ui(check, check, n - 1);
            fmpz_mul_2exp(check, check, 99);
            if (!fmpz_equal(total, check) || !_fmpz_is_canonical(total))
            {
                flint_printf("FAIL:\n");
                flint_printf("total: "); fmpz_print(total); flint_printf("\n");
                flint_printf("check: "); fmpz_print(check); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        _fmpz_vec_clear(v, n);

        fmpz_clear(total);
        fmpz_clear(check);
    }

    flint_free(wargs);

    TEST_FUNCTION_END(state);
}
