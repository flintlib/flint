/*
    Copyright (C) 2026 KainanP

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
#include "fmpz.h"

/*
    Check that the free list of mpz's survives memory functions which call
    _fmpz_clear_mpz before they return, as a garbage collector running
    finalizers does (e.g. Julia's allocator in Nemo).
*/

#define NUM_TRACKED 64
#define NUM_PENDING 100

static void * (* sys_alloc)(size_t);
static void * (* sys_calloc)(size_t, size_t);
static void * (* sys_realloc)(void *, size_t);
static void (* sys_free)(void *);

/* "Finalizers": integers cleared by the next allocation */
static fmpz * pending;
static slong num_pending;
static int in_finalizers;

/* Blocks returned by gc_realloc, with their requested sizes */
static void * tracked_ptr[NUM_TRACKED];
static size_t tracked_size[NUM_TRACKED];
static slong num_tracked;

/* Blocks released by gc_realloc. They stay allocated until the end of the
   test so that stale writes to them do not corrupt the heap. */
static void * released[NUM_TRACKED];
static slong num_released;

static void
run_finalizers(void)
{
    if (in_finalizers)
        return;

    in_finalizers = 1;
    while (num_pending > 0)
    {
        num_pending--;
        fmpz_clear(pending + num_pending);
    }
    in_finalizers = 0;
}

static void *
gc_alloc(size_t n)
{
    void * p = sys_alloc(n);
    run_finalizers();
    return p;
}

static void *
gc_calloc(size_t n, size_t m)
{
    void * p = sys_calloc(n, m);
    run_finalizers();
    return p;
}

/* Moving realloc: the new block is zeroed beyond the old contents, and has
   room for twice the requested size so that a stale write past the end of
   the old size stays within the block. */
static void *
gc_realloc(void * p, size_t n)
{
    slong i = 0;
    void * q;

    if (p != NULL)
        for (i = 0; i < num_tracked && tracked_ptr[i] != p; i++) ;

    if ((p != NULL && i == num_tracked)
        || (p == NULL && num_tracked == NUM_TRACKED)
        || (p != NULL && num_released == NUM_TRACKED))
    {
        q = sys_realloc(p, n);
    }
    else
    {
        q = sys_calloc(2 * n, 1);

        if (p == NULL)
        {
            i = num_tracked++;
        }
        else
        {
            memcpy(q, p, FLINT_MIN(n, tracked_size[i]));
            released[num_released++] = p;
        }

        tracked_ptr[i] = q;
        tracked_size[i] = n;
    }

    run_finalizers();
    return q;
}

static void
gc_free(void * p)
{
    slong i;

    for (i = 0; i < num_tracked && tracked_ptr[i] != p; i++) ;

    if (i < num_tracked)
    {
        tracked_ptr[i] = tracked_ptr[num_tracked - 1];
        tracked_size[i] = tracked_size[num_tracked - 1];
        num_tracked--;
    }

    sys_free(p);
}

TEST_FUNCTION_START(fmpz_clear_mpz_reentrant, state)
{
#if FLINT_USES_GC
    TEST_FUNCTION_END_SKIPPED(state);
#else
    slong iter;

    __flint_get_memory_functions(&sys_alloc, &sys_calloc, &sys_realloc, &sys_free);

    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong i, n;
        fmpz * a;
        mpz_ptr * z;

        n = 4000 + n_randint(state, 12000);

        a = flint_malloc(sizeof(fmpz) * n);
        pending = flint_malloc(sizeof(fmpz) * NUM_PENDING);
        z = flint_malloc(sizeof(mpz_ptr) * (n + NUM_PENDING));

        /* start with an empty free list */
        _fmpz_cleanup();

        __flint_set_memory_functions(gc_alloc, gc_calloc, gc_realloc, gc_free);

        for (i = 0; i < n; i++)
        {
            fmpz_init(a + i);
            fmpz_randbits(a + i, state, 2 * FLINT_BITS);
            fmpz_setbit(a + i, 2 * FLINT_BITS);
        }

        for (i = 0; i < NUM_PENDING; i++)
        {
            fmpz_init(pending + i);
            fmpz_randbits(pending + i, state, 2 * FLINT_BITS);
            fmpz_setbit(pending + i, 2 * FLINT_BITS);
        }

        /* the free list grows while these are cleared, so the integers in
           pending get cleared from inside the memory functions */
        num_pending = NUM_PENDING;
        for (i = 0; i < n; i++)
            fmpz_clear(a + i);

        /* in case the free list did not grow (large pages) */
        run_finalizers();

        /* every mpz on the free list must be valid */
        for (i = 0; i < n + NUM_PENDING; i++)
        {
            z[i] = _fmpz_new_mpz();

            if (z[i] == NULL)
                TEST_FUNCTION_FAIL("free list entry %wd is NULL (n = %wd)\n", i, n);
        }

        for (i = 0; i < n + NUM_PENDING; i++)
            _fmpz_clear_mpz(PTR_TO_COEFF(z[i]));

        _fmpz_cleanup();

        __flint_set_memory_functions(sys_alloc, sys_calloc, sys_realloc, sys_free);

        for (i = 0; i < num_released; i++)
            sys_free(released[i]);
        num_released = 0;
        num_tracked = 0;

        flint_free(a);
        flint_free(pending);
        flint_free(z);
    }

    TEST_FUNCTION_END(state);
#endif
}
