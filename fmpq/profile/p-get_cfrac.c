/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "profiler.h"

int
main(void)
{
    int i;

    flint_printf("\n");
    fflush(stdout);

    for (i = 0; i <= 24; i++)
    {
        fmpq_t x, r, l, m;
        fmpz *c1, *c2;
        slong n1, /*n2, */bound;
        slong gcd_time;
        timeit_t timer;
        slong expected_length;

        fmpq_init(x);
        fmpq_init(r);
        fmpq_init(l);
        fmpq_init(m);

        fmpz_fib_ui(fmpq_numref(x), (1 << i) + 1);
        fmpz_fib_ui(fmpq_denref(x), (1 << i));

flint_printf("\n--- fib(1+2^%wd)/fib(2^%wd) (numerator bits = %wu) ---\n", i, i, fmpz_bits(fmpq_numref(x)));

        expected_length = FLINT_MAX(1, (1 << i) - 1);

        timeit_start(timer);
        fmpq_canonicalise(x);
        timeit_stop(timer);
        gcd_time = timer->wall;
        flint_printf("gcd: %wd\n", timer->wall);
        gcd_time = FLINT_MAX(1, gcd_time);


        bound = fmpq_cfrac_bound(x);
        if (bound < expected_length)
        {
            flint_printf("FAIL0: bound = %wd, expected = %wd\n", bound, expected_length);
            flint_abort();
        }

        c1 = _fmpz_vec_init(bound);
        c2 = _fmpz_vec_init(bound);

        timeit_start(timer);
        n1 = fmpq_get_cfrac(c1, r, x, bound);
        timeit_stop(timer);
        if (!fmpq_is_zero(r) || n1 != expected_length)
        {
            flint_printf("FAIL1: n1 = %wd, expected = %wd\n", n1, expected_length);
            flint_abort();
        }
        flint_printf("new: %wd  (new/gcd: %f)\n", timer->wall, (double)(timer->wall)/(double)(gcd_time));

        fmpz_fib_ui(fmpq_numref(r), (1 << i) + 1);
        fmpz_fib_ui(fmpq_denref(r), (1 << i));
        fmpz_fib_ui(fmpq_numref(l), (1 << i) + 2);
        fmpz_fib_ui(fmpq_denref(l), (1 << i) + 1);
        if (i == 0)
            fmpq_swap(l, r);

        timeit_start(timer);
        _fmpq_simplest_between(fmpq_numref(m), fmpq_denref(m),
                               fmpq_numref(l), fmpq_denref(l),
                               fmpq_numref(r), fmpq_denref(r));
        timeit_stop(timer);
        if (!fmpq_equal(m, x))
        {
            flint_printf("FAIL: between\n");
            flint_abort();
        }
        flint_printf("bet: %wd  (bet/gcd: %f)\n", timer->wall, (double)(timer->wall)/(double)(gcd_time));

#if 0
        timeit_start(timer);
        n2 = fmpq_get_cfrac_naive(c2, r, x, bound);
        timeit_stop(timer);
        if (!fmpq_is_zero(r) || n2 != expected_length)
        {
            flint_printf("FAIL2: n2 = %wd, expected = %wd\n", n2, expected_length);
            flint_abort();
        }

        flint_printf("old: %wd\n", timer->wall);
#endif
        _fmpz_vec_clear(c1, bound);
        _fmpz_vec_clear(c2, bound);
        fmpq_clear(x);
        fmpq_clear(r);
        fmpq_clear(l);
        fmpq_clear(m);
    }

    flint_cleanup_master();
    return 0;
}

