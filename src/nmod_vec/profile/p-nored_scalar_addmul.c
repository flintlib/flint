/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Timings of _nmod_vec_nored_scalar_addmul_halflimb, res += c vec without
    reduction for entries and c below 2^(FLINT_BITS/2), which has vector
    code with AVX2 and NEON, against the plain loop compiled with the same
    flags. Times in nanoseconds per entry (the vectors stay in L1 or L2).
*/

#include "profiler.h"
#include "ulong_extras.h"
#include "nmod_vec.h"

static void
plain_loop(nn_ptr res, nn_srcptr vec, slong len, ulong c)
{
    slong i;

    for (i = 0; i < len; i++)
        res[i] += vec[i] * c;
}

static void
vector_code(nn_ptr res, nn_srcptr vec, slong len, ulong c)
{
    _nmod_vec_nored_scalar_addmul_halflimb(res, vec, len, c);
}

/* called through pointers, so that the repeated calls are not merged */
static void (* volatile fn[2])(nn_ptr, nn_srcptr, slong, ulong) =
    { vector_code, plain_loop };

int main(void)
{
    static const slong lens[] = { 8, 16, 24, 32, 64, 128, 256, 1024, 4096, 0 };
    slong il, len, i, reps;
    ulong c;
    nn_ptr res, vec;
    timeit_t timer;
    double t[2];
    int which, round;
    FLINT_TEST_INIT(state);

    flint_printf("res += c vec, entries and c below 2^%d, ns per entry\n",
                 FLINT_BITS / 2);
    flint_printf("   len     vector      plain\n");

    for (il = 0; (len = lens[il]) != 0; il++)
    {
        res = _nmod_vec_init(len);
        vec = _nmod_vec_init(len);
        for (i = 0; i < len; i++)
        {
            res[i] = n_randlimb(state);
            vec[i] = n_randlimb(state) >> (FLINT_BITS / 2);
        }
        c = n_randlimb(state) >> (FLINT_BITS / 2);

        reps = FLINT_MAX(1, 20000000 / (len + 10));
        t[0] = t[1] = -1;

        for (round = 0; round < 5; round++)
        {
            for (which = 0; which < 2; which++)
            {
                double tt;

                timeit_start_us(timer);
                for (i = 0; i < reps; i++)
                    fn[which](res, vec, len, c);
                timeit_stop_us(timer);

                tt = 1e3 * timer->wall / ((double) reps * len);
                if (t[which] < 0 || tt < t[which])
                    t[which] = tt;
            }
        }

        flint_printf("%6wd %10.3f %10.3f\n", len, t[0], t[1]);

        _nmod_vec_clear(res);
        _nmod_vec_clear(vec);
    }

    FLINT_TEST_CLEAR(state);
    flint_cleanup_master();
    return 0;
}
