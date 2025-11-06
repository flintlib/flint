/*
    Copyright (C) 2017 Apoorv Mishra
    Copyright (C) 2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "test_helpers.h"
#include "ulong_extras.h"

#define N               8
#define NR_SAMPLES      N * (1 << 8)
#define NR_SAMPLES_MAX  NR_SAMPLES
#define LIMIT_SMALL     (2 + n_randint(state, 1 << 14))
#define LIMIT_BIG \
    (UWORD_MAX / 2 + n_randint(state, UWORD(1) << FLINT_BITS / 2))
#define SQR(x)          ((x) * (x))
#define ABS(x)          ((x) > 0 ? (x) : -(x))
#define MU(limit)       (((double) (limit) - 1) / 2)
#define VAR(limit)      ((SQR((double) (limit)) - 1) / 12)
#if FLINT64
# define FMT_LIM        "0x%016" _WORD_FMT "x"
#else
# define FMT_LIM        "0x%08" _WORD_FMT "x"
#endif

static inline int check_mu(double mu, ulong lim, slong nr_samples)
{
    double z = (mu - MU(lim)) / sqrt(VAR(lim) / nr_samples);
    return fabs(z) <= 4.475328424654207; /* alpha = 2^{-16} */
}

static inline int check_var(double var, ulong lim, slong nr_samples)
{
    /* Approximate chi-squared through normal distribution */
    double chi2 = nr_samples * var / VAR(lim);
    double mu_chi = nr_samples - 1;
    double z_chi = (chi2 - mu_chi) / sqrt(2 * mu_chi);
    return fabs(z_chi) <= 4.475328424654207; /* alpha = 2^{-16} */
}

TEST_FUNCTION_START(n_urandint, state)
{
    int result;
    ulong * ts = flint_malloc(NR_SAMPLES_MAX * sizeof(ulong));

    for (slong ix = 0; ix < 1000 * flint_test_multiplier(); ix++)
    {
            int is_big = n_randint(state, 2) == 0;
        ulong lim = is_big ? LIMIT_BIG : LIMIT_SMALL;
        slong nr_samples = NR_SAMPLES;
        double mu, var = 0;
        double up[N] = {0.0};

        if (is_big)
        {
            ulong t0 = 0, t1 = 0;
            for (slong ix = 0; ix < nr_samples; ix++)
            {
                ts[ix] = n_urandint(state, lim);
                add_ssaaaa(t1, t0, t1, t0, 0, ts[ix]);
            }
#if FLINT64
            mu = (t0 + 0x1.0p64 * t1) / nr_samples;
#else
            mu = (t0 + 0x1.0p32 * t1) / nr_samples;
#endif
        }
        else
        {
            ulong t0 = 0;
            for (slong ix = 0; ix < nr_samples; ix++)
            {
                ts[ix] = n_urandint(state, lim);
                t0 += ts[ix];
            }
            mu = (double) t0 / nr_samples;
        }

        for (slong ix = 0; ix < nr_samples / N; ix++)
            for (slong jx = 0; jx < N; jx++)
                up[jx] += SQR(ts[N * ix + jx] - mu);

        for (slong jx = 0; jx < N; jx++) var += up[jx];
        var /= nr_samples;

        result = check_mu(mu, lim, nr_samples)
            && check_var(var, lim, nr_samples);
        if (!result)
            TEST_FUNCTION_FAIL("ix = %wd\n"
                               "lim = " FMT_LIM "\n"
                               "nr_samples = %wd\n"
                               "mu fail: %d, var fail: %d\n"
                               "MU = %10e, VAR = %10e\n"
                               "mu = %10e, var = %10e\n",
                               ix, lim, nr_samples,
                               !check_mu(mu, lim, nr_samples),
                               !check_var(var, lim, nr_samples),
                               MU(lim), VAR(lim), mu, var);
    }
    flint_free(ts);

    TEST_FUNCTION_END(state);
}

#undef N
#undef NR_SAMPLES
#undef NR_SAMPLES_MAX
#undef LIMIT_SMALL
#undef LIMIT_BIG
#undef SQR
#undef ABS
#undef MU
#undef VAR
#undef FMT_LIM
