/*
    Copyright (C) 2024 Fredrik Johansson
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "tune.h"

#undef FLINT_MPN_MULHIGH_K_TAB_SIZE
#define FLINT_MPN_MULHIGH_K_TAB_SIZE FLINT_MPN_MULHIGH_K_TAB_MAX_SIZE

#define flint_mpn_mulhigh_k_tab flint_mpn_mulhigh_k_tab_0
#define _flint_mpn_mulhigh_n_mulders _flint_mpn_mulhigh_n_mulders_0

FLINT_DLL extern short flint_mpn_mulhigh_k_tab[FLINT_MPN_MULHIGH_K_TAB_SIZE];
mp_limb_t _flint_mpn_mulhigh_n_mulders(mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);

/* IDEA:

    1. Set all entries in flint_mpn_mulhigh_k_tab to zeroes.
    2. Skip any further altering of those entries having a corresponding
       hardcoded version.
    3. When reasonable, check if basecase is faster.
    4. The sequence of k should be weakly increasing. We only need to consider k
       between n / 2 and n.
       TODO: How can we prune the search with this information?
    5. Run some warmups and run some hotlaps for each k. Only use the fastest
       consistent time (see tune.c).
    6. Use the k that corresponds to the fastest time and push that to
       flint_mpn_mulhigh_k_tab so that it can be used in consecutive runs.
    7. For large enough n, check if _flint_mpn_mulhigh_n_mul is faster. When we
       have 50 consecutive runs of full multiplication that are faster high
       multiplication, we exit.
*/

#define BASECASE_REASONABLE(n) ((n) < 200)

double measure_func(tune_func_t, void *, int, int);

void _tune_flint_mpn_mulhigh_n(void * vparam)
{
    struct mulhigh_param_0 * param = vparam;
    nn_ptr ap, bp, xp, yp;
    slong len;
    flint_time_t t0, t1;
    slong ix;

    ap = param->ap;
    bp = param->bp;
    xp = param->xp;
    yp = param->yp;
    len = param->len;

    flint_time_get(t0);
    for (ix = 0; ix < len; ix++)
        func(ap, bp, xp[ix], yp[ix]);
    flint_time_get(t1);

    return flint_time_nsec_diff(t1, t0);
}

void
tune_flint_mpn_mulhigh_n(int FLINT_UNUSED(warmups), int FLINT_UNUSED(min_runs))
{
    slong n;
    mp_ptr rp, xp, yp;
    flint_rand_t state;

    /* Initialize flint_mpn_mulhigh_k_tab */
    for (n = 0; n < FLINT_MPN_MULHIGH_K_TAB_MAX_SIZE; n++)
        flint_mpn_mulhigh_k_tab[n] = 0;

    rp = flint_malloc(2 * sizeof(mp_limb_t) * FLINT_MPN_MULHIGH_K_TAB_MAX_SIZE);
    xp = flint_malloc(sizeof(mp_limb_t) * FLINT_MPN_MULHIGH_K_TAB_MAX_SIZE);
    yp = flint_malloc(sizeof(mp_limb_t) * FLINT_MPN_MULHIGH_K_TAB_MAX_SIZE);

    flint_rand_init(state);
    flint_mpn_rrandom(xp, state, FLINT_MPN_MULHIGH_K_TAB_MAX_SIZE);
    flint_mpn_rrandom(yp, state, FLINT_MPN_MULHIGH_K_TAB_MAX_SIZE);
    flint_rand_clear(state);

    for (n = 1; n < FLINT_MPN_MULHIGH_K_TAB_MAX_SIZE; n++)
    {
        if (FLINT_HAVE_MULHIGH_FUNC(n))
            continue;

        if (BASECASE_REASONABLE(n))
        {

        }
        else
        {
        }
    }

    flint_free(rp);
    flint_free(xp);
    flint_free(yp);
}
