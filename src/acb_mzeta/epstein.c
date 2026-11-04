/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"
#include "impl.h"

int
acb_mzeta_epstein_generic(acb_t res, const acb_t nu, const arb_mat_t A,
    arb_srcptr x, arb_srcptr y, int flags, acb_mzeta_cache_t cache, slong prec)
{
    acb_mzeta_cache_t tmp_cache;
    int ok;

    if (cache == NULL)
    {
        acb_mzeta_cache_init(tmp_cache);
        ok = acb_mzeta_epstein_generic(res, nu, A, x, y, flags, tmp_cache, prec);
        acb_mzeta_cache_clear(tmp_cache);
        return ok;
    }

    /* the shell algorithm (exact rational input) above double precision */
    if (prec > 95 && acb_mzeta_epstein_shell(res, nu, A, x, y, flags, cache, prec))
        return 1;

    /* double-precision kernels at low precision (unless the result lost
       too much accuracy, e.g. by cancellation) */
    if (prec <= 95 && acb_mzeta_epstein_crandall_dfloat(res, nu, A, x, y, flags, cache, prec)
        && acb_rel_accuracy_bits(res) >= prec - 12)
        return 1;

    /* theta-Mellin at high precision, where the number of lattice points
       of the Crandall sums (~ prec^(d/2)) dominates (measured crossovers) */
    if (((arb_mat_nrows(A) == 3 && prec >= 800) || (arb_mat_nrows(A) == 4 && prec >= 400)
            || (arb_mat_nrows(A) >= 5 && prec >= 250))
        && arb_mat_nrows(A) <= 8
        && acb_mzeta_epstein_theta(res, nu, A, x, y, flags, cache, prec) && acb_is_finite(res))
        return 1;

    return acb_mzeta_epstein_crandall(res, nu, A, x, y, flags, cache, prec);
}

/* the user-facing function: a thread-local cache keeps the data of the
   most recently used lattices (and tables), for the common use of many
   evaluations with a fixed lattice */
FLINT_TLS_PREFIX acb_mzeta_cache_struct _acb_mzeta_tls_cache;
FLINT_TLS_PREFIX int _acb_mzeta_tls_cache_initialized = 0;

static void
_acb_mzeta_tls_cache_cleanup(void)
{
    if (_acb_mzeta_tls_cache_initialized)
    {
        acb_mzeta_cache_clear(&_acb_mzeta_tls_cache);
        _acb_mzeta_tls_cache_initialized = 0;
    }
}

int
acb_mzeta_epstein(acb_t res, const acb_t nu, const arb_mat_t A,
    arb_srcptr x, arb_srcptr y, int flags, slong prec)
{
    if (!_acb_mzeta_tls_cache_initialized)
    {
        acb_mzeta_cache_init(&_acb_mzeta_tls_cache);
        _acb_mzeta_tls_cache_initialized = 1;
        flint_register_cleanup_function(_acb_mzeta_tls_cache_cleanup);
    }
    return acb_mzeta_epstein_generic(res, nu, A, x, y, flags, &_acb_mzeta_tls_cache, prec);
}
