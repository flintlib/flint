/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "impl.h"

void
_acb_mzeta_fixacc_init(_acb_mzeta_fixacc_struct * A, slong E)
{
#if FLINT_USES_PTHREAD
    A->lock = flint_malloc(sizeof(pthread_mutex_t));
    pthread_mutex_init(A->lock, NULL);
#else
    A->lock = NULL;
#endif
    fmpz_init(A->mid);
    fmpz_init(A->rad);
    A->E = E;
}

void
_acb_mzeta_fixacc_clear(_acb_mzeta_fixacc_struct * A)
{
#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(A->lock);
    flint_free(A->lock);
#endif
    fmpz_clear(A->mid);
    fmpz_clear(A->rad);
}

/* adds floor(x 2^E) (or of -x) to mid and ceil(rad(x) 2^E) + 1 to rad */
void
_acb_mzeta_fixacc_add(_acb_mzeta_fixacc_struct * A, const arb_t x, int negate)
{
    arf_t t;
    fmpz_t y, r;

    arf_init(t);
    fmpz_init(y);
    fmpz_init(r);

    arf_mul_2exp_si(t, arb_midref(x), A->E);
    if (negate)
        arf_neg(t, t);
    arf_get_fmpz(y, t, ARF_RND_FLOOR);
    arf_set_mag(t, arb_radref(x));
    arf_mul_2exp_si(t, t, A->E);
    arf_get_fmpz(r, t, ARF_RND_CEIL);
    fmpz_add_ui(r, r, 1);

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(A->lock);
#endif
    fmpz_add(A->mid, A->mid, y);
    fmpz_add(A->rad, A->rad, r);
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(A->lock);
#endif

    arf_clear(t);
    fmpz_clear(y);
    fmpz_clear(r);
}

void
_acb_mzeta_fixacc_get(arb_t res, const _acb_mzeta_fixacc_struct * A)
{
    mag_t m;
    mag_init(m);
    arb_set_fmpz(res, A->mid);
    mag_set_fmpz(m, A->rad);
    arb_add_error_mag(res, m);
    arb_mul_2exp_si(res, res, -A->E);
    mag_clear(m);
}
