/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"
#include "acb_hypgeom/impl.h"

/* n terms of the series with the tail bound; C0 is the bound of
   acb_hypgeom_pfq_bound_factor for n terms, or NULL */
static void
_acb_hypgeom_pfq_direct_bound(acb_t res, acb_srcptr a, slong p, acb_srcptr b, slong q,
    const acb_t z, slong n, const mag_t C0, slong prec)
{
    acb_t s, t;
    mag_t err, C;
    int terminates = 0;
    slong n0 = n;

    acb_init(s);
    acb_init(t);
    mag_init(err);
    mag_init(C);

    /* A series that terminates within n terms needs no tail bound.
       Otherwise the bound is computed first: if it is infinite, there
       is no point in summing the series, unless it terminates later
       (a polynomial outside of the unit disc, or truncated where its
       terms have not started to decrease for good): all its terms are
       summed then, up to the number of terms that is allowed for a
       series at this precision. */
    {
        slong i, last = WORD_MAX;

        for (i = 0; i < p; i++)
        {
            if (acb_is_int(a + i) && arb_is_nonpositive(acb_realref(a + i))
                && arf_cmpabs_2exp_si(arb_midref(acb_realref(a + i)), FLINT_BITS - 3) < 0)
            {
                slong m = 1 - arf_get_si(arb_midref(acb_realref(a + i)), ARF_RND_DOWN);
                last = FLINT_MIN(last, m);
            }
        }

        /* (a few more terms cost less than the bound) */
        if (last > n && last <= 50 + 10 * prec && last <= n + FLINT_MAX(32, n / 4))
            n = last;

        terminates = (last <= n);

        if (terminates)
        {
            mag_zero(C);
        }
        else
        {
            if (C0 != NULL && n == n0)
                mag_set(C, C0);
            else
                acb_hypgeom_pfq_bound_factor(C, a, p, b, q, z, n);

            if (mag_is_inf(C))
            {
                if (last <= 50 + 10 * prec)
                {
                    n = last;
                    terminates = 1;
                    mag_zero(C);
                }
                else
                {
                    acb_indeterminate(res);
                    acb_clear(s);
                    acb_clear(t);
                    mag_clear(err);
                    mag_clear(C);
                    return;
                }
            }
        }
    }

    acb_hypgeom_pfq_sum(s, t, a, p, b, q, z, n, prec);

    if (!acb_is_zero(t))
    {
        /* (not expected: the term after the last one of a terminating
           series is zero) */
        if (terminates)
            acb_hypgeom_pfq_bound_factor(C, a, p, b, q, z, n);

        acb_get_mag(err, t);
        mag_mul(err, err, C);

        if (_acb_vec_is_real(a, p) && _acb_vec_is_real(b, q) && acb_is_real(z))
            arb_add_error_mag(acb_realref(s), err);
        else
            acb_add_error_mag(s, err);
    }

    acb_set(res, s);

    acb_clear(s);
    acb_clear(t);
    mag_clear(err);
    mag_clear(C);
}

/* The series with the number of terms of acb_hypgeom_pfq_choose_n_max,
   and with the tail bound that was computed for the choice. */
void
_acb_hypgeom_pfq_direct_max(acb_t res, acb_srcptr a, slong p, acb_srcptr b, slong q,
    const acb_t z, slong prec, slong n_max)
{
    mag_t C;
    slong n;
    int have_bound;

    mag_init(C);
    n = _acb_hypgeom_pfq_choose_n_bound(C, &have_bound, a, p, b, q, z, prec, n_max);
    _acb_hypgeom_pfq_direct_bound(res, a, p, b, q, z, n, have_bound ? C : NULL, prec);
    mag_clear(C);
}

void
acb_hypgeom_pfq_direct(acb_t res, acb_srcptr a, slong p, acb_srcptr b, slong q,
    const acb_t z, slong n, slong prec)
{
    if (n < 0)
        _acb_hypgeom_pfq_direct_max(res, a, p, b, q, z, prec, FLINT_MIN(WORD_MAX / 2, 50 + 10 * prec));
    else
        _acb_hypgeom_pfq_direct_bound(res, a, p, b, q, z, n, NULL, prec);
}
