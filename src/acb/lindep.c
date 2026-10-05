/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpz_lll.h"
#include "acb.h"

/* The precision at which the combinations of vec with coefficients of
   at most coeff_bits bits are evaluated for validation: exact (the
   midpoints are then exact sums and only the radii of vec contribute)
   when the midpoints span a moderate range of exponents, otherwise
   enough bits beyond the working precision that rounding errors stay
   below the resolution of the lattice. */
static slong
_acb_lindep_check_prec(acb_srcptr vec, slong len, slong coeff_bits, slong prec)
{
    fmpz_t emin, emax, t;
    slong i, j, wp;
    int any = 0;

    fmpz_init(emin);
    fmpz_init(emax);
    fmpz_init(t);

    for (i = 0; i < len; i++)
    {
        for (j = 0; j < 2; j++)
        {
            const arf_struct * x = (j == 0) ? arb_midref(acb_realref(vec + i)) : arb_midref(acb_imagref(vec + i));

            if (arf_is_zero(x))
                continue;

            /* the exponent of the least significant bit */
            fmpz_sub_ui(t, ARF_EXPREF(x), arf_bits(x));
            if (!any || fmpz_cmp(t, emin) < 0)
                fmpz_set(emin, t);
            if (!any || fmpz_cmp(ARF_EXPREF(x), emax) > 0)
                fmpz_set(emax, ARF_EXPREF(x));
            any = 1;
        }
    }

    wp = prec + 10 + coeff_bits + FLINT_BIT_COUNT(len);

    if (any)
    {
        fmpz_sub(t, emax, emin);
        if (fmpz_cmp_si(t, 16 * FLINT_MIN(prec, WORD(1) << 24) + 4096) <= 0)
            wp = ARF_PREC_EXACT;
    }

    fmpz_clear(emin);
    fmpz_clear(emax);
    fmpz_clear(t);

    return wp;
}

void
_acb_lindep_combination(acb_t s, acb_srcptr vec, const fmpz * c, slong len, slong prec)
{
    acb_t t;
    slong j, wp;

    wp = _acb_lindep_check_prec(vec, len, FLINT_ABS(_fmpz_vec_max_bits(c, len)), prec);

    acb_init(t);
    acb_zero(s);
    for (j = 0; j < len; j++)
    {
        acb_mul_fmpz(t, vec + j, c + j, wp);
        acb_add(s, s, t, wp);
    }
    acb_clear(t);
}

slong
acb_lindep(fmpz_mat_t rel, acb_srcptr vec, slong len, slong prec)
{
    fmpz_mat_t A;
    fmpz_lll_t ctx;
    fmpz_t scale_exp;
    arf_t tmpr, halfr;
    mag_t max_size, max_rad, tmpmag;
    acb_t s;
    slong i, j, num = 0, accuracy;
    int nonreal;

    fmpz_mat_zero(rel);

    if (len < 1 || !_acb_vec_is_finite(vec, len))
        return 0;

    nonreal = 0;
    for (i = 0; i < len; i++)
        if (!arb_contains_zero(acb_imagref(vec + i)))
            nonreal = 1;

    mag_init(max_size);
    mag_init(max_rad);
    mag_init(tmpmag);

    for (i = 0; i < len; i++)
    {
        arf_get_mag(tmpmag, arb_midref(acb_realref(vec + i)));
        mag_max(max_size, max_size, tmpmag);
        arf_get_mag(tmpmag, arb_midref(acb_imagref(vec + i)));
        mag_max(max_size, max_size, tmpmag);
        mag_max(max_rad, max_rad, arb_radref(acb_realref(vec + i)));
        mag_max(max_rad, max_rad, arb_radref(acb_imagref(vec + i)));
    }

    if (mag_is_zero(max_size))
    {
        /* every entry has midpoint zero: the unit vectors */
        for (i = 0; i < len; i++)
            fmpz_one(fmpz_mat_entry(rel, i, i));
        mag_clear(max_size);
        mag_clear(max_rad);
        mag_clear(tmpmag);
        return len;
    }

    fmpz_mat_init(A, len, len + 1 + nonreal);
    fmpz_init(scale_exp);
    arf_init(tmpr);
    arf_init(halfr);
    acb_init(s);
    arf_set_d(halfr, 0.5);

    prec = FLINT_MAX(prec, 2);
    if (!mag_is_zero(max_rad))
    {
        accuracy = _fmpz_sub_small(MAG_EXPREF(max_size), MAG_EXPREF(max_rad));
        accuracy = FLINT_MAX(accuracy, 10);
        prec = FLINT_MIN(prec, accuracy);
    }

    fmpz_neg(scale_exp, MAG_EXPREF(max_size));
    fmpz_add_ui(scale_exp, scale_exp, prec);
    /* some of the bits are kept for checking (protection against
       spurious relations) */
    fmpz_sub_ui(scale_exp, scale_exp, FLINT_MAX(10, prec / 20));

    for (i = 0; i < len; i++)
        fmpz_one(fmpz_mat_entry(A, i, i));

    for (i = 0; i < len; i++)
    {
        arf_mul_2exp_fmpz(tmpr, arb_midref(acb_realref(vec + i)), scale_exp);
        arf_add(tmpr, tmpr, halfr, prec, ARF_RND_NEAR);
        arf_floor(tmpr, tmpr);
        arf_get_fmpz(fmpz_mat_entry(A, i, len), tmpr, ARF_RND_NEAR);

        if (nonreal)
        {
            arf_mul_2exp_fmpz(tmpr, arb_midref(acb_imagref(vec + i)), scale_exp);
            arf_add(tmpr, tmpr, halfr, prec, ARF_RND_NEAR);
            arf_floor(tmpr, tmpr);
            arf_get_fmpz(fmpz_mat_entry(A, i, len + 1), tmpr, ARF_RND_NEAR);
        }
    }

    fmpz_lll_context_init(ctx, 0.75, 0.51, 1, 0);
    fmpz_lll(A, NULL, ctx);

    /* the rows of the reduced basis whose combinations are numerically
       zero, in order */
    for (i = 0; i < len; i++)
    {
        if (_fmpz_vec_is_zero(fmpz_mat_row(A, i), len))
            continue;

        _acb_lindep_combination(s, vec, fmpz_mat_row(A, i), len, prec);

        if (acb_contains_zero(s))
        {
            for (j = 0; j < len; j++)
                fmpz_set(fmpz_mat_entry(rel, num, j), fmpz_mat_entry(A, i, j));
            num++;
        }
    }

    fmpz_mat_clear(A);
    fmpz_clear(scale_exp);
    arf_clear(tmpr);
    arf_clear(halfr);
    acb_clear(s);
    mag_clear(max_size);
    mag_clear(max_rad);
    mag_clear(tmpmag);

    return num;
}
