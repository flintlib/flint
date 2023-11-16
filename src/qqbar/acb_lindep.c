/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpz_lll.h"
#include "qqbar.h"

int
_qqbar_acb_lindep(fmpz * rel, acb_srcptr vec, slong len, int check, slong prec)
{
    arf_t tmpr, halfr;
    fmpz_mat_t A;
    fmpz_lll_t ctx;
    fmpz_t scale_exp;
    int nonreal, found;
    slong i, accuracy;
    acb_t z2;
    mag_t max_size, max_rad, tmpmag;

    for (i = 0; i < len; i++)
        if (!acb_is_finite(vec + i))
            return 0;

    found = 0;

    nonreal = 0;
    for (i = 0; i < len; i++)
        if (!arb_contains_zero(acb_imagref(vec + i)))
            nonreal = 1;

    fmpz_mat_init(A, len, len + 1 + nonreal);
    fmpz_init(scale_exp);
    acb_init(z2);
    arf_init(tmpr);
    arf_init(halfr);
    mag_init(max_size);
    mag_init(max_rad);
    mag_init(tmpmag);

    arf_set_d(halfr, 0.5);

    for (i = 0; i < len; i++)
    {
        arf_get_mag(tmpmag, arb_midref(acb_realref(vec + i)));
        mag_max(max_size, max_size, tmpmag);
        arf_get_mag(tmpmag, arb_midref(acb_imagref(vec + i)));
        mag_max(max_size, max_size, tmpmag);

        mag_max(max_rad, max_rad, arb_radref(acb_realref(vec + i)));
        mag_max(max_rad, max_rad, arb_radref(acb_imagref(vec + i)));
    }

    prec = FLINT_MAX(prec, 2);
    if (!mag_is_zero(max_size) && !mag_is_zero(max_rad))
    {
        accuracy = _fmpz_sub_small(MAG_EXPREF(max_size), MAG_EXPREF(max_rad));
        accuracy = FLINT_MAX(accuracy, 10);
        prec = FLINT_MIN(prec, accuracy);
    }

    if (mag_is_zero(max_size))
    {
        fmpz_zero(scale_exp);  /* todo: quick return? */
    }
    else
    {
        fmpz_neg(scale_exp, MAG_EXPREF(max_size));
        fmpz_add_ui(scale_exp, scale_exp, prec);
    }

    /* Using 5% of the bits for checking will provide some protection
       against spurious relations */
    fmpz_sub_ui(scale_exp, scale_exp, FLINT_MAX(10, prec * 0.05));

    /* Create matrix */
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

    /* LLL reduction */
    fmpz_lll_context_init(ctx, 0.75, 0.51, 1, 0);
    fmpz_lll(A, NULL, ctx);

    for (i = 0; i < len; i++)
        fmpz_set(rel + i, fmpz_mat_entry(A, 0, i));

#if 0
    flint_printf("rel ");
    for (i = 0; i < len; i++)
    {
        fmpz_print(rel + i);
        flint_printf(" ");
    }
    flint_printf("\n");
#endif

    /* Heuristic check */
    if (check)
    {
        for (i = 0; i < len; i++)
            acb_addmul_fmpz(z2, vec + i, rel + i, prec + 10);
        found = !_fmpz_vec_is_zero(rel, len) && acb_contains_zero(z2);
    }
    else
    {
        found = !_fmpz_vec_is_zero(rel, len);
    }

/*
    if (found)
    {
        flint_printf("REL:\n");
        fmpz_mat_print_pretty(A); flint_printf("\n\n");
    }
*/

    fmpz_mat_clear(A);
    fmpz_clear(scale_exp);
    acb_clear(z2);
    arf_clear(tmpr);
    arf_clear(halfr);
    mag_clear(max_size);
    mag_clear(max_rad);
    mag_clear(tmpmag);

    return found;
}

