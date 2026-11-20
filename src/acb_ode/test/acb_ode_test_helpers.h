/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Generated using Claude Code (Fable 5, Opus 5.5). */

#ifndef ACB_ODE_TEST_HELPERS_H
#define ACB_ODE_TEST_HELPERS_H

#include "acb.h"
#include "gr.h"

static void
check_acb(const acb_t got, const acb_t ref, slong min_acc, const char * name)
{
    slong acc;

    FLINT_ASSERT(acb_is_finite(ref));

    if (!acb_overlaps(got, ref))
    {
        flint_printf("FAIL (incorrect result): %s\n", name);
        flint_printf("got = %{acb}\n", got);
        flint_printf("ref = %{acb}\n", ref);
        flint_abort();
    }

    if (!acb_contains_zero(ref))
    {
        acc = acb_rel_accuracy_bits(got);
    }
    else
    {
        mag_t m;

        FLINT_ASSERT(mag_cmp_2exp_si(arb_radref(acb_realref(ref)), -30) <= 0);
        FLINT_ASSERT(mag_cmp_2exp_si(arb_radref(acb_imagref(ref)), -30) <= 0);

        mag_init(m);
        acb_get_mag(m, got);
        if (mag_is_zero(m))
            acc = ARF_PREC_EXACT;
        else if (!mag_is_finite(m))
            acc = -ARF_PREC_EXACT;
        else
            acc = -fmpz_get_si(MAG_EXPREF(m));
        mag_clear(m);
    }

    if (acc < min_acc)
    {
        flint_printf("FAIL (inaccurate result): %s\n", name);
        flint_printf("got = %{acb}\n", got);
        flint_printf("ref = %{acb}\n", ref);
        flint_printf("accuracy = %wd bits, required %wd\n", acc, min_acc);
        flint_abort();
    }
}

static void
check_status(int status, const char * name)
{
    if (status != GR_SUCCESS)
    {
        flint_printf("FAIL (unexpected failure): %s\n", name);
        flint_abort();
    }
}

static void
acb_set_str(acb_t res, const char * s, slong prec)
{
    gr_ctx_t CC;

    gr_ctx_init_complex_acb(CC, prec);
    GR_MUST_SUCCEED(gr_set_str(res, s, CC));
    gr_ctx_clear(CC);
}

#endif
