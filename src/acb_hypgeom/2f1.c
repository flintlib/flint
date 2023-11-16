/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"
#include "arb_hypgeom.h"

static void
_acb_hypgeom_2f1r_reduced(acb_t res,
    const acb_t b, const acb_t c, const acb_t z, slong prec)
{
    acb_t t, u;
    acb_init(t);
    acb_init(u);
    acb_sub_ui(t, z, 1, prec);
    acb_neg(t, t);
    acb_neg(u, b);
    acb_pow(t, t, u, prec);
    acb_rgamma(u, c, prec);
    acb_mul(t, t, u, prec);
    acb_set(res, t);
    acb_clear(t);
    acb_clear(u);
    return;
}

void
acb_hypgeom_2f1_nointegration(acb_t res, const acb_t a, const acb_t b,
        const acb_t c, const acb_t z, int flags, slong prec)
{
    int algorithm, regularized;

    regularized = flags & ACB_HYPGEOM_2F1_REGULARIZED;

    if (!acb_is_finite(a) || !acb_is_finite(b) || !acb_is_finite(c) || !acb_is_finite(z))
    {
        acb_indeterminate(res);
        return;
    }

    if (acb_is_zero(z))
    {
        if (regularized)
            acb_rgamma(res, c, prec);
        else
            acb_one(res);
        return;
    }

    if (regularized && acb_is_int(c) && arb_is_nonpositive(acb_realref(c)))
    {
        if ((acb_is_int(a) && arb_is_nonpositive(acb_realref(a)) &&
            arf_cmp(arb_midref(acb_realref(a)), arb_midref(acb_realref(c))) >= 0) ||
            (acb_is_int(b) && arb_is_nonpositive(acb_realref(b)) &&
            arf_cmp(arb_midref(acb_realref(b)), arb_midref(acb_realref(c))) >= 0))
        {
            acb_zero(res);
            return;
        }
    }

    if (regularized && acb_eq(a, c))
    {
        _acb_hypgeom_2f1r_reduced(res, b, c, z, prec);
        return;
    }

    if (regularized && acb_eq(b, c))
    {
        _acb_hypgeom_2f1r_reduced(res, a, c, z, prec);
        return;
    }

    /* polynomial */
    if (acb_is_int(a) && arf_sgn(arb_midref(acb_realref(a))) <= 0 &&
         arf_cmpabs_ui(arb_midref(acb_realref(a)), prec) < 0)
    {
        acb_hypgeom_2f1_direct(res, a, b, c, z, regularized, prec);
        return;
    }

    /* polynomial */
    if (acb_is_int(b) && arf_sgn(arb_midref(acb_realref(b))) <= 0 &&
         arf_cmpabs_ui(arb_midref(acb_realref(b)), prec) < 0)
    {
        acb_hypgeom_2f1_direct(res, a, b, c, z, regularized, prec);
        return;
    }

    /* Try to reduce to a polynomial case using the Pfaff transformation */
    /* TODO: look at flags for integer c-b, c-a here, even when c is nonexact */
    if (acb_is_exact(c))
    {
        acb_t t;
        acb_init(t);

        acb_sub(t, c, b, prec);

        if (acb_is_int(t) && arb_is_nonpositive(acb_realref(t)))
        {
            acb_hypgeom_2f1_transform(res, a, b, c, z, flags, 1, prec);
            acb_clear(t);
            return;
        }

        acb_sub(t, c, a, prec);

        if (acb_is_int(t) && arb_is_nonpositive(acb_realref(t)))
        {
            int f1, f2;

            /* When swapping a, b, also swap the flags. */
            f1 = flags & ACB_HYPGEOM_2F1_AC;
            f2 = flags & ACB_HYPGEOM_2F1_BC;

            flags &= ~ACB_HYPGEOM_2F1_AC;
            flags &= ~ACB_HYPGEOM_2F1_BC;

            if (f1) flags |= ACB_HYPGEOM_2F1_BC;
            if (f2) flags |= ACB_HYPGEOM_2F1_AC;

            acb_hypgeom_2f1_transform(res, b, a, c, z, flags, 1, prec);
            acb_clear(t);
            return;
        }

        acb_clear(t);
    }

    /* special value at z = 1 */
    if (acb_is_one(z))
    {
        acb_t t, u, v;

        acb_init(t);
        acb_init(u);
        acb_init(v);

        acb_sub(t, c, a, prec);
        acb_sub(u, c, b, prec);
        acb_sub(v, t, b, prec);

        if (arb_is_positive(acb_realref(v)))
        {
            acb_rgamma(t, t, prec);
            acb_rgamma(u, u, prec);
            acb_mul(t, t, u, prec);
            acb_gamma(v, v, prec);
            acb_mul(t, t, v, prec);

            if (!regularized)
            {
                acb_gamma(v, c, prec);
                acb_mul(t, t, v, prec);
            }

            acb_set(res, t);
        }
        else
        {
            acb_indeterminate(res);
        }

        acb_clear(t);
        acb_clear(u);
        acb_clear(v);

        return;
    }

    algorithm = acb_hypgeom_2f1_choose(z);

    if (algorithm == 0)
    {
        acb_hypgeom_2f1_direct(res, a, b, c, z, regularized, prec);
    }
    else if (algorithm >= 1 && algorithm <= 5)
    {
        acb_hypgeom_2f1_transform(res, a, b, c, z, flags, algorithm, prec);
    }
    else
    {
        acb_hypgeom_2f1_corner(res, a, b, c, z, regularized, prec);
    }
}

void
acb_hypgeom_2f1(acb_t res, const acb_t a, const acb_t b,
        const acb_t c, const acb_t z, int flags, slong prec)
{
    acb_t res2;
    slong acc, max, t;

    acb_init(res2);

    acb_hypgeom_2f1_nointegration(res2, a, b, c, z, flags, prec);

    acc = acb_rel_accuracy_bits(res2);

    if (acc < 0.5 * prec)
    {
        max = prec;
        t = acb_rel_accuracy_bits(z);
        max = FLINT_MIN(max, t);
        t = acb_rel_accuracy_bits(a);
        max = FLINT_MIN(max, t);
        t = acb_rel_accuracy_bits(b);
        max = FLINT_MIN(max, t);
        t = acb_rel_accuracy_bits(c);
        max = FLINT_MIN(max, t);

        if (max > 2 && acc < 0.5 * max)
        {
            if (acb_is_real(a) && acb_is_real(b) && acb_is_real(c) && acb_is_real(z) &&
                arf_cmpabs_2exp_si(arb_midref(acb_realref(a)), 60) < 0 &&
                arf_cmpabs_2exp_si(arb_midref(acb_realref(b)), 60) < 0 &&
                arf_cmpabs_2exp_si(arb_midref(acb_realref(c)), 60) < 0 &&
                arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 60) < 0)
            {

                arb_hypgeom_2f1_integration(acb_realref(res),
                        acb_realref(a), acb_realref(b), acb_realref(c), acb_realref(z), flags, prec);
                arb_zero(acb_imagref(res));

                if (acb_rel_accuracy_bits(res) > acb_rel_accuracy_bits(res2) ||
                    (acb_is_finite(res) && !acb_is_finite(res2)))
                {
                    acb_swap(res, res2);
                }
            }
        }
    }

    acb_swap(res, res2);
    acb_clear(res2);
}
