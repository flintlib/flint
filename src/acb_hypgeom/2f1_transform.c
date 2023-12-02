/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_hypgeom.h"

/*

which == 0 -- z
which == 1 -- z/(z-1)
which == 2 -- 1/z
which == 3 -- 1/(1-z)
which == 4 -- 1-z
which == 5 -- 1-1/z

*/

void
_acb_hypgeom_2f1_transform_limit(acb_t res, const acb_poly_t a, const acb_poly_t b,
    const acb_poly_t c, const acb_poly_t z, int which, slong prec)
{
    acb_poly_t ba, ca, cb, cab, ac1, bc1, ab1, ba1, w, t, u, v, s;
    acb_t tt;

    acb_poly_init(ba);
    acb_poly_init(ca); acb_poly_init(cb); acb_poly_init(cab);
    acb_poly_init(ac1); acb_poly_init(bc1);
    acb_poly_init(ab1); acb_poly_init(ba1);
    acb_poly_init(w); acb_poly_init(t);
    acb_poly_init(u); acb_poly_init(v);
    acb_poly_init(s);
    acb_init(tt);

    acb_poly_add_si(s, z, -1, prec);   /* s = 1 - z */
    acb_poly_neg(s, s);
    acb_poly_sub(ba, b, a, prec);      /* ba = b - a */
    acb_poly_sub(ca, c, a, prec);      /* ca = c - a */
    acb_poly_sub(cb, c, b, prec);      /* cb = c - b */
    acb_poly_sub(cab, ca, b, prec);    /* cab = c - a - b */
    acb_poly_add_si(ac1, ca, -1, prec); acb_poly_neg(ac1, ac1); /* ac1 = a - c + 1 */
    acb_poly_add_si(bc1, cb, -1, prec); acb_poly_neg(bc1, bc1); /* bc1 = b - c + 1 */
    acb_poly_add_si(ab1, ba, -1, prec); acb_poly_neg(ab1, ab1); /* ab1 = a - b + 1 */
    acb_poly_add_si(ba1, ba, 1, prec);                          /* ba1 = b - a + 1 */

    /* t = left term, u = right term (DLMF 15.8.1 - 15.8.5) */
    if (which == 2)
    {
        acb_poly_inv_series(w, z, 2, prec);  /* w = 1/z */
        acb_hypgeom_2f1_series_direct(t, a, ac1, ab1, w, 1, 2, prec);
        acb_hypgeom_2f1_series_direct(u, b, bc1, ba1, w, 1, 2, prec);
    }
    else if (which == 3)
    {
        acb_poly_inv_series(w, s, 2, prec);  /* w = 1/(1-z) */
        acb_hypgeom_2f1_series_direct(t, a, cb, ab1, w, 1, 2, prec);
        acb_hypgeom_2f1_series_direct(u, b, ca, ba1, w, 1, 2, prec);
    }
    else if (which == 4)
    {
        acb_poly_set(w, s);                  /* w = 1-z */
        acb_poly_add(v, ac1, b, prec);       /* v = a+b-c+1 */
        acb_hypgeom_2f1_series_direct(t, a, b, v, w, 1, 2, prec);
        acb_poly_add_si(v, cab, 1, prec);    /* v = c-a-b+1 */
        acb_hypgeom_2f1_series_direct(u, ca, cb, v, w, 1, 2, prec);
    }
    else if (which == 5)
    {
        acb_poly_inv_series(w, z, 2, prec);  /* w = 1-1/z */
        acb_poly_neg(w, w);
        acb_poly_add_si(w, w, 1, prec);
        acb_poly_add(v, ac1, b, prec);       /* v = a+b-c+1 */
        acb_hypgeom_2f1_series_direct(t, a, ac1, v, w, 1, 2, prec);
        acb_poly_add_si(v, cab, 1, prec);    /* v = c-a-b+1 */
        acb_poly_add_si(u, a, -1, prec);     /* u = 1-a */
        acb_poly_neg(u, u);
        acb_hypgeom_2f1_series_direct(u, ca, u, v, w, 1, 2, prec);
    }
    else
    {
        flint_throw(FLINT_ERROR, "invalid transformation!\n");
    }

    /* gamma factors */
    acb_poly_rgamma_series(v, a, 2, prec);
    acb_poly_mullow(u, u, v, 2, prec);
    acb_poly_rgamma_series(v, ca, 2, prec);
    acb_poly_mullow(t, t, v, 2, prec);

    acb_poly_rgamma_series(v, b, 2, prec);
    if (which == 2 || which == 3)
        acb_poly_mullow(t, t, v, 2, prec);
    else
        acb_poly_mullow(u, u, v, 2, prec);

    acb_poly_rgamma_series(v, cb, 2, prec);
    if (which == 2 || which == 3)
        acb_poly_mullow(u, u, v, 2, prec);
    else
        acb_poly_mullow(t, t, v, 2, prec);

    if (which == 2 || which == 3)
    {
        if (which == 2)
            acb_poly_neg(s, z);  /* -z, otherwise 1-z since before */

        acb_poly_neg(v, a);
        acb_poly_pow_series(v, s, v, 2, prec);
        acb_poly_mullow(t, t, v, 2, prec);

        acb_poly_neg(v, b);
        acb_poly_pow_series(v, s, v, 2, prec);
        acb_poly_mullow(u, u, v, 2, prec);
    }
    else
    {
        acb_poly_pow_series(v, s, cab, 2, prec);
        acb_poly_mullow(u, u, v, 2, prec);

        if (which == 5)
        {
            acb_poly_neg(v, a);
            acb_poly_pow_series(v, z, v, 2, prec);
            acb_poly_mullow(t, t, v, 2, prec);

            acb_poly_neg(v, ca);
            acb_poly_pow_series(v, z, v, 2, prec);
            acb_poly_mullow(u, u, v, 2, prec);
        }
    }

    acb_poly_sub(t, t, u, prec);

    if (which == 2 || which == 3)
        acb_poly_sin_pi_series(v, ba, 2, prec);
    else
        acb_poly_sin_pi_series(v, cab, 2, prec);

    acb_poly_get_coeff_acb(tt, t, 1);
    acb_poly_get_coeff_acb(res, v, 1);
    acb_div(res, tt, res, prec);
    acb_const_pi(tt, prec);
    acb_mul(res, res, tt, prec);

    acb_poly_clear(ba);
    acb_poly_clear(ca); acb_poly_clear(cb); acb_poly_clear(cab);
    acb_poly_clear(ac1); acb_poly_clear(bc1);
    acb_poly_clear(ab1); acb_poly_clear(ba1);
    acb_poly_clear(w); acb_poly_clear(t);
    acb_poly_clear(u); acb_poly_clear(v);
    acb_poly_clear(s);
    acb_clear(tt);
}

void
acb_hypgeom_2f1_transform_limit(acb_t res, const acb_t a, const acb_t b,
    const acb_t c, const acb_t z, int regularized, int which, slong prec)
{
    acb_poly_t aa, bb, cc, zz;
    acb_t t;

    if (acb_contains_zero(z) || !acb_is_finite(z))
    {
        acb_indeterminate(res);
        return;
    }

    if (arb_contains_si(acb_realref(z), 1) && arb_contains_zero(acb_imagref(z)))
    {
        acb_indeterminate(res);
        return;
    }

    if (!regularized)
    {
        acb_init(t);
        acb_gamma(t, c, prec);
        acb_hypgeom_2f1_transform_limit(res, a, b, c, z, 1, which, prec);
        acb_mul(res, res, t, prec);
        acb_clear(t);
        return;
    }

    acb_poly_init(aa);
    acb_poly_init(bb);
    acb_poly_init(cc);
    acb_poly_init(zz);
    acb_init(t);

    acb_poly_set_acb(aa, a);
    acb_poly_set_acb(bb, b);
    acb_poly_set_acb(cc, c);
    acb_poly_set_acb(zz, z);

    if (which == 2 || which == 3)
    {
        acb_sub(t, b, a, prec);
        acb_poly_set_coeff_si(aa, 1, 1);

        /* prefer b-a nonnegative (either is correct) to avoid
           expensive operations in the hypergeometric series */
        if (arb_is_nonnegative(acb_realref(t)))
            _acb_hypgeom_2f1_transform_limit(res, aa, bb, cc, zz, which, prec);
        else
            _acb_hypgeom_2f1_transform_limit(res, bb, aa, cc, zz, which, prec);
    }
    else
    {
        acb_poly_set_coeff_si(aa, 1, 1);
        _acb_hypgeom_2f1_transform_limit(res, aa, bb, cc, zz, which, prec);
    }

    acb_poly_clear(aa);
    acb_poly_clear(bb);
    acb_poly_clear(cc);
    acb_poly_clear(zz);
    acb_clear(t);
}

void
acb_hypgeom_2f1_transform_nolimit(acb_t res, const acb_t a, const acb_t b,
    const acb_t c, const acb_t z, int regularized, int which, slong prec)
{
    acb_t ba, ca, cb, cab, ac1, bc1, ab1, ba1, w, t, u, v, s;

    if (acb_contains_zero(z) || !acb_is_finite(z))
    {
        acb_indeterminate(res);
        return;
    }

    if (arb_contains_si(acb_realref(z), 1) && arb_contains_zero(acb_imagref(z)))
    {
        acb_indeterminate(res);
        return;
    }

    if (!regularized)
    {
        acb_init(t);
        acb_gamma(t, c, prec);
        acb_hypgeom_2f1_transform_nolimit(res, a, b, c, z, 1, which, prec);
        acb_mul(res, res, t, prec);
        acb_clear(t);
        return;
    }

    acb_init(ba);
    acb_init(ca); acb_init(cb); acb_init(cab);
    acb_init(ac1); acb_init(bc1);
    acb_init(ab1); acb_init(ba1);
    acb_init(w); acb_init(t);
    acb_init(u); acb_init(v);
    acb_init(s);

    acb_add_si(s, z, -1, prec);   /* s = 1 - z */
    acb_neg(s, s);

    acb_sub(ba, b, a, prec);      /* ba = b - a */
    acb_sub(ca, c, a, prec);      /* ca = c - a */
    acb_sub(cb, c, b, prec);      /* cb = c - b */
    acb_sub(cab, ca, b, prec);    /* cab = c - a - b */

    acb_add_si(ac1, ca, -1, prec); acb_neg(ac1, ac1); /* ac1 = a - c + 1 */
    acb_add_si(bc1, cb, -1, prec); acb_neg(bc1, bc1); /* bc1 = b - c + 1 */
    acb_add_si(ab1, ba, -1, prec); acb_neg(ab1, ab1); /* ab1 = a - b + 1 */
    acb_add_si(ba1, ba, 1, prec);                     /* ba1 = b - a + 1 */

    /* t = left term, u = right term (DLMF 15.8.1 - 15.8.5) */
    if (which == 2)
    {
        acb_inv(w, z, prec);  /* w = 1/z */
        acb_hypgeom_2f1_direct(t, a, ac1, ab1, w, 1, prec);
        acb_hypgeom_2f1_direct(u, b, bc1, ba1, w, 1, prec);
    }
    else if (which == 3)
    {
        acb_inv(w, s, prec);  /* w = 1/(1-z) */
        acb_hypgeom_2f1_direct(t, a, cb, ab1, w, 1, prec);
        acb_hypgeom_2f1_direct(u, b, ca, ba1, w, 1, prec);
    }
    else if (which == 4)
    {
        acb_set(w, s);                  /* w = 1-z */
        acb_add(v, ac1, b, prec);       /* v = a+b-c+1 */
        acb_hypgeom_2f1_direct(t, a, b, v, w, 1, prec);
        acb_add_si(v, cab, 1, prec);    /* v = c-a-b+1 */
        acb_hypgeom_2f1_direct(u, ca, cb, v, w, 1, prec);
    }
    else if (which == 5)
    {
        acb_inv(w, z, prec);  /* w = 1-1/z */
        acb_neg(w, w);
        acb_add_si(w, w, 1, prec);
        acb_add(v, ac1, b, prec);       /* v = a+b-c+1 */
        acb_hypgeom_2f1_direct(t, a, ac1, v, w, 1, prec);
        acb_add_si(v, cab, 1, prec);    /* v = c-a-b+1 */
        acb_add_si(u, a, -1, prec);     /* u = 1-a */
        acb_neg(u, u);
        acb_hypgeom_2f1_direct(u, ca, u, v, w, 1, prec);
    }
    else
    {
        flint_throw(FLINT_ERROR, "invalid transformation!\n");
    }

    /* gamma factors */
    acb_rgamma(v, a, prec);
    acb_mul(u, u, v, prec);
    acb_rgamma(v, ca, prec);
    acb_mul(t, t, v, prec);

    acb_rgamma(v, b, prec);
    if (which == 2 || which == 3)
        acb_mul(t, t, v, prec);
    else
        acb_mul(u, u, v, prec);

    acb_rgamma(v, cb, prec);
    if (which == 2 || which == 3)
        acb_mul(u, u, v, prec);
    else
        acb_mul(t, t, v, prec);

    if (which == 2 || which == 3)
    {
        if (which == 2)
            acb_neg(s, z);  /* -z, otherwise 1-z since before */

        acb_neg(v, a);
        acb_pow(v, s, v, prec);
        acb_mul(t, t, v, prec);

        acb_neg(v, b);
        acb_pow(v, s, v, prec);
        acb_mul(u, u, v, prec);
    }
    else
    {
        acb_pow(v, s, cab, prec);
        acb_mul(u, u, v, prec);

        if (which == 5)
        {
            acb_neg(v, a);
            acb_pow(v, z, v, prec);
            acb_mul(t, t, v, prec);

            acb_neg(v, ca);
            acb_pow(v, z, v, prec);
            acb_mul(u, u, v, prec);
        }
    }

    acb_sub(t, t, u, prec);

    if (which == 2 || which == 3)
        acb_sin_pi(v, ba, prec);
    else
        acb_sin_pi(v, cab, prec);

    acb_div(t, t, v, prec);
    acb_const_pi(v, prec);
    acb_mul(t, t, v, prec);
    acb_set(res, t);

    acb_clear(ba);
    acb_clear(ca); acb_clear(cb); acb_clear(cab);
    acb_clear(ac1); acb_clear(bc1);
    acb_clear(ab1); acb_clear(ba1);
    acb_clear(w); acb_clear(t);
    acb_clear(u); acb_clear(v);
    acb_clear(s);
}

void
acb_hypgeom_2f1_transform(acb_t res, const acb_t a, const acb_t b,
    const acb_t c, const acb_t z, int flags, int which, slong prec)
{
    int regularized;

    regularized = flags & ACB_HYPGEOM_2F1_REGULARIZED;

    if (which == 1)
    {
        acb_t t, u, v;

        acb_init(t);
        acb_init(u);
        acb_init(v);

        acb_sub_ui(t, z, 1, prec); /* t = z-1 */
        acb_div(u, z, t, prec); /* u = z/(z-1) */
        acb_neg(t, t);
        acb_neg(v, a);
        acb_pow(t, t, v, prec); /* t = (1-z)^-a */
        acb_sub(v, c, b, prec); /* v = c-b */

        /* We cannot use regularized=1 directly, since if c is a nonnegative
           integer, the transformation formula reads (lhs) * 0 = (rhs) * 0. */
        acb_hypgeom_2f1_direct(u, a, v, c, u, 1, prec);

        if (!regularized)
        {
            acb_gamma(v, c, prec);
            acb_mul(u, u, v, prec);
        }

        acb_mul(res, u, t, prec);


        acb_clear(t);
        acb_clear(u);
        acb_clear(v);
    }
    else
    {
        acb_t d;
        int limit;

        acb_init(d);

        if (which == 2 || which == 3)
        {
            if (flags & ACB_HYPGEOM_2F1_AB)
            {
                limit = 1;
            }
            else
            {
                acb_sub(d, b, a, prec);
                limit = acb_is_int(d);
            }
        }
        else
        {
            if (flags & ACB_HYPGEOM_2F1_ABC)
            {
                limit = 1;
            }
            else
            {
                acb_sub(d, c, a, prec);
                acb_sub(d, d, b, prec);
                limit = acb_is_int(d);
            }
        }

        if (limit)
            acb_hypgeom_2f1_transform_limit(res, a, b, c, z, regularized, which, prec);
        else
            acb_hypgeom_2f1_transform_nolimit(res, a, b, c, z, regularized, which, prec);

        acb_clear(d);
    }

    if (!acb_is_finite(res))
        acb_indeterminate(res);
}

