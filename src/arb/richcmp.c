/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

int arb_eq(const arb_t x, const arb_t y)
{
    if (arf_is_nan(arb_midref(x)) || arf_is_nan(arb_midref(y)))
        return 0;

    if (mag_is_inf(arb_radref(x)) || mag_is_inf(arb_radref(y)))
        return 0;

    if (arf_is_inf(arb_midref(x)) || arf_is_inf(arb_midref(y)) ||
        (mag_is_zero(arb_radref(x)) && mag_is_zero(arb_radref(y))))
        return arf_equal(arb_midref(x), arb_midref(y));

    return 0;
}

int arb_ne(const arb_t x, const arb_t y)
{
    if (arf_is_nan(arb_midref(x)) || arf_is_nan(arb_midref(y)))
        return 0;

    if (mag_is_inf(arb_radref(x)) || mag_is_inf(arb_radref(y)))
        return 0;

    if (arf_is_inf(arb_midref(x)) || arf_is_inf(arb_midref(y)) ||
        (mag_is_zero(arb_radref(x)) && mag_is_zero(arb_radref(y))))
        return !arf_equal(arb_midref(x), arb_midref(y));

    return !arb_overlaps(x, y);
}

int arb_lt(const arb_t x, const arb_t y)
{
    if (arf_is_nan(arb_midref(x)) || arf_is_nan(arb_midref(y)))
        return 0;

    if (mag_is_inf(arb_radref(x)) || mag_is_inf(arb_radref(y)))
        return 0;

    if (arf_is_inf(arb_midref(x)) || arf_is_inf(arb_midref(y)) ||
        (mag_is_zero(arb_radref(x)) && mag_is_zero(arb_radref(y))))
        return arf_cmp(arb_midref(x), arb_midref(y)) < 0;

    if (arb_overlaps(x, y))
        return 0;

    return arf_cmp(arb_midref(x), arb_midref(y)) < 0;
}

int arb_le(const arb_t x, const arb_t y)
{
    arf_struct u[4];
    arf_t t;
    int res;

    if (arf_is_nan(arb_midref(x)) || arf_is_nan(arb_midref(y)))
        return 0;

    if (mag_is_inf(arb_radref(x)) || mag_is_inf(arb_radref(y)))
        return (arf_is_neg_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)))
            || (arf_is_pos_inf(arb_midref(y)) && mag_is_finite(arb_radref(y)));

    if (arf_is_inf(arb_midref(x)) || arf_is_inf(arb_midref(y)) ||
        (mag_is_zero(arb_radref(x)) && mag_is_zero(arb_radref(y))))
        return arf_cmp(arb_midref(x), arb_midref(y)) <= 0;

    if (!arb_overlaps(x, y))
        return arf_cmp(arb_midref(x), arb_midref(y)) < 0;

    /* need to check for the case where a single point overlaps */
    /* todo: optimize this */
    /* xm + xr <= ym - yr  <=>  xm + xr - ym + yr <= 0 */
    arf_init_set_shallow(u + 0, arb_midref(x));
    arf_init_neg_shallow(u + 1, arb_midref(y));
    arf_init_set_mag_shallow(u + 2, arb_radref(x));
    arf_init_set_mag_shallow(u + 3, arb_radref(y));

    arf_init(t);
    arf_sum(t, u, 4, MAG_BITS, ARF_RND_DOWN);
    res = (arf_sgn(t) <= 0);
    arf_clear(t);

    return res;
}

int arb_gt(const arb_t x, const arb_t y)
{
    if (arf_is_nan(arb_midref(x)) || arf_is_nan(arb_midref(y)))
        return 0;

    if (mag_is_inf(arb_radref(x)) || mag_is_inf(arb_radref(y)))
        return 0;

    if (arf_is_inf(arb_midref(x)) || arf_is_inf(arb_midref(y)) ||
        (mag_is_zero(arb_radref(x)) && mag_is_zero(arb_radref(y))))
        return arf_cmp(arb_midref(x), arb_midref(y)) > 0;

    if (arb_overlaps(x, y))
        return 0;

    return arf_cmp(arb_midref(x), arb_midref(y)) > 0;
}

int arb_ge(const arb_t x, const arb_t y)
{
    arf_struct u[4];
    arf_t t;
    int res;

    if (arf_is_nan(arb_midref(x)) || arf_is_nan(arb_midref(y)))
        return 0;

    if (mag_is_inf(arb_radref(x)) || mag_is_inf(arb_radref(y)))
        return (arf_is_pos_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)))
            || (arf_is_neg_inf(arb_midref(y)) && mag_is_finite(arb_radref(y)));

    if (arf_is_inf(arb_midref(x)) || arf_is_inf(arb_midref(y)) ||
        (mag_is_zero(arb_radref(x)) && mag_is_zero(arb_radref(y))))
        return arf_cmp(arb_midref(x), arb_midref(y)) >= 0;

    if (!arb_overlaps(x, y))
        return arf_cmp(arb_midref(x), arb_midref(y)) > 0;

    /* need to check for the case where a single point overlaps */
    /* todo: optimize this */
    /* xm - xr >= ym + yr  <=>  xm - xr - ym - yr >= 0 */
    arf_init_set_shallow(u + 0, arb_midref(x));
    arf_init_neg_shallow(u + 1, arb_midref(y));
    arf_init_neg_mag_shallow(u + 2, arb_radref(x));
    arf_init_neg_mag_shallow(u + 3, arb_radref(y));

    arf_init(t);
    arf_sum(t, u, 4, MAG_BITS, ARF_RND_DOWN);
    res = (arf_sgn(t) >= 0);
    arf_clear(t);

    return res;
}

int
arb_contains_zero(const arb_t x)
{
    return arf_cmpabs_mag(arb_midref(x), arb_radref(x)) <= 0;
}

int
arb_is_nonzero(const arb_t x)
{
    return !arb_contains_zero(x);
}

int
arb_is_positive(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) > 0) &&
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) < 0) &&
         !arf_is_nan(arb_midref(x));
}

int
arb_is_nonnegative(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) >= 0) &&
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) <= 0) &&
         !arf_is_nan(arb_midref(x));
}

int
arb_is_negative(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) < 0) &&
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) < 0) &&
         !arf_is_nan(arb_midref(x));
}

int
arb_is_nonpositive(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) <= 0) &&
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) <= 0) &&
         !arf_is_nan(arb_midref(x));
}

int
arb_contains_negative(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) < 0) ||
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) > 0)
        || arf_is_nan(arb_midref(x));
}

int
arb_contains_nonpositive(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) <= 0) ||
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) >= 0)
        || arf_is_nan(arb_midref(x));
}

int
arb_contains_positive(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) > 0) ||
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) > 0)
        || arf_is_nan(arb_midref(x));
}

int
arb_contains_nonnegative(const arb_t x)
{
    return (arf_sgn(arb_midref(x)) >= 0) ||
        (arf_mag_cmpabs(arb_radref(x), arb_midref(x)) >= 0)
        || arf_is_nan(arb_midref(x));
}

