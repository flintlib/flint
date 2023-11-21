/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

#define a acb_realref(x)
#define b acb_imagref(x)
#define c acb_realref(y)
#define d acb_imagref(y)
#define e acb_realref(z)
#define f acb_imagref(z)
#define ar arb_radref(a)
#define br arb_radref(b)
#define cr arb_radref(c)
#define dr arb_radref(d)

static void
_acb_sqr_fast(acb_t z, const acb_t x, slong prec)
{
    int inexact;

    mag_t am, bm, er, fr;

    mag_fast_init_set_arf(am, arb_midref(a));
    mag_fast_init_set_arf(bm, arb_midref(b));

    mag_init(er);
    mag_init(fr);

    mag_fast_addmul(er, am, ar);
    mag_fast_addmul(er, bm, br);
    mag_fast_mul_2exp_si(er, er, 1);
    mag_fast_addmul(er, ar, ar);
    mag_fast_addmul(er, br, br);

    mag_fast_addmul(fr, bm, ar);
    mag_fast_addmul(fr, am, br);
    mag_fast_addmul(fr, ar, br);
    mag_fast_mul_2exp_si(fr, fr, 1);

    inexact = arf_complex_sqr(arb_midref(e), arb_midref(f),
                    arb_midref(a), arb_midref(b), prec, ARB_RND);

    if (inexact & 1)
        arf_mag_add_ulp(arb_radref(e), er, arb_midref(e), prec);
    else
        mag_set(arb_radref(e), er);

    if (inexact & 2)
        arf_mag_add_ulp(arb_radref(f), fr, arb_midref(f), prec);
    else
        mag_set(arb_radref(f), fr);
}

static void
_acb_sqr_slow(acb_t z, const acb_t x, slong prec)
{
    int inexact;

    mag_t am, bm, er, fr;

    mag_init_set_arf(am, arb_midref(a));
    mag_init_set_arf(bm, arb_midref(b));

    mag_init(er);
    mag_init(fr);

    mag_addmul(er, am, ar);
    mag_addmul(er, bm, br);
    mag_mul_2exp_si(er, er, 1);
    mag_addmul(er, ar, ar);
    mag_addmul(er, br, br);

    mag_addmul(fr, bm, ar);
    mag_addmul(fr, am, br);
    mag_addmul(fr, ar, br);
    mag_mul_2exp_si(fr, fr, 1);

    inexact = arf_complex_sqr(arb_midref(e), arb_midref(f),
                    arb_midref(a), arb_midref(b), prec, ARB_RND);

    if (inexact & 1)
        arf_mag_add_ulp(arb_radref(e), er, arb_midref(e), prec);
    else
        mag_swap(arb_radref(e), er);

    if (inexact & 2)
        arf_mag_add_ulp(arb_radref(f), fr, arb_midref(f), prec);
    else
        mag_swap(arb_radref(f), fr);

    mag_clear(am);
    mag_clear(bm);

    mag_clear(er);
    mag_clear(fr);
}

static void
_acb_mul_fast(acb_t z, const acb_t x, const acb_t y, slong prec)
{
    int inexact;

    mag_t am, bm, cm, dm, er, fr;

    mag_fast_init_set_arf(am, arb_midref(a));
    mag_fast_init_set_arf(bm, arb_midref(b));
    mag_fast_init_set_arf(cm, arb_midref(c));
    mag_fast_init_set_arf(dm, arb_midref(d));

    mag_init(er);
    mag_init(fr);

    mag_fast_addmul(er, am, cr);
    mag_fast_addmul(er, bm, dr);
    mag_fast_addmul(er, cm, ar);
    mag_fast_addmul(er, dm, br);
    mag_fast_addmul(er, ar, cr);
    mag_fast_addmul(er, br, dr);

    mag_fast_addmul(fr, am, dr);
    mag_fast_addmul(fr, bm, cr);
    mag_fast_addmul(fr, cm, br);
    mag_fast_addmul(fr, dm, ar);
    mag_fast_addmul(fr, br, cr);
    mag_fast_addmul(fr, ar, dr);

    inexact = arf_complex_mul(arb_midref(e), arb_midref(f),
                    arb_midref(a), arb_midref(b),
                    arb_midref(c), arb_midref(d), prec, ARB_RND);

    if (inexact & 1)
        arf_mag_add_ulp(arb_radref(e), er, arb_midref(e), prec);
    else
        mag_set(arb_radref(e), er);

    if (inexact & 2)
        arf_mag_add_ulp(arb_radref(f), fr, arb_midref(f), prec);
    else
        mag_set(arb_radref(f), fr);
}

static void
_acb_mul_slow(acb_t z, const acb_t x, const acb_t y, slong prec)
{
    int inexact;

    mag_t am, bm, cm, dm, er, fr;

    mag_init_set_arf(am, arb_midref(a));
    mag_init_set_arf(bm, arb_midref(b));
    mag_init_set_arf(cm, arb_midref(c));
    mag_init_set_arf(dm, arb_midref(d));

    mag_init(er);
    mag_init(fr);

    mag_addmul(er, am, cr);
    mag_addmul(er, bm, dr);
    mag_addmul(er, cm, ar);
    mag_addmul(er, dm, br);
    mag_addmul(er, ar, cr);
    mag_addmul(er, br, dr);

    mag_addmul(fr, am, dr);
    mag_addmul(fr, bm, cr);
    mag_addmul(fr, cm, br);
    mag_addmul(fr, dm, ar);
    mag_addmul(fr, br, cr);
    mag_addmul(fr, ar, dr);

    inexact = arf_complex_mul(arb_midref(e), arb_midref(f),
                    arb_midref(a), arb_midref(b),
                    arb_midref(c), arb_midref(d), prec, ARB_RND);

    if (inexact & 1)
        arf_mag_add_ulp(arb_radref(e), er, arb_midref(e), prec);
    else
        mag_swap(arb_radref(e), er);

    if (inexact & 2)
        arf_mag_add_ulp(arb_radref(f), fr, arb_midref(f), prec);
    else
        mag_swap(arb_radref(f), fr);

    mag_clear(am);
    mag_clear(bm);
    mag_clear(cm);
    mag_clear(dm);

    mag_clear(er);
    mag_clear(fr);
}

void
acb_mul(acb_t z, const acb_t x, const acb_t y, slong prec)
{
    if (arb_is_zero(b))
    {
        arb_mul(f, d, a, prec);
        arb_mul(e, c, a, prec);
    }
    else if (arb_is_zero(d))
    {
        arb_mul(f, b, c, prec);
        arb_mul(e, a, c, prec);
    }
    else if (arb_is_zero(a))
    {
        arb_mul(e, c, b, prec);
        arb_mul(f, d, b, prec);
        acb_mul_onei(z, z);
    }
    else if (arb_is_zero(c))
    {
        arb_mul(e, a, d, prec);
        arb_mul(f, b, d, prec);
        acb_mul_onei(z, z);
    }
    /* squaring = a^2-b^2, 2ab */
    else if (x == y)
    {
        if (ARB_IS_LAGOM(a) && ARB_IS_LAGOM(b))
            _acb_sqr_fast(z, x, prec);
        else
            _acb_sqr_slow(z, x, prec);
    }
    else
    {
        if (ARB_IS_LAGOM(a) && ARB_IS_LAGOM(b) &&
            ARB_IS_LAGOM(c) && ARB_IS_LAGOM(d))
            _acb_mul_fast(z, x, y, prec);
        else
            _acb_mul_slow(z, x, y, prec);
    }
}

#undef a
#undef b
#undef c
#undef d
#undef e
#undef f
#undef ar
#undef br
#undef cr
#undef dr

