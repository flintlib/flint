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

static void
_acb_hypgeom_li(acb_t res, const acb_t z, slong prec)
{
    if (acb_is_zero(z))
    {
        acb_zero(res);
    }
    else
    {
        acb_log(res, z, prec);
        acb_hypgeom_ei(res, res, prec);
    }
}

void
_acb_hypgeom_const_li2_eval(arb_t s, slong prec)
{
    acb_t t;
    acb_init(t);
    acb_set_ui(t, 2);
    _acb_hypgeom_li(t, t, prec);
    arb_set(s, acb_realref(t));
    acb_clear(t);
}

ARB_DEF_CACHED_CONSTANT(_acb_hypgeom_const_li2, _acb_hypgeom_const_li2_eval)

static void
_acb_hypgeom_li_offset(acb_t res, const acb_t z, slong prec)
{
    if (acb_is_int(z) && arf_cmp_2exp_si(arb_midref(acb_realref(z)), 1) == 0)
    {
        acb_zero(res);
    }
    else
    {
        arb_t t;
        arb_init(t);
        _acb_hypgeom_const_li2(t, prec);
        _acb_hypgeom_li(res, z, prec);
        arb_sub(acb_realref(res), acb_realref(res), t, prec);
        arb_clear(t);
    }
}

void
acb_hypgeom_li(acb_t res, const acb_t z, int offset, slong prec)
{
    if (offset)
        _acb_hypgeom_li_offset(res, z, prec);
    else
        _acb_hypgeom_li(res, z, prec);
}

