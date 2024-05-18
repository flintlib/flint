/*
    Copyright (C) 2014, 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_modular.h"
#include "acb_elliptic.h"

void
acb_modular_elliptic_e(acb_t res, const acb_t m, slong prec)
{
    acb_elliptic_e(res, m, prec);
}

void
acb_modular_elliptic_k(acb_t k, const acb_t m, slong prec)
{
    acb_elliptic_k(k, m, prec);
}

void
acb_modular_elliptic_k_cpx(acb_ptr w, const acb_t m, slong len, slong prec)
{
    acb_elliptic_k_jet(w, m, len, prec);
}

void
acb_modular_elliptic_p(acb_t r, const acb_t z, const acb_t tau, slong prec)
{
    acb_elliptic_p(r, z, tau, prec);
}

void
acb_modular_elliptic_p_zpx(acb_ptr r, const acb_t z, const acb_t tau, slong len, slong prec)
{
    acb_elliptic_p_jet(r, z, tau, len, prec);
}
