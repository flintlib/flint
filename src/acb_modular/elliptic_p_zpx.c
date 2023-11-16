/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_modular.h"
#include "acb_elliptic.h"

void
acb_modular_elliptic_p_zpx(acb_ptr r, const acb_t z, const acb_t tau, slong len, slong prec)
{
    acb_elliptic_p_jet(r, z, tau, len, prec);
}

