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
acb_modular_elliptic_k(acb_t k, const acb_t m, slong prec)
{
    acb_elliptic_k(k, m, prec);
}

