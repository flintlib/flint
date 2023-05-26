/*
    Copyright (C) 2020 Fredrik Johansson
    Copyright (C) 2022 Raoul Bourquin

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_phi(ca_t res, ca_ctx_t ctx)
{
    qqbar_t phi;
    qqbar_init(phi);
    qqbar_phi(phi);

    ca_set_qqbar(res, phi, ctx);

    qqbar_clear(phi);
}
