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
    ca_sqrt_ui(res, 5, ctx);
    ca_add_ui(res, res, 1, ctx);
    ca_div_ui(res, res, 2, ctx);
}
