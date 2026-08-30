/*
    Copyright (C) 2022 Raoul Bourquin

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_tetranacci_constant(ca_t res, ca_ctx_t ctx)
{
    qqbar_t tc;
    qqbar_init(tc);
    qqbar_tetranacci_constant(tc);

    ca_set_qqbar(res, tc, ctx);

    qqbar_clear(tc);
}
