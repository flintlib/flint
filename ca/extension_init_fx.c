/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_extension_init_fx(ca_extension_t ext, ulong func, const ca_t x)
{
    ext->type = CA_EXT_FUNCTION;
    ext->data.function.func = func;
    ext->data.function.num_args = 1;
/*
    ext->data.function.args = ca_vec_init(1);
    ca_set(ext->data.function.args, x);
*/
}
