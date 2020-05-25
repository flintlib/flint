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
ca_extension_init_const(ca_extension_t ext, ulong func)
{
    ext->type = CA_EXT_FUNCTION;
    ext->data.function.func = func;
    ext->data.function.num_args = 0;
    ext->data.function.args = NULL;
    acb_init(&ext->data.function.enclosure);

    if (func == CA_Pi)
        acb_const_pi(&ext->data.function.enclosure, 128);
    else
        flint_abort();
}
