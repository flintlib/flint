/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_vec.h"

truth_t
_ca_vec_check_is_zero(ca_srcptr vec, slong len, ca_ctx_t ctx)
{
    slong i;
    int have_unknown;
    truth_t is_zero;

    have_unknown = 0;
    for (i = 0; i < len; i++)
    {
        is_zero = ca_check_is_zero(vec + i, ctx);

        if (is_zero == T_FALSE)
            return T_FALSE;

        if (is_zero == T_UNKNOWN)
            have_unknown = 1;
    }

    if (have_unknown)
        return T_UNKNOWN;
    else
        return T_TRUE;
}
