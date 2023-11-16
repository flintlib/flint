/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

int ca_ext_can_evaluate_qqbar(const ca_ext_t x, ca_ctx_t ctx);

/* todo: move, rename */
static truth_t
ca_ext_is_algebraic(const ca_ext_t x, ca_ctx_t ctx)
{
    if (CA_EXT_IS_QQBAR(x))
        return T_TRUE;

    if (ca_ext_can_evaluate_qqbar(x, ctx))
        return T_TRUE;

    return T_UNKNOWN;
}

truth_t
ca_check_is_algebraic(const ca_t x, ca_ctx_t ctx)
{
    ca_field_srcptr field;

    if (CA_IS_SPECIAL(x))
    {
        if (ca_is_unknown(x, ctx))
            return T_UNKNOWN;

        return T_FALSE;
    }

    field = CA_FIELD(x, ctx);

    if (CA_IS_QQ(x, ctx) || CA_FIELD_IS_NF(field))
    {
        return T_TRUE;
    }
    else
    {
        slong len, i;

        len = CA_FIELD_LENGTH(field);

        /* todo: handle simple transcendental numbers, e.g. Q(i,pi) */
        /* need to verify that some generator is used in the poly */
        /* for Q(a,b,pi) we don't know, because a, b could cancel out pi */
        for (i = 0; i < len; i++)
        {
            if (ca_ext_is_algebraic(CA_FIELD_EXT_ELEM(field, i), ctx) != T_TRUE)
                return T_UNKNOWN;
        }

        return T_TRUE;
    }
}

