/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

truth_t
ca_check_is_neg_i(const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        if (ca_is_unknown(x, ctx))
            return T_UNKNOWN;

        return T_FALSE;
    }

    if (x->field == CA_FIELD_ID_QQ)
        return T_FALSE;

    if (x->field == CA_FIELD_ID_QQ_I)
    {
        const fmpz *n, *d;

        n = QNF_ELEM_NUMREF(CA_NF_ELEM(x));
        d = QNF_ELEM_DENREF(CA_NF_ELEM(x));

        if (fmpz_is_one(d) && fmpz_is_zero(n) && fmpz_equal_si(n + 1, -1))
            return T_TRUE;

        return T_FALSE;
    }

    return T_UNKNOWN;
}

