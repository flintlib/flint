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
ca_check_is_neg_one(const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        if (ca_is_unknown(x, ctx))
            return T_UNKNOWN;

        return T_FALSE;
    }
    else if (CA_IS_QQ(x, ctx))
    {
        if (fmpz_is_one(fmpq_denref(CA_FMPQ(x))) && fmpz_equal_si(fmpq_numref(CA_FMPQ(x)), -1))
            return T_TRUE;
        else
            return T_FALSE;
    }
    else if (CA_IS_QQ_I(x, ctx))
    {
        const fmpz *n, *d;

        n = QNF_ELEM_NUMREF(CA_NF_ELEM(x));
        d = QNF_ELEM_NUMREF(CA_NF_ELEM(x));

        if (fmpz_is_one(d) && fmpz_equal_si(n, -1) && fmpz_is_zero(n + 1))
            return T_TRUE;

        return T_FALSE;
    }
    else
    {
        truth_t res;
        ca_t t;
        ca_init(t, ctx);
        ca_set_si(t, -1, ctx);
        res = ca_check_equal(x, t, ctx);
        ca_clear(t, ctx);
        return res;
    }
}
