/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

/* todo: better impl */
/* may want to return symbolic csgn (can still bound numerically) */
void
ca_csgn(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    ca_t re, im, zero;
    truth_t is_zero;

    if (CA_IS_SPECIAL(x))
    {
        if (ca_check_is_signed_inf(x, ctx) == T_TRUE)
        {
            ca_sgn(res, x, ctx);
            ca_csgn(res, res, ctx);
        }
        else if (ca_check_is_uinf(x, ctx) == T_TRUE || ca_check_is_undefined(x, ctx) == T_TRUE)
        {
            ca_undefined(res, ctx);
        }
        else
        {
            ca_unknown(res, ctx);
        }
        return;
    }

    is_zero = ca_check_is_zero(x, ctx);

    if (is_zero == T_TRUE)
    {
        ca_zero(res, ctx);
        return;
    }

    if (is_zero == T_UNKNOWN)
    {
        ca_unknown(res, ctx);
        return;
    }

    ca_init(re, ctx);
    ca_init(zero, ctx);

    ca_re(re, x, ctx);

    if (ca_check_gt(re, zero, ctx) == T_TRUE)
    {
        ca_one(res, ctx);
    }
    else if (ca_check_lt(re, zero, ctx) == T_TRUE)
    {
        ca_neg_one(res, ctx);
    }
    else if (ca_check_is_zero(re, ctx) == T_TRUE)
    {
        ca_init(im, ctx);
        ca_im(im, x, ctx);

        if (ca_check_gt(im, zero, ctx) == T_TRUE)
        {
            ca_one(res, ctx);
        }
        else if (ca_check_lt(im, zero, ctx) == T_TRUE)
        {
            ca_neg_one(res, ctx);
        }
        else
        {
            ca_unknown(res, ctx);
        }
        ca_clear(im, ctx);
    }
    else
    {
        ca_unknown(res, ctx);
    }

    ca_clear(re, ctx);
}
