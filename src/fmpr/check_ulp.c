/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

int
fmpr_check_ulp(const fmpr_t result, slong r, slong prec)
{
    if (r == FMPR_RESULT_EXACT)
    {
        return 1;
    }
    else
    {

        fmpr_t err;
        fmpr_t ulp;
        int ok;

        fmpr_init(err);
        fmpr_init(ulp);

        fmpr_ulp(ulp, result, prec);
        fmpr_set_error_result(err, result, r);

        ok = fmpr_equal(err, ulp);

        fmpr_clear(err);
        fmpr_clear(ulp);

        return ok;
    }
}

