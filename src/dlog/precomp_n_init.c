/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "dlog.h"

/* group of order n modulo mod, mod a prime and no information on n */
void
dlog_precomp_n_init(dlog_precomp_t pre, ulong a, ulong mod, ulong n, ulong num)
{
    if (n % 2 && n_is_probabprime(n))
        dlog_precomp_p_init(pre, a, mod, n, num);
    else
    {
        if (n < DLOG_TABLE_N_LIM)
        {
           dlog_precomp_small_init(pre, a, mod, n, num);
        }
        else
        {
            if (n < DLOG_BSGS_LIM)
            {
                ulong m;
                m = dlog_bsgs_size(n, num);
                pre->type = DLOG_BSGS;
                pre->cost = dlog_bsgs_init(pre->t.bsgs, a, mod, n, m);
            } else {
                pre->type = DLOG_CRT;
                pre->cost = dlog_crt_init(pre->t.crt, a, mod, n, num);
            }
        }
    }
}
