/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fexpr.h"

void
fexpr_call_vec(fexpr_t res, const fexpr_t f, fexpr_srcptr args, slong len)
{
    if (len == 0)
    {
        fexpr_call0(res, f);
    }
    else if (len == 1)
    {
        fexpr_call1(res, f, args);
    }
    else if (len == 2)
    {
        fexpr_call2(res, f, args, args + 1);
    }
    else if (len == 3)
    {
        fexpr_call3(res, f, args, args + 1, args + 2);
    }
    else if (len == 4)
    {
        fexpr_call4(res, f, args, args + 1, args + 2, args + 3);
    }
    else
    {
        slong i, f_size, args_size, index_size, size, pos, arg_size;
        mp_ptr out;

        f_size = fexpr_size(f);

        args_size = 0;
        for (i = 0; i < len; i++)
            args_size += fexpr_size(args + i);

        /* write index:
            data[1] = nargs
            data[2] = position of f
            data[3], data[4], ..., positions of every 1/4 args for random access */
        index_size = 2 + (len + 4 - 1) / 4;

        size = 1 + index_size + f_size + args_size;

        fexpr_fit_size(res, size);
        out = res->data;

        out[0] = FEXPR_TYPE_CALLN | (size << FEXPR_TYPE_BITS);
        out[1] = len;

        pos = 1 + index_size;
        out[2] = pos;
        flint_mpn_copyi(out + pos, f->data, f_size);
        pos += f_size;

        for (i = 0; i < len; i++)
        {
            if (i % 4 == 0)
                out[3 + i / 4] = pos;

            arg_size = fexpr_size(args + i);
            flint_mpn_copyi(out + pos, args[i].data, arg_size);
            pos += arg_size;
        }
    }
}
