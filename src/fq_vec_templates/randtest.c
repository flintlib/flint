/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"


void
_TEMPLATE(T, vec_randtest) (TEMPLATE(T, struct) * f,
                            flint_rand_t state,
                            slong len, const TEMPLATE(T, ctx_t) ctx)
{
    slong i, sparseness;

    if (n_randint(state, 2))
    {
        for (i = 0; i < len; i++)
            TEMPLATE(T, randtest) (f + i, state, ctx);
    }
    else
    {
        sparseness = 1 + n_randint(state, FLINT_MAX(2, len));

        for (i = 0; i < len; i++)
        {
            if (n_randint(state, sparseness))
                TEMPLATE(T, zero) (f + i, ctx);
            else
                TEMPLATE(T, randtest) (f + i, state, ctx);
        }
    }
}

#endif
