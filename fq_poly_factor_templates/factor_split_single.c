/*
    Copyright (C) 2019 Edouard Rousseau

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
TEMPLATE(T, poly_factor_split_single) (TEMPLATE(T, poly_t) linfactor,
                          const TEMPLATE(T, poly_t) input,
                          const TEMPLATE(T, ctx_t) ctx)
{
    if (input->length == 2)
    {
        TEMPLATE(T, poly_set) (linfactor, input, ctx);
    }
    else
    {
        flint_rand_t state;
        ulong deflation;
        TEMPLATE(T, poly_t) pol;

        flint_randinit(state);
        TEMPLATE(T, poly_init) (pol, ctx);
        TEMPLATE(T, poly_set) (linfactor, input, ctx);
        deflation = TEMPLATE(T, poly_deflation) (input, ctx);

        if (deflation == 1 || deflation == TEMPLATE(T, poly_degree)(input, ctx))
        {
            TEMPLATE(T, poly_set) (pol, input, ctx);

            while (TEMPLATE(T, poly_degree)(linfactor, ctx) != 1) {
                while (!TEMPLATE(T, poly_factor_equal_deg_prob)
                       (linfactor, state, pol, 1, ctx))
                {
                };
                TEMPLATE(T, poly_set) (pol, linfactor, ctx);
            }   
        }
        else
        {
            TEMPLATE(T, poly_deflate) (pol, input, deflation, ctx);

            while (TEMPLATE(T, poly_degree)(pol, ctx) != 1) {
                while (!TEMPLATE(T, poly_factor_equal_deg_prob)
                       (linfactor, state, pol, 1, ctx))
                {
                };
                TEMPLATE(T, poly_set) (pol, linfactor, ctx);
            }   
            
            TEMPLATE(T, poly_inflate) (pol, linfactor, deflation, ctx);

            while (TEMPLATE(T, poly_degree)(pol, ctx) != 1) {
                while (!TEMPLATE(T, poly_factor_equal_deg_prob)
                       (linfactor, state, pol, 1, ctx))
                {
                };
                TEMPLATE(T, poly_set) (pol, linfactor, ctx);
            }   
        }

        flint_randclear(state);
        TEMPLATE(T, poly_clear) (pol, ctx);
    }
}
#endif
