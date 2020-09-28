/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "utils_flint.h"

pair_t
fmpz_mpoly_select_pop_pair(pairs_t pairs, const fmpz_mpoly_vec_t G, const fmpz_mpoly_ctx_t ctx)
{
    slong len, choice, nvars;
    pair_t result;

    nvars = ctx->minfo->nvars;
    len = pairs->length;
    choice = 0;

    if (len > 1)
    {
        slong i, j;
        ulong * tmp;
        ulong * lcm;
        ulong * best_lcm;
        int best;

        lcm = flint_malloc(sizeof(ulong) * (nvars + 1));
        tmp = flint_malloc(sizeof(ulong) * (nvars + 1));
        best_lcm = flint_malloc(sizeof(ulong) * (nvars + 1));

        for (i = 0; i < len; i++)
        {
            fmpz_mpoly_get_term_exp_ui(lcm, G->p + pairs->pairs[i].a, 0, ctx);
            fmpz_mpoly_get_term_exp_ui(tmp, G->p + pairs->pairs[i].b, 0, ctx);

            lcm[nvars] = 0;

            for (j = 0; j < nvars; j++)
            {
                lcm[j] = FLINT_MAX(lcm[j], tmp[j]);
                lcm[nvars] += lcm[j];  /* total degree */
            }

            if (i == 0)
            {
                for (j = 0; j < nvars + 1; j++)
                    best_lcm[j] = lcm[j];
            }
            else
            {
                /* todo: support all orderings */
                best = 1;

                if (ctx->minfo->ord == ORD_LEX)
                {
                    for (j = 0; j < nvars && best; j++)
                        if (lcm[j] > best_lcm[j])
                            best = 0;
                }
                else
                {
                    for (j = 0; j < nvars + 1 && best; j++)
                        if (lcm[j] > best_lcm[j])
                            best = 0;
                }

                if (best)
                {
                    for (j = 0; j < nvars + 1; j++)
                        best_lcm[j] = lcm[j];
                    choice = i;
                }
            }
        }

        flint_free(tmp);
        flint_free(lcm);
        flint_free(best_lcm);
    }

    result = pairs->pairs[choice];
    pairs->pairs[choice] = pairs->pairs[pairs->length - 1];
    pairs->length--;

    return result;
}
