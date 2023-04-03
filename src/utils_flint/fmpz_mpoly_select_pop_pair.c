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
        slong i, j, a, b;
        ulong * exp;
        ulong * tmp;
        ulong * lcm;
        ulong * best_lcm;
        ulong l, total;
        int best;

        exp = flint_malloc(sizeof(ulong) * G->length * nvars);
        lcm = flint_malloc(sizeof(ulong) * (nvars + 1));
        tmp = flint_malloc(sizeof(ulong) * (nvars + 1));
        best_lcm = flint_malloc(sizeof(ulong) * (nvars + 1));

        for (i = 0; i <= nvars; i++)
            best_lcm[i] = UWORD_MAX;

        for (i = 0; i < G->length; i++)
            fmpz_mpoly_get_term_exp_ui(exp + i * nvars, G->p + i, 0, ctx);

        for (i = 0; i < len; i++)
        {
            a = pairs->pairs[i].a;
            b = pairs->pairs[i].b;
            total = 0;
            best = 1;

            if (ctx->minfo->ord == ORD_LEX)
            {
                for (j = 0; j < nvars; j++)
                {
                    l = FLINT_MAX(exp[a * nvars + j], exp[b * nvars + j]);

                    if (l > best_lcm[j])
                    {
                        best = 0;
                        break;
                    }

                    lcm[j] = l;
                    total += l;  /* total degree */
                }
            }
            else  /* todo: correct order */
            {
                for (j = 0; j < nvars; j++)
                {
                    l = FLINT_MAX(exp[a * nvars + j], exp[b * nvars + j]);
                    total += l;  /* total degree */

                    if (total >= best_lcm[j])
                    {
                        best = 0;
                        break;
                    }

                    lcm[j] = l;
                }
            }

            if (best)
            {
                for (j = 0; j < nvars; j++)
                    best_lcm[j] = lcm[j];

                best_lcm[nvars] = total;
                choice = i;
            }
        }

        flint_free(exp);
        flint_free(tmp);
        flint_free(lcm);
        flint_free(best_lcm);
    }

    result = pairs->pairs[choice];
    pairs->pairs[choice] = pairs->pairs[pairs->length - 1];
    pairs->length--;

    return result;
}
