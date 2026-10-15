/*
    Copyright (C) 2018 Daniel Schultz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "mpoly.h"
#include "fmpz_mpoly.h"

/* evaluate B(x_1,...,x_n) at x_i = y_c[i], y_j are vars of ctxAC (or
   x_i = 0 when c[i] < 0): a renaming of the variables, done term by
   term in time linear in the number of variables */
void fmpz_mpoly_compose_fmpz_mpoly_gen(fmpz_mpoly_t A,
                             const fmpz_mpoly_t B, const slong * c,
                     const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC)
{
    slong i, v, Blen = B->length;
    slong nB = ctxB->minfo->nvars, nA = ctxAC->minfo->nvars;
    fmpz_mpoly_t T;
    int monotone;

    if (Blen == 0)
    {
        fmpz_mpoly_zero(A, ctxAC);
        return;
    }

    /* an injective order-preserving renaming between contexts with the
       same ordering keeps the terms sorted and distinct */
    monotone = (ctxB->minfo->ord == ctxAC->minfo->ord);
    for (v = 0; v < nB && monotone; v++)
        if (c[v] < 0 || (v > 0 && c[v] <= c[v - 1]))
            monotone = 0;

    fmpz_mpoly_init(T, ctxAC);
    fmpz_mpoly_fit_length(T, Blen, ctxAC);

    if (B->bits <= FLINT_BITS)
    {
        ulong * eB = flint_malloc(sizeof(ulong) * nB);
        ulong * eA = flint_malloc(sizeof(ulong) * nA);

        for (i = 0; i < Blen; i++)
        {
            int zero = 0;

            fmpz_mpoly_get_term_exp_ui(eB, B, i, ctxB);
            for (v = 0; v < nA; v++)
                eA[v] = 0;
            for (v = 0; v < nB; v++)
            {
                if (eB[v] != 0)
                {
                    if (c[v] < 0)
                        zero = 1;
                    else
                        eA[c[v]] += eB[v];
                }
            }

            if (!zero)
                fmpz_mpoly_push_term_fmpz_ui(T, B->coeffs + i, eA, ctxAC);
        }

        flint_free(eB);
        flint_free(eA);
    }
    else
    {
        fmpz * eB = _fmpz_vec_init(nB);
        fmpz * eA = _fmpz_vec_init(nA);
        fmpz ** pB = flint_malloc(sizeof(fmpz *) * nB);
        fmpz ** pA = flint_malloc(sizeof(fmpz *) * nA);

        for (v = 0; v < nB; v++)
            pB[v] = eB + v;
        for (v = 0; v < nA; v++)
            pA[v] = eA + v;

        for (i = 0; i < Blen; i++)
        {
            int zero = 0;

            fmpz_mpoly_get_term_exp_fmpz(pB, B, i, ctxB);
            _fmpz_vec_zero(eA, nA);
            for (v = 0; v < nB; v++)
            {
                if (!fmpz_is_zero(eB + v))
                {
                    if (c[v] < 0)
                        zero = 1;
                    else
                        fmpz_add(eA + c[v], eA + c[v], eB + v);
                }
            }

            if (!zero)
                fmpz_mpoly_push_term_fmpz_fmpz(T, B->coeffs + i, pA, ctxAC);
        }

        _fmpz_vec_clear(eB, nB);
        _fmpz_vec_clear(eA, nA);
        flint_free(pB);
        flint_free(pA);
    }

    if (!monotone)
    {
        fmpz_mpoly_sort_terms(T, ctxAC);
        fmpz_mpoly_combine_like_terms(T, ctxAC);
    }

    fmpz_mpoly_swap(A, T, ctxAC);
    fmpz_mpoly_clear(T, ctxAC);
}
