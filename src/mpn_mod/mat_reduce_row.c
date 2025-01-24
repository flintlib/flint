/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_mod.h"

/* todo: optimize for when 2n or 2n-1 rather than 2n+1 limbs suffice */
int mpn_mod_mat_reduce_row(slong * column, gr_mat_t A, slong * P, slong * L,
                                         slong m, gr_ctx_t ctx)
{
    slong n = A->c, i, j, r, res = -WORD(1);
    int status = GR_SUCCESS;
    slong Astride = A->stride;
    slong nlimbs = MPN_MOD_CTX_NLIMBS(ctx);
    slong tnlimbs = 2 * nlimbs + 1;
    nn_ptr aa = A->entries;
    nn_ptr tmp, h, u;
    slong alloc;

    TMP_INIT;
    TMP_START;

    alloc = (n * tnlimbs) + (2 * nlimbs) + (nlimbs);

    tmp = TMP_ALLOC(alloc * sizeof(ulong));
    u = tmp + (n * tnlimbs);
    h = u + (2 * nlimbs);

#define ENTRY(ii, jj) (aa + (((ii) * Astride + (jj)) * nlimbs))
#define TMP(ii) (tmp + tnlimbs * (ii))

    for (i = 0; i < n; i++)
    {
        flint_mpn_copyi(TMP(i), ENTRY(m, i), nlimbs);
        flint_mpn_zero(TMP(i) + nlimbs, tnlimbs - nlimbs);
    }

    for (i = 0; i < n; i++)
    {
        if (i != 0)
        {
            mpn_mod_set_mpn(ENTRY(m, i), TMP(i), tnlimbs, ctx);
        }

        if (mpn_mod_is_zero(ENTRY(m, i), ctx) == T_FALSE)
        {
            r = P[i];

            if (r != -WORD(1))
            {
                mpn_mod_neg(h, ENTRY(m, i), ctx);
                mpn_mod_zero(ENTRY(m, i), ctx);

                if (nlimbs == 2)
                {
                    for (j = i + 1; j < L[r]; j++)
                    {
                        ulong t[4];
                        nn_srcptr ap = ENTRY(r, j);
                        nn_ptr s = TMP(j);

                        FLINT_MPN_MUL_2X2(t[3], t[2], t[1], t[0], ap[1], ap[0], h[1], h[0]);
                        add_sssssaaaaaaaaaa(s[4], s[3], s[2], s[1], s[0],
                                            s[4], s[3], s[2], s[1], s[0],
                                            0, t[3], t[2], t[1], t[0]);
                    }
                }
                else
                {
                    for (j = i + 1; j < L[r]; j++)
                    {
                        flint_mpn_mul_n(u, ENTRY(r, j), h, nlimbs);
                        TMP(j)[tnlimbs - 1] += mpn_add_n(TMP(j), TMP(j), u, 2 * nlimbs);
                    }
                }
            }
            else
            {
                status = mpn_mod_inv(h, ENTRY(m, i), ctx);
                if (status != GR_SUCCESS)
                    goto cleanup;

                mpn_mod_one(ENTRY(m, i), ctx);

                for (j = i + 1; j < L[m]; j++)
                {
                    mpn_mod_set_mpn(ENTRY(m, j), TMP(j), tnlimbs, ctx);
                    mpn_mod_mul(ENTRY(m, j), ENTRY(m, j), h, ctx);
                }

                P[i] = m;
                res = i;
                break;

            }
        }
    }

cleanup:

    TMP_END;

    *column = res;
    return status;
}
