/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"
#include "ca_ext.h"
#include "ca_field.h"
#include "fmpz_mpoly.h"

void _ca_field_ideal_insert_clear_mpoly(ca_field_t K, fmpz_mpoly_t poly, fmpz_mpoly_ctx_t mctx, ca_ctx_t ctx);
void _nf_elem_get_fmpz_poly_den_shallow(fmpz_poly_t pol, fmpz_t den, const nf_elem_t a, const nf_t nf);

void
fmpz_mpoly_set_coeff_si_x(fmpz_mpoly_t poly,
        slong c,
        slong x_var, slong x_exp,
        const fmpz_mpoly_ctx_t ctx);

void
ca_field_build_ideal_gamma(ca_field_t K, ca_ctx_t ctx)
{
    calcium_func_code Fi, Fj;
    slong i, j, len, num_gamma;

    len = CA_FIELD_LENGTH(K);

    if (len < 2)
        return;

    num_gamma = 0;
    for (i = 0; i < len; i++)
    {
        Fi = CA_EXT_HEAD(CA_FIELD_EXT_ELEM(K, i));

        if (Fi == CA_Gamma)
        {
            num_gamma++;
        }
    }

    if (num_gamma < 2)
        return;

    for (i = 0; i < len; i++)
    {
        Fi = CA_EXT_HEAD(CA_FIELD_EXT_ELEM(K, i));

        if (Fi == CA_Gamma)
        {
            for (j = i + 1; j < len; j++)
            {
                Fj = CA_EXT_HEAD(CA_FIELD_EXT_ELEM(K, j));

                if (Fj == CA_Gamma)
                {
                    ca_ptr xi, xj;
                    ca_t t, u;
                    fmpz_t N;
                    slong n;

                    xi = CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(K, i));
                    xj = CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(K, j));

                    ca_init(t, ctx);
                    ca_init(u, ctx);
                    fmpz_init(N);

                    ca_sub(t, xi, xj, ctx);

                    if (ca_get_fmpz(N, t, ctx) && fmpz_cmp_si(N, -10) >= 0 && fmpz_cmp_si(N, 10) <= 0)
                    {
                        n = fmpz_get_si(N);

                        if (n == 0)
                        {
                            fmpz_mpoly_t poly;
                            fmpz_mpoly_init(poly, CA_FIELD_MCTX(K, ctx));
                            fmpz_mpoly_set_coeff_si_x(poly, 1, i, 1, CA_FIELD_MCTX(K, ctx));
                            fmpz_mpoly_set_coeff_si_x(poly, -1, j, 1, CA_FIELD_MCTX(K, ctx));
                            _ca_field_ideal_insert_clear_mpoly(K, poly, CA_FIELD_MCTX(K, ctx), ctx);
                        }
                        else
                        {
                            ca_field_struct * L;   /* Field of t */
                            slong L_len;
                            slong * tgen_map;
                            slong k, m;
                            int success;

                            /* gamma(x+3) = (x+2)*(x+1)*x * gamma(x) */
                            /* (x-1)*(x-2)*(x-3) * gamma(x-3) = gamma(x) */
                            if (n > 0)
                            {
                                ca_set(t, xj, ctx);
                                for (k = 1; k < n; k++)
                                {
                                    ca_add_ui(u, xj, k, ctx);
                                    ca_mul(t, t, u, ctx);
                                }
                            }
                            else
                            {
                                ca_sub_ui(t, xj, 1, ctx);
                                for (k = 1; k < (-n); k++)
                                {
                                    ca_sub_ui(u, xj, k + 1, ctx);
                                    ca_mul(t, t, u, ctx);
                                }
                            }

                            L = CA_FIELD(t, ctx);
                            L_len = CA_FIELD_LENGTH(L);

                            success = 1;
                            tgen_map = flint_malloc(L_len * sizeof(slong));

                            for (m = 0; m < L_len; m++)
                            {
                                for (k = 0; k < len; k++)
                                {
                                    if (CA_FIELD_EXT_ELEM(L, m) == CA_FIELD_EXT_ELEM(K, k))
                                    {
                                        tgen_map[m] = k;
                                        break;
                                    }

                                    if (k == len - 1)
                                        success = 0;
                                }
                            }

                            if (success)
                            {
                                fmpz_mpoly_t p, q, pxi, pxj;

                                fmpz_mpoly_init(p, CA_FIELD_MCTX(K, ctx));
                                fmpz_mpoly_init(q, CA_FIELD_MCTX(K, ctx));
                                fmpz_mpoly_init(pxi, CA_FIELD_MCTX(K, ctx));
                                fmpz_mpoly_init(pxj, CA_FIELD_MCTX(K, ctx));

                                /* todo: factor out */
                                if (L_len == 0)
                                {
                                    fmpz_mpoly_set_fmpz(p, CA_FMPQ_NUMREF(t), CA_FIELD_MCTX(K, ctx));
                                    fmpz_mpoly_set_fmpz(q, CA_FMPQ_DENREF(t), CA_FIELD_MCTX(K, ctx));
                                }
                                else if (CA_FIELD_IS_NF(L))
                                {
                                    fmpz_poly_t pol;
                                    fmpz_t den;

                                    _nf_elem_get_fmpz_poly_den_shallow(pol, den, CA_NF_ELEM(t), CA_FIELD_NF(L));

                                    fmpz_mpoly_set_gen_fmpz_poly(p, tgen_map[0], pol, CA_FIELD_MCTX(K, ctx));
                                    fmpz_mpoly_set_fmpz(q, den, CA_FIELD_MCTX(K, ctx));
                                }
                                else
                                {
                                    fmpz_mpoly_compose_fmpz_mpoly_gen(p,
                                                              fmpz_mpoly_q_numref(CA_MPOLY_Q(t)),
                                                                tgen_map,
                                                                CA_FIELD_MCTX(L, ctx),
                                                                CA_FIELD_MCTX(K, ctx));
                                    fmpz_mpoly_compose_fmpz_mpoly_gen(q,
                                                              fmpz_mpoly_q_denref(CA_MPOLY_Q(t)),
                                                                tgen_map,
                                                                CA_FIELD_MCTX(L, ctx),
                                                                CA_FIELD_MCTX(K, ctx));
                                }

                                fmpz_mpoly_gen(pxi, i, CA_FIELD_MCTX(K, ctx));
                                fmpz_mpoly_gen(pxj, j, CA_FIELD_MCTX(K, ctx));

                                if (n < 0)
                                    fmpz_mpoly_swap(p, q, CA_FIELD_MCTX(K, ctx));

                                /* pxi = (p/q) * pxj */
                                /* q * pxi - p * pxj = 0*/
                                fmpz_mpoly_mul(q, q, pxi, CA_FIELD_MCTX(K, ctx));
                                fmpz_mpoly_mul(p, p, pxj, CA_FIELD_MCTX(K, ctx));
                                fmpz_mpoly_sub(q, q, p, CA_FIELD_MCTX(K, ctx));

                                _ca_field_ideal_insert_clear_mpoly(K, q, CA_FIELD_MCTX(K, ctx), ctx);

                                fmpz_mpoly_clear(p, CA_FIELD_MCTX(K, ctx));
                                fmpz_mpoly_clear(pxi, CA_FIELD_MCTX(K, ctx));
                                fmpz_mpoly_clear(pxj, CA_FIELD_MCTX(K, ctx));
                            }

                            flint_free(tgen_map);
                        }
                    }

                    ca_clear(t, ctx);
                    ca_clear(u, ctx);
                    fmpz_clear(N);
                }
            }
        }
    }
}
