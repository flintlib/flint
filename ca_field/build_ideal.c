/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"
#include "ca_ext.h"
#include "ca_field.h"

void _nf_elem_get_fmpz_poly_den_shallow(fmpz_poly_t pol, fmpz_t den, const nf_elem_t a, const nf_t nf);
void fmpz_mpoly_set_gen_fmpz_poly(fmpz_mpoly_t res, slong var, const fmpz_poly_t pol, const fmpz_mpoly_ctx_t ctx);


        /* Find log relations. Todo: this needs to be done before building
           the final field, because we may be introducing extra elements
           (pi, i). */
#if 0
        if (0)
        {
            slong * logs;
            slong num_logs;

            num_logs = 0;
            logs = flint_malloc(sizeof(slong) * fields_len);

            /* todo: find linear combinations of logarithms */
            for (i = 0; i < fields_len; i++)
            {
                if ((ctx->fields + fields[i])->type == CA_FIELD_TYPE_FUNC &&
                    (ctx->fields + fields[i])->data.func.func == CA_Log)
                {
                    logs[num_logs] = i;
                    num_logs++;
                }
            }

            if (num_logs >= 2)
            {
                acb_ptr z;
                slong j;
                fmpz * rel;

                z = _acb_vec_init(num_logs + 1);
                rel = _fmpz_vec_init(num_logs + 1);
                /* todo: pi * i */

                for (j = 0; j < num_logs; j++)
                {
                    ca_get_acb(z + j, (ctx->fields + fields[logs[j]])->data.func.args, 256, ctx);
                    acb_log(z + j, z + j, 256);
                }

                if (_qqbar_acb_lindep(rel, z, num_logs, 1, 256))
                {
                    ca_t prod, upow;

                    ca_init(prod, ctx);
                    ca_init(upow, ctx);

                    ca_one(prod, ctx);

                    for (j = 0; j < num_logs; j++)
                    {
                        if (!fmpz_is_zero(rel + j))
                        {
                            ca_pow_fmpz(upow, (ctx->fields + fields[logs[j]])->data.func.args, rel + j, ctx);

/*
                            ca_print((ctx->fields + fields[logs[j]])->data.func.args, ctx); printf(" ^ "); fmpz_print(rel + j); printf(" = "); ca_print(upow, ctx); printf("\n");
*/

                            ca_mul(prod, prod, upow, ctx);
                        }
                    }

                    if (ca_check_is_one(prod, ctx) == T_TRUE)
                    {
                        printf("proved log relation!\n");

/*
                        for (j = 0; j < num_logs; j++)
                        {
                            acb_printn(z + j, 10, 0); printf("   ");
                        }
                        printf("\n");
*/

                        for (j = 0; j < num_logs; j++)
                        {
                            if (!fmpz_is_zero(rel + j) || 1)
                            {
                                fmpz_print(rel + j);
                                flint_printf(" * log(");
                                ca_print((ctx->fields + fields[logs[j]])->data.func.args, ctx);
                                flint_printf(")    ");
                            }
                        }
                        flint_printf("\n");
                    }

                    ca_clear(prod, ctx);
                    ca_clear(upow, ctx);
                }

                _acb_vec_clear(z, num_logs + 1);
                _fmpz_vec_clear(rel, num_logs + 1);
            }

            flint_free(logs);
        }
#endif

void
ca_field_build_ideal(ca_field_t K, ca_ctx_t ctx)
{
    slong i, len;

    len = CA_FIELD_LENGTH(K);

    if (len <= 1)
        return;

    for (i = 0; i < len; i++)
    {
        ca_ext_struct * x = CA_FIELD_EXT_ELEM(K, i);

        /* Absolute annihilating polynomial for number field */
        if (CA_EXT_IS_QQBAR(x))
        {
            if (CA_FIELD_IDEAL_LENGTH(K) == 0)
                CA_FIELD_IDEAL(K) = flint_malloc(sizeof(fmpz_mpoly_struct));
            else
                CA_FIELD_IDEAL(K) = flint_realloc(CA_FIELD_IDEAL(K), (CA_FIELD_IDEAL_LENGTH(K) + 1) * sizeof(fmpz_mpoly_struct));

            fmpz_mpoly_init(CA_FIELD_IDEAL(K) + CA_FIELD_IDEAL_LENGTH(K), CA_FIELD_MCTX(K, ctx));
            fmpz_mpoly_set_gen_fmpz_poly(CA_FIELD_IDEAL(K) + CA_FIELD_IDEAL_LENGTH(K), i, QQBAR_POLY(CA_EXT_QQBAR(x)), CA_FIELD_MCTX(K, ctx));

            CA_FIELD_IDEAL_LENGTH(K)++;
            continue;
        }

        /* x = sqrt(t) -> x^2 - t = 0 */
        /* A relation will only be added if t can be expressed in the
           present field. */
        if (CA_EXT_HEAD(x) == CA_Sqrt)
        {
            ca_srcptr t;
            ca_field_struct * L;   /* Field of t */
            slong L_len;
            slong * tgen_map;
            slong j, k;
            int success;

            if (CA_EXT_FUNC_NARGS(x) != 1)
                flint_abort();

            t = CA_EXT_FUNC_ARGS(x);

            if (CA_IS_SPECIAL(t))
                flint_abort();

            L = CA_FIELD(t, ctx);
            L_len = CA_FIELD_LENGTH(L);

            success = 1;
            tgen_map = flint_malloc(L_len * sizeof(slong));

            for (j = 0; j < L_len; j++)
            {
                for (k = 0; k < len; k++)
                {
                    if (CA_FIELD_EXT_ELEM(L, j) == CA_FIELD_EXT_ELEM(K, k))
                    {
                        tgen_map[j] = k;
                        break;
                    }

                    if (k == len - 1)
                        success = 0;
                }
            }

            if (success)
            {
                /* u^2 - p/q  -->  q u^2 - p */
                fmpz_mpoly_t p, q, u2;

                fmpz_mpoly_init(p, CA_FIELD_MCTX(K, ctx));
                fmpz_mpoly_init(q, CA_FIELD_MCTX(K, ctx));
                fmpz_mpoly_init(u2, CA_FIELD_MCTX(K, ctx));

                if (L_len == 0)  /* todo: use CA_FIELD_IS_QQ when renamed */
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

                fmpz_mpoly_gen(u2, i, CA_FIELD_MCTX(K, ctx));
                fmpz_mpoly_pow_ui(u2, u2, 2, CA_FIELD_MCTX(K, ctx));

                fmpz_mpoly_mul(u2, u2, q, CA_FIELD_MCTX(K, ctx));
                fmpz_mpoly_sub(u2, u2, p, CA_FIELD_MCTX(K, ctx));

                /* todo: some kind of ideal fit_length method... */
                if (CA_FIELD_IDEAL_LENGTH(K) == 0)
                    CA_FIELD_IDEAL(K) = flint_malloc(sizeof(fmpz_mpoly_struct));
                else
                    CA_FIELD_IDEAL(K) = flint_realloc(CA_FIELD_IDEAL(K), (CA_FIELD_IDEAL_LENGTH(K) + 1) * sizeof(fmpz_mpoly_struct));

                /* todo: avoid a copy */
                fmpz_mpoly_init(CA_FIELD_IDEAL_ELEM(K, CA_FIELD_IDEAL_LENGTH(K)), CA_FIELD_MCTX(K, ctx));
                fmpz_mpoly_set(CA_FIELD_IDEAL_ELEM(K, CA_FIELD_IDEAL_LENGTH(K)), u2, CA_FIELD_MCTX(K, ctx));

                CA_FIELD_IDEAL_LENGTH(K)++;

                fmpz_mpoly_clear(p, CA_FIELD_MCTX(K, ctx));
                fmpz_mpoly_clear(q, CA_FIELD_MCTX(K, ctx));
                fmpz_mpoly_clear(u2, CA_FIELD_MCTX(K, ctx));
            }

            flint_free(tgen_map);
        }
    }
}
