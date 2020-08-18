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



void
ca_field_build_ideal(ca_field_t K, ca_ctx_t ctx)
{
    slong i, len;

    len = CA_FIELD_LENGTH(K);

    if (len <= 0)
        return;

    if (len == 1 && CA_EXT_IS_QQBAR(CA_FIELD_EXT_ELEM(K, 0)))
        return;

    /* Find direct algebraic relations. */
    if (len >= 2)
    {
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

                    if (L_len == 0)  /* todo: use CA_IS_QQ when renamed */
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

    /* Find log relations. */
    if (len >= 2)
    {
        slong * logs;
        slong num_logs, num_logs_with_pi_i;
        slong prec;
        int have_pi, have_i, have_pi_i;
        slong index_pi, index_i;

        num_logs = 0;
        logs = flint_malloc(sizeof(slong) * len);

        for (i = 0; i < len; i++)
        {
            if (CA_EXT_HEAD(CA_FIELD_EXT_ELEM(K, i)) == CA_Log)
            {
                logs[num_logs] = i;
                num_logs++;
            }
        }

        have_i = have_pi = 0;
        index_i = index_pi = -1;
        for (i = 0; i < len; i++)
        {
            if (CA_FIELD_EXT_ELEM(K, i) == CA_FIELD_EXT_ELEM(ctx->field_qq_i, 0))
            {
                index_i = i;
                have_i = 1;
                break;
            }
        }

        for (i = 0; i < len; i++)
        {
            if (CA_EXT_HEAD(CA_FIELD_EXT_ELEM(K, i)) == CA_Pi)
            {
                index_pi = i;
                have_pi = 1;
                break;
            }
        }

        have_pi_i = have_pi && have_i;
        num_logs_with_pi_i = num_logs + have_pi_i;

        if (num_logs_with_pi_i >= 2)
        {
            acb_ptr z;
            acb_t t;
            mag_t tm;
            slong j, alloc;
            fmpz * rel;
            int found_relation = 0;

            /* todo: dynamic precision determined by context */
            prec = 128;

            alloc = num_logs_with_pi_i;
            z = _acb_vec_init(alloc);
            rel = _fmpz_vec_init(alloc);
            acb_init(t);
            mag_init(tm);

            for (j = 0; j < num_logs; j++)
            {
                /* todo: take advantage of cached enclosure */
                ca_get_acb(z + j, CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(K, logs[j])), prec, ctx);
                acb_log(z + j, z + j, prec);
            }

            /* Add 2 pi i as another logarithm */
            if (have_pi_i)
            {
                acb_const_pi(z + num_logs, prec);
                acb_mul_onei(z + num_logs, z + num_logs);
                acb_mul_2exp_si(z + num_logs, z + num_logs, 1);
            }

            while (num_logs_with_pi_i >= 2)
            {
                found_relation = 0;

                if (_qqbar_acb_lindep(rel, z, num_logs_with_pi_i, 1, prec))
                {
                    ca_t prod, upow;

                    /*  a^m * b^n = 1 => m*log(a) + n*log(b) = 2 pi i k   */
                    /* Verify that (m*log(a) + n*log(b)) / (2 pi i) contains unique integer. */
                    /* It is enough to show that |... + ...| < 2^1. */
                    acb_zero(t);
                    for (j = 0; j < num_logs_with_pi_i; j++)
                        if (!fmpz_is_zero(rel + j))
                            acb_addmul_fmpz(t, z + j, rel + j, prec);
                    acb_get_mag(tm, t);

                    if (mag_cmp_2exp_si(tm, 1) < 0)
                    {
                        ca_init(prod, ctx);
                        ca_init(upow, ctx);

                        /* Check product of arguments to log. */
                        /* Todo: separate positive/negative exponents, thus avoiding inverses? */
                        ca_one(prod, ctx);
                        for (j = 0; j < num_logs; j++)
                        {
                            if (!fmpz_is_zero(rel + j))
                            {
                                ca_pow_fmpz(upow, CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(K, logs[j])), rel + j, ctx);
                                ca_mul(prod, prod, upow, ctx);
                            }
                        }

                        if (ca_check_is_one(prod, ctx) == T_TRUE)
                        {
                            found_relation = 1;

        /*
                            printf("proved log relation!\n");
                            for (j = 0; j < num_logs; j++)
                            {
                                fmpz_print(rel + j); printf(" ");
                                acb_printn(z + j, 10, 0); printf("   ");
                            }
                            printf("\n");

                            for (j = 0; j < num_logs; j++)
                            {
                                if (!fmpz_is_zero(rel + j) || 1)
                                {
                                    fmpz_print(rel + j);
                                    flint_printf(" * log(");
                                    ca_print(CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(K, logs[j])), ctx);
                                    flint_printf(")    ");
                                }
                            }
                            flint_printf("\n");
    */

                            /* todo: some kind of ideal fit_length method... */
                            if (CA_FIELD_IDEAL_LENGTH(K) == 0)
                                CA_FIELD_IDEAL(K) = flint_malloc(sizeof(fmpz_mpoly_struct));
                            else
                                CA_FIELD_IDEAL(K) = flint_realloc(CA_FIELD_IDEAL(K), (CA_FIELD_IDEAL_LENGTH(K) + 1) * sizeof(fmpz_mpoly_struct));

                            fmpz_mpoly_init(CA_FIELD_IDEAL_ELEM(K, CA_FIELD_IDEAL_LENGTH(K)), CA_FIELD_MCTX(K, ctx));

                            {
                                ulong * exp;

                                exp = flint_malloc(sizeof(ulong) * len);

                                for (j = 0; j < num_logs_with_pi_i; j++)
                                {
                                    slong k;

                                    if (fmpz_is_zero(rel + j))
                                        continue;

                                    for (k = 0; k < len; k++)
                                        exp[k] = 0;

                                    /* 2 pi i */
                                    if (j == num_logs)
                                    {
                                        exp[index_i] = 1;
                                        exp[index_pi] = 1;
                                        fmpz_mul_2exp(rel + j, rel + j, 1);
                                    }
                                    else
                                    {
                                        exp[logs[j]] = 1;
                                    }

                                    fmpz_mpoly_set_coeff_fmpz_ui(CA_FIELD_IDEAL_ELEM(K, CA_FIELD_IDEAL_LENGTH(K)),
                                        rel + j,
                                        exp,
                                        CA_FIELD_MCTX(K, ctx));
                                }

                                flint_free(exp);
                            }

                            CA_FIELD_IDEAL_LENGTH(K)++;
                        }

                        ca_clear(prod, ctx);
                        ca_clear(upow, ctx);
                    }
                }

                if (!found_relation)
                    break;

                for (j = 0; j < num_logs - 1; j++)
                    logs[j] = logs[j + 1];

                for (j = 0; j < num_logs_with_pi_i - 1; j++)
                    acb_swap(z + j, z + j + 1);

                num_logs--;
                num_logs_with_pi_i--;
            }

            _acb_vec_clear(z, alloc);
            _fmpz_vec_clear(rel, alloc);
            acb_clear(t);
            mag_clear(tm);
        }

        flint_free(logs);
    }

    /* Find relations involving exponentials, powers and roots. */
    /* Todo: involve algebraic numbers (including roots of unity) */
    if (len >= 1)
    {
        slong * powers;
        slong num_powers;
        slong prec;

        num_powers = 0;
        powers = flint_malloc(sizeof(slong) * len);

        for (i = 0; i < len; i++)
        {
            if (CA_EXT_HEAD(CA_FIELD_EXT_ELEM(K, i)) == CA_Sqrt)
            {
                powers[num_powers] = i;
                num_powers++;
            }
            else if (CA_EXT_HEAD(CA_FIELD_EXT_ELEM(K, i)) == CA_Pow)
            {
                powers[num_powers] = i;
                num_powers++;
            }
            else if (CA_EXT_HEAD(CA_FIELD_EXT_ELEM(K, i)) == CA_Exp)
            {
                powers[num_powers] = i;
                num_powers++;
            }
            else if (CA_EXT_IS_QQBAR(CA_FIELD_EXT_ELEM(K, i)))
            {
                powers[num_powers] = i;
                num_powers++;
            }
        }

        if (num_powers >= 1)
        {
            acb_ptr z;
            fmpz * rel;
            slong alloc, j;
            int found_relation = 0;

            /* todo: dynamic precision determined by context */
            prec = 128;

            alloc = num_powers + 1;

            z = _acb_vec_init(alloc);
            rel = _fmpz_vec_init(alloc);

            for (i = 0; i < num_powers; i++)
            {
                ca_ext_srcptr ext = CA_FIELD_EXT_ELEM(K, powers[i]);

                if (CA_EXT_HEAD(ext) == CA_Sqrt)
                {
                    ca_get_acb(z + i, CA_EXT_FUNC_ARGS(ext), prec, ctx);
                    acb_log(z + i, z + i, prec);
                    acb_mul_2exp_si(z + i, z + i, -1);
                }
                else if (CA_EXT_HEAD(ext) == CA_Pow)
                {
                    ca_get_acb(z + i, CA_EXT_FUNC_ARGS(ext), prec, ctx);
                    acb_log(z + i, z + i, prec);
                    ca_get_acb(z + i + 1, CA_EXT_FUNC_ARGS(ext) + 1, prec, ctx);   /* xxx */
                    acb_mul(z + i, z + i, z + i + 1, prec);
                }
                else if (CA_EXT_HEAD(ext) == CA_Exp)
                {
                    ca_get_acb(z + i, CA_EXT_FUNC_ARGS(ext), prec, ctx);
                }
                else if (CA_EXT_IS_QQBAR(ext))
                {
                    qqbar_get_acb(z + i, CA_EXT_QQBAR(ext), prec);
                    acb_log(z + i, z + i, prec);
                }
                else
                {
                    flint_abort();
                }
            }

            /* todo: exp(pi i) = -1 --- look for other roots of unity! */
            /* todo: what else to include -- log(2), log(3), log(5) ? */

            acb_const_pi(z + num_powers, prec);
            acb_mul_onei(z + num_powers, z + num_powers);

            /* todo: verify that logs are not evaluated at zero... */

            while (num_powers >= 1)
            {
                found_relation = 0;

/*
                for (i = 0; i < num_powers + 1; i++)
                {
                    acb_printn(z + i, 30, 0); printf("\n");
                }
*/

                if (_qqbar_acb_lindep(rel, z, num_powers + 1, 1, prec) && FLINT_ABS(_fmpz_vec_max_bits(rel, num_powers + 1)) <= 10)
                {
                    ca_t t, u;
                    ca_init(t, ctx);
                    ca_init(u, ctx);

/*
                    flint_printf("Candidate: ");
                    for (i = 0; i < num_powers + 1; i++)
                    {
                        fmpz_print(rel + i); printf(" * ");

                        if (i == num_powers)
                            printf("Pi*I  ");
                        else
                        {
                            ca_ext_print(CA_FIELD_EXT_ELEM(K, powers[i]), ctx); flint_printf("  ");
                        }
                    }
                    flint_printf("\n\n");
*/

                    for (i = 0; i < num_powers; i++)
                    {
                        if (!fmpz_is_zero(rel + i))
                        {
                            ca_ext_srcptr ext = CA_FIELD_EXT_ELEM(K, powers[i]);

                            if (CA_EXT_HEAD(ext) == CA_Sqrt)
                            {
                                ca_log(u, CA_EXT_FUNC_ARGS(ext), ctx);
                                ca_div_ui(u, u, 2, ctx);
                            }
                            else if (CA_EXT_HEAD(ext) == CA_Pow)
                            {
                                ca_log(u, CA_EXT_FUNC_ARGS(ext), ctx);
                                ca_mul(u, u, CA_EXT_FUNC_ARGS(ext) + 1, ctx);
                            }
                            else if (CA_EXT_HEAD(ext) == CA_Exp)
                            {
                                ca_set(u, CA_EXT_FUNC_ARGS(ext), ctx);
                            }
                            else if (CA_EXT_IS_QQBAR(ext))
                            {
                                ca_set_qqbar(u, CA_EXT_QQBAR(ext), ctx);
                                ca_log(u, u, ctx);
                            }
                            else
                            {
                                flint_abort();
                            }

                            ca_mul_fmpz(u, u, rel + i, ctx);
                            ca_add(t, t, u, ctx);
                        }
                    }

                    if (!fmpz_is_zero(rel + num_powers))
                    {
                        ca_pi_i(u, ctx);
                        ca_mul_fmpz(u, u, rel + num_powers, ctx);
                        ca_add(t, t, u, ctx);
                    }

/*
                    printf("prove: "); ca_print(t, ctx); printf("\n");
*/

                    if (ca_check_is_zero(t, ctx) == T_TRUE)
                    {
                        found_relation = 1;

/*
                        flint_printf("Found: ");
                        for (i = 0; i < num_powers + 1; i++)
                        {
                            fmpz_print(rel + i); printf(" * ");

                            if (i == num_powers)
                                printf("Pi*I  ");
                            else
                            {
                                ca_ext_print(CA_FIELD_EXT_ELEM(K, powers[i]), ctx); flint_printf("  ");
                            }
                        }
                        flint_printf("\n\n");
*/

                        /* todo: some kind of ideal fit_length method... */
                        if (CA_FIELD_IDEAL_LENGTH(K) == 0)
                            CA_FIELD_IDEAL(K) = flint_malloc(sizeof(fmpz_mpoly_struct));
                        else
                            CA_FIELD_IDEAL(K) = flint_realloc(CA_FIELD_IDEAL(K), (CA_FIELD_IDEAL_LENGTH(K) + 1) * sizeof(fmpz_mpoly_struct));

                        fmpz_mpoly_init(CA_FIELD_IDEAL_ELEM(K, CA_FIELD_IDEAL_LENGTH(K)), CA_FIELD_MCTX(K, ctx));

                        /* x^a y^b +/- z^d w^e = 0 */

                        {
                            ulong * exp1;
                            ulong * exp2;
                            int neg;

                            exp1 = flint_calloc(len, sizeof(ulong));
                            exp2 = flint_calloc(len, sizeof(ulong));

                            neg = fmpz_is_odd(rel + num_powers);

                            for (j = 0; j < num_powers; j++)
                            {
                                if (fmpz_is_zero(rel + j))
                                    continue;

                                if (fmpz_sgn(rel + j) > 0)
                                    exp1[powers[j]] = rel[j];
                                else
                                    exp2[powers[j]] = -rel[j];
                            }

                            /* todo: normalise sign? */
                            fmpz_mpoly_set_coeff_si_ui(CA_FIELD_IDEAL_ELEM(K, CA_FIELD_IDEAL_LENGTH(K)),
                                1,
                                exp1,
                                CA_FIELD_MCTX(K, ctx));

                            fmpz_mpoly_set_coeff_si_ui(CA_FIELD_IDEAL_ELEM(K, CA_FIELD_IDEAL_LENGTH(K)),
                                neg ? 1 : -1,
                                exp2,
                                CA_FIELD_MCTX(K, ctx));

                            flint_free(exp1);
                            flint_free(exp2);
                        }

                        CA_FIELD_IDEAL_LENGTH(K)++;
                    }

                    ca_clear(t, ctx);
                    ca_clear(u, ctx);
                }

                if (!found_relation)
                    break;

                for (j = 0; j < num_powers - 1; j++)
                    powers[j] = powers[j + 1];

                for (j = 0; j < num_powers; j++)
                    acb_swap(z + j, z + j + 1);

                num_powers--;
            }

            _acb_vec_clear(z, alloc);
            _fmpz_vec_clear(rel, alloc);
        }

        flint_free(powers);
    }
}
