/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "ca.h"
#include "ca_ext.h"
#include "ca_field.h"

#include "fmpz_mat.h"
#include "fmpz_lll.h"
#include "qqbar.h"
#include "fmpz_mpoly.h"

int qqbar_mul_checked(qqbar_t res, const qqbar_t x, const qqbar_t y, slong deg_limit, slong bits_limit);
int qqbar_pow_fmpz_checked(qqbar_t res, const qqbar_t x, const fmpz_t y, slong deg_limit, slong bits_limit);

slong
acb_multi_lindep(fmpz_mat_t rel, acb_srcptr vec, slong len, int check, slong prec)
{
    fmpz_mat_t A;
    arf_t tmpr, halfr;
    fmpz_lll_t ctx;
    fmpz_t scale_exp;
    int nonreal, found;
    slong i, row, accuracy;
    acb_t z2;
    mag_t max_size, max_rad, tmpmag;

    if (fmpz_mat_nrows(rel) != 0 || fmpz_mat_ncols(rel) != 0)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    fmpz_mat_clear(rel);

    for (i = 0; i < len; i++)
    {
        if (!acb_is_finite(vec + i))
        {
            fmpz_mat_init(rel, 0, len);
            return 0;
        }
    }

    found = 0;

    nonreal = 0;
    for (i = 0; i < len; i++)
        if (!arb_contains_zero(acb_imagref(vec + i)))
            nonreal = 1;

    fmpz_mat_init(A, len, len + 1 + nonreal);
    fmpz_init(scale_exp);
    acb_init(z2);
    arf_init(tmpr);
    arf_init(halfr);
    mag_init(max_size);
    mag_init(max_rad);
    mag_init(tmpmag);

    arf_set_d(halfr, 0.5);

    for (i = 0; i < len; i++)
    {
        arf_get_mag(tmpmag, arb_midref(acb_realref(vec + i)));
        mag_max(max_size, max_size, tmpmag);
        arf_get_mag(tmpmag, arb_midref(acb_imagref(vec + i)));
        mag_max(max_size, max_size, tmpmag);

        mag_max(max_rad, max_rad, arb_radref(acb_realref(vec + i)));
        mag_max(max_rad, max_rad, arb_radref(acb_imagref(vec + i)));
    }

    prec = FLINT_MAX(prec, 2);
    if (!mag_is_zero(max_size) && !mag_is_zero(max_rad))
    {
        accuracy = _fmpz_sub_small(MAG_EXPREF(max_size), MAG_EXPREF(max_rad));
        accuracy = FLINT_MAX(accuracy, 10);
        prec = FLINT_MIN(prec, accuracy);
    }

    if (mag_is_zero(max_size))
    {
        fmpz_zero(scale_exp);  /* todo: quick return? */
    }
    else
    {
        fmpz_neg(scale_exp, MAG_EXPREF(max_size));
        fmpz_add_ui(scale_exp, scale_exp, prec);
    }

    /* Using 5% of the bits for checking will provide some protection
       against spurious relations */
    fmpz_sub_ui(scale_exp, scale_exp, FLINT_MAX(10, prec * 0.05));

    /* Create matrix */
    for (i = 0; i < len; i++)
        fmpz_one(fmpz_mat_entry(A, i, i));

    for (i = 0; i < len; i++)
    {
        arf_mul_2exp_fmpz(tmpr, arb_midref(acb_realref(vec + i)), scale_exp);
        arf_add(tmpr, tmpr, halfr, prec, ARF_RND_NEAR);
        arf_floor(tmpr, tmpr);
        arf_get_fmpz(fmpz_mat_entry(A, i, len), tmpr, ARF_RND_NEAR);

        if (nonreal)
        {
            arf_mul_2exp_fmpz(tmpr, arb_midref(acb_imagref(vec + i)), scale_exp);
            arf_add(tmpr, tmpr, halfr, prec, ARF_RND_NEAR);
            arf_floor(tmpr, tmpr);
            arf_get_fmpz(fmpz_mat_entry(A, i, len + 1), tmpr, ARF_RND_NEAR);
        }
    }

    /* LLL reduction */
    fmpz_lll_context_init(ctx, 0.75, 0.51, 1, 0);
#if 1
    fmpz_lll(A, NULL, ctx);
#else
    {
        fmpz_t gb;
        fmpz_init(gb);
        fmpz_one(gb);
        fmpz_mul_2exp(gb, gb, prec / 2);
        fmpz_lll_with_removal(A, NULL, gb, ctx);
        fmpz_clear(gb);
    }
#endif

    /* Heuristic check */
    for (row = 0; row < len; row++)
    {
        acb_zero(z2);
        for (i = 0; i < len; i++)
            acb_addmul_fmpz(z2, vec + i, fmpz_mat_entry(A, row, i), prec + 10);

        if (!_fmpz_vec_is_zero(A->rows[row], len) && acb_contains_zero(z2))
        {
            found++;
        }
        else
        {
            _fmpz_vec_zero(A->rows[row], fmpz_mat_ncols(A));
        }
    }

    fmpz_mat_init(rel, found, len);

    i = 0;
    for (row = 0; row < len; row++)
    {
        if (!_fmpz_vec_is_zero(A->rows[row], len))
        {
            _fmpz_vec_set(rel->rows[i], A->rows[row], len);
            i++;
        }
    }

    /* todo: move outside? */
    if (found > 1)
    {
        fmpz_mat_hnf(rel, rel);
    }

    fmpz_mat_clear(A);
    fmpz_clear(scale_exp);
    acb_clear(z2);
    arf_clear(tmpr);
    arf_clear(halfr);
    mag_clear(max_size);
    mag_clear(max_rad);
    mag_clear(tmpmag);

    /* fixme: bogus */
    return found;
}

void _nf_elem_get_fmpz_poly_den_shallow(fmpz_poly_t pol, fmpz_t den, const nf_elem_t a, const nf_t nf);

void
_ca_field_ideal_insert_clear_mpoly(ca_field_t K, fmpz_mpoly_t poly, fmpz_mpoly_ctx_t mctx, ca_ctx_t ctx)
{
    if (poly->length == 0)
    {
        flint_throw(FLINT_ERROR, "ERROR: inserting the zero polynomial into ideal\n");
    }

    if (fmpz_sgn(poly->coeffs) < 0)
        fmpz_mpoly_neg(poly, poly, mctx);

    /* todo: move!!! */
    fmpz_mpoly_vec_insert_unique(CA_FIELD_IDEAL(K), poly, mctx);
    fmpz_mpoly_clear(poly, mctx);
}

int ext_as_pow_pq(slong *p, slong *q, const ca_ext_t x, ca_ctx_t ctx)
{
    if (CA_EXT_HEAD(x) == CA_Sqrt)
    {
        *p = 1;
        *q = 2;
        return 1;
    }

    if (CA_EXT_HEAD(x) == CA_Pow && CA_IS_QQ(CA_EXT_FUNC_ARGS(x) + 1, ctx))
    {
        fmpz pp, qq;

        pp = *CA_FMPQ_NUMREF(CA_EXT_FUNC_ARGS(x) + 1);
        qq = *CA_FMPQ_DENREF(CA_EXT_FUNC_ARGS(x) + 1);

        if (fmpz_bits(&pp) <= 6 && fmpz_bits(&qq) <= 6)
        {
            *p = fmpz_get_si(&pp);
            *q = fmpz_get_si(&qq);
            return 1;
        }
    }

    return 0;
}

int
ca_field_prove_log_relation(ca_field_t K, const fmpz * rel,
    acb_srcptr z,
    const slong * logs,
    slong num_logs, slong num_logs_with_pi_i, slong prec, ca_ctx_t ctx)
{
    ca_t prod, upow;
    slong j;
    acb_t t;
    mag_t tm;
    int success = 0;

    /*
    printf("possible log relation!\n");
    for (j = 0; j < num_logs; j++)
    {
        fmpz_print(rel + j); printf(" ");
        acb_printn(z + j, 10, 0); printf("   ");
    }
    printf("\n");
    */

    acb_init(t);
    mag_init(tm);

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
        /* Or: otherwise balance the product. */
        ca_one(prod, ctx);
        for (j = 0; j < num_logs; j++)
        {
            if (!fmpz_is_zero(rel + j))
            {
                ca_pow_fmpz(upow, CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(K, logs[j])), rel + j, ctx);
                ca_mul(prod, prod, upow, ctx);
            }
        }

        /* printf("product: "); ca_print(prod, ctx); printf("\n\n"); */

        success = (ca_check_is_one(prod, ctx) == T_TRUE);

        ca_clear(prod, ctx);
        ca_clear(upow, ctx);
    }

    if (0 && success)
    {
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
    }

    acb_clear(t);
    mag_clear(tm);

    return success;
}

slong
ca_field_insert_log_relation(ca_field_t K,
    fmpz * rel,
    const slong * logs,
    slong index_i,
    slong index_pi,
    slong num_logs, slong num_logs_with_pi_i, ca_ctx_t ctx)
{
    fmpz_mpoly_t poly;
    ulong * exp;
    slong which_removed;
    slong j, len;

    len = CA_FIELD_LENGTH(K);
    exp = flint_malloc(sizeof(ulong) * len);

    fmpz_mpoly_init(poly, CA_FIELD_MCTX(K, ctx));

    which_removed = -1;

    for (j = 0; j < num_logs_with_pi_i; j++)
    {
        slong k;

        if (fmpz_is_zero(rel + j))
            continue;

        if (which_removed == -1)
            which_removed = j;

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

        fmpz_mpoly_set_coeff_fmpz_ui(poly,
            rel + j,
            exp,
            CA_FIELD_MCTX(K, ctx));
    }

    flint_free(exp);

    _ca_field_ideal_insert_clear_mpoly(K, poly, CA_FIELD_MCTX(K, ctx), ctx);

    return which_removed;
}


/* Find log relations. */
void
ca_field_build_ideal_logs(ca_field_t K, ca_ctx_t ctx)
{
    slong * logs;
    slong num_logs, num_logs_with_pi_i;
    slong prec;
    int have_pi, have_i, have_pi_i;
    slong index_pi, index_i;
    slong i, len;

    len = CA_FIELD_LENGTH(K);

    if (len < 2)
        return;

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
        fmpz_mat_t A;
        acb_ptr z;
        slong j, alloc;
        fmpz * rel;
        int found_relation = 0;
        slong which_removed = 0;
        slong row;

        prec = ctx->options[CA_OPT_LLL_PREC];

        alloc = num_logs_with_pi_i;
        z = _acb_vec_init(alloc);

        for (j = 0; j < num_logs; j++)
        {
            /* todo: take advantage of cached enclosure */
            ca_get_acb(z + j, CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(K, logs[j])), prec, ctx);
            acb_log(z + j, z + j, prec);
        }

        /* Add 2*pi*i as another logarithm */
        if (have_pi_i)
        {
            acb_const_pi(z + num_logs, prec);
            acb_mul_onei(z + num_logs, z + num_logs);
            acb_mul_2exp_si(z + num_logs, z + num_logs, 1);
        }

        if (1)   /* Find all integer relations at once. */
        {
            fmpz_mat_init(A, 0, 0);
            acb_multi_lindep(A, z, num_logs_with_pi_i, 1, prec);

            for (row = 0; row < fmpz_mat_nrows(A); row++)
            {
                rel = A->rows[row];

                if (!_fmpz_vec_is_zero(rel, num_logs_with_pi_i) && FLINT_ABS(_fmpz_vec_max_bits(rel, num_logs_with_pi_i)) <= 10)
                {
                    if (ca_field_prove_log_relation(K, rel, z, logs, num_logs, num_logs_with_pi_i, prec, ctx))
                    {
                        found_relation = 1;
                        which_removed = ca_field_insert_log_relation(K, rel, logs, index_i, index_pi, num_logs, num_logs_with_pi_i, ctx);
                    }
                }
            }

            fmpz_mat_clear(A);
        }
        else
        {
            rel = _fmpz_vec_init(alloc);

            while (num_logs_with_pi_i >= 2)
            {
                found_relation = 0;

                if (_qqbar_acb_lindep(rel, z, num_logs_with_pi_i, 1, prec))
                {
                    if (ca_field_prove_log_relation(K, rel, z, logs, num_logs, num_logs_with_pi_i, prec, ctx))
                    {
                        found_relation = 1;
                        which_removed = ca_field_insert_log_relation(K, rel, logs, index_i, index_pi, num_logs, num_logs_with_pi_i, ctx);
                    }
                }

                if (!found_relation)
                    break;

                for (j = which_removed; j < num_logs - 1; j++)
                    logs[j] = logs[j + 1];

                for (j = which_removed; j < num_logs_with_pi_i - 1; j++)
                    acb_swap(z + j, z + j + 1);

                num_logs--;
                num_logs_with_pi_i--;
            }

            _fmpz_vec_clear(rel, alloc);
        }

        _acb_vec_clear(z, alloc);
    }

    flint_free(logs);
}

int
ca_field_prove_multiplicative_relation(ca_field_t K, const fmpz * rel,
    acb_srcptr z,
    const slong * powers,
    slong num_powers, slong prec, ca_ctx_t ctx)
{
    ca_t t, u;
    slong i;
    int success = 0;
    int all_qqbar = 1;

    ca_init(t, ctx);
    ca_init(u, ctx);

    if (ctx->options[CA_OPT_VERBOSE])
    {
        flint_printf("Attempt to prove multiplicative relation:\n");
        for (i = 0; i < num_powers + 1; i++)
        {
            flint_printf("    [ ^");
            fmpz_print(rel + i);
            flint_printf("] ");

            if (i == num_powers)
                printf("(-1)  ");
            else
            {
                ca_ext_print(CA_FIELD_EXT_ELEM(K, powers[i]), ctx); flint_printf("  ");
            }
            flint_printf("\n");
        }
        flint_printf("\n");
    }

    /* Special-case all qqbars; todo -- separate qqbars, exps; better algorithm */
    for (i = 0; i < num_powers && all_qqbar; i++)
    {
        if (!fmpz_is_zero(rel + i))
        {
            ca_ext_srcptr ext = CA_FIELD_EXT_ELEM(K, powers[i]);

            if (!CA_EXT_IS_QQBAR(ext))
                all_qqbar = 0;
        }
    }

    if (all_qqbar)
    {
        qqbar_t a, b;

        qqbar_init(a);
        qqbar_init(b);

        qqbar_one(a);

        for (i = 0; i < num_powers; i++)
        {
            if (!fmpz_is_zero(rel + i))
            {
                ca_ext_srcptr ext = CA_FIELD_EXT_ELEM(K, powers[i]);

                /* xxx: bogus limits */
                if (!qqbar_pow_fmpz_checked(b, CA_EXT_QQBAR(ext), rel + i, ctx->options[CA_OPT_QQBAR_DEG_LIMIT], ctx->options[CA_OPT_PREC_LIMIT] * 10))
                {
                    success = 0;
                    goto qqbar_end;
                }

                if (!qqbar_mul_checked(a, a, b, ctx->options[CA_OPT_QQBAR_DEG_LIMIT], ctx->options[CA_OPT_PREC_LIMIT] * 10))
                {
                    success = 0;
                    goto qqbar_end;
                }
            }
        }

        /* (-1)^ */
        if (fmpz_is_odd(rel + num_powers))
            qqbar_neg(a, a);

        success = qqbar_is_one(a);

qqbar_end:

/*
        if (!success)
        {
            printf("qqbar failed!\n");
        }
*/

        qqbar_clear(a);
        qqbar_clear(b);
    }
    else
    {
        /* Todo: don't duplicate these computations */
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
                    flint_throw(FLINT_ERROR, "(%s)\n", __func__);
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

        success = (ca_check_is_zero(t, ctx) == T_TRUE);
    }

    if (ctx->options[CA_OPT_VERBOSE])
    {
        flint_printf("    Success = %d\n\n", success);
    }

    ca_clear(t, ctx);
    ca_clear(u, ctx);

    return success;
}

slong
ca_field_insert_multiplicative_relation(ca_field_t K,
    fmpz * rel,
    const slong * powers,
    slong num_powers, ca_ctx_t ctx)
{
    /* x^a y^b +/- z^d w^e = 0 */
    fmpz_mpoly_t poly;
    ulong * exp1;
    ulong * exp2;
    int neg;
    slong len, j, which_removed;

    len = CA_FIELD_LENGTH(K);
    which_removed = -1;

    fmpz_mpoly_init(poly, CA_FIELD_MCTX(K, ctx));

    exp1 = flint_calloc(len, sizeof(ulong));
    exp2 = flint_calloc(len, sizeof(ulong));

    neg = fmpz_is_odd(rel + num_powers);

    for (j = 0; j < num_powers; j++)
    {
        if (fmpz_is_zero(rel + j))
            continue;

        if (which_removed == -1)
            which_removed = j;

        if (fmpz_sgn(rel + j) > 0)
            exp1[powers[j]] = rel[j];
        else
            exp2[powers[j]] = -rel[j];
    }

    /* todo: normalise sign? */
    fmpz_mpoly_set_coeff_si_ui(poly,
        1,
        exp1,
        CA_FIELD_MCTX(K, ctx));

    fmpz_mpoly_set_coeff_si_ui(poly,
        neg ? 1 : -1,
        exp2,
        CA_FIELD_MCTX(K, ctx));

    flint_free(exp1);
    flint_free(exp2);

    _ca_field_ideal_insert_clear_mpoly(K, poly, CA_FIELD_MCTX(K, ctx), ctx);

    return which_removed;
}


/* Find relations involving exponentials, powers and roots. */
/* Todo: more special cases (including roots of unity) */
void
ca_field_build_ideal_multiplicative(ca_field_t K, ca_ctx_t ctx)
{
    slong i, len;
    slong * powers;
    slong num_powers;
    slong prec;

    len = CA_FIELD_LENGTH(K);

    if (len == 0)
        return;

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
        slong which_removed = 0;
        fmpz_mat_t A;
        slong row;

        prec = ctx->options[CA_OPT_LLL_PREC];

        alloc = num_powers + 1;

        z = _acb_vec_init(alloc);

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
                flint_throw(FLINT_ERROR, "(%s)\n", __func__);
            }
        }

        /* todo: exp(pi i) = -1 --- look for other roots of unity! */
        /* todo: what else to include -- log(2), log(3), log(5) ? */

        acb_const_pi(z + num_powers, prec);
        acb_mul_onei(z + num_powers, z + num_powers);

        /* todo: verify that logs are not evaluated at zero... */

        if (1)   /* Find all integer relations at once. */
        {
            fmpz_mat_init(A, 0, 0);
            acb_multi_lindep(A, z, num_powers + 1, 1, prec);

            for (row = 0; row < fmpz_mat_nrows(A); row++)
            {
                rel = A->rows[row];

                if (!_fmpz_vec_is_zero(rel, num_powers + 1) && FLINT_ABS(_fmpz_vec_max_bits(rel, num_powers + 1)) <= 10)
                {
                    if (ca_field_prove_multiplicative_relation(K, rel, z, powers, num_powers, prec, ctx))
                    {
                        ca_field_insert_multiplicative_relation(K, rel, powers, num_powers, ctx);
                    }
                }
            }

            fmpz_mat_clear(A);
        }
        else
        {
            rel = _fmpz_vec_init(alloc);

            while (num_powers >= 1)
            {
                found_relation = 0;

                if (_qqbar_acb_lindep(rel, z, num_powers + 1, 1, prec) && FLINT_ABS(_fmpz_vec_max_bits(rel, num_powers + 1)) <= 10)
                {
                    if (!_fmpz_vec_is_zero(rel, num_powers + 1) && FLINT_ABS(_fmpz_vec_max_bits(rel, num_powers + 1)) <= 10)
                    {
                        if (ca_field_prove_multiplicative_relation(K, rel, z, powers, num_powers, prec, ctx))
                        {
                            ca_field_insert_multiplicative_relation(K, rel, powers, num_powers, ctx);
                        }
                    }
                }

                if (!found_relation)
                    break;

                for (j = which_removed; j < num_powers - 1; j++)
                    powers[j] = powers[j + 1];

                for (j = which_removed; j < num_powers; j++)
                    acb_swap(z + j, z + j + 1);

                num_powers--;
            }

            _fmpz_vec_clear(rel, alloc);
        }

        _acb_vec_clear(z, alloc);
    }

    flint_free(powers);
}

/* todo: move to utils */
void
fmpz_mpoly_set_coeff_si_x(fmpz_mpoly_t poly,
        slong c,
        slong x_var, slong x_exp,
        const fmpz_mpoly_ctx_t ctx);

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
        slong deg, vieta_limit;

        vieta_limit = ctx->options[CA_OPT_VIETA_LIMIT];

        for (i = 0; vieta_limit > 0 && i < len; i++)
        {
            ca_ext_struct * x = CA_FIELD_EXT_ELEM(K, i);

            /* If all conjugates of an algebraic number are present,
               add Vieta's formulas. */
            if (CA_EXT_IS_QQBAR(x) && ((deg = qqbar_degree(CA_EXT_QQBAR(x))) <= vieta_limit))
            {
                slong j, num_conjugates;

                num_conjugates = 0;

                for (j = i + 1; j < len; j++)
                {
                    ca_ext_struct * y = CA_FIELD_EXT_ELEM(K, j);

                    if (CA_EXT_IS_QQBAR(y) && fmpz_poly_equal(QQBAR_POLY(CA_EXT_QQBAR(x)), QQBAR_POLY(CA_EXT_QQBAR(y))))
                        num_conjugates++;
                    else
                        break;
                }

                if (num_conjugates + 1 == deg)
                {
                    fmpz_mpoly_t r;
                    slong k;
                    slong * vars;
                    ulong binom;
                    /* an * (x1 + x2 + ... + xn) + a[n-1] = 0, etc. */
                    /* todo: remove coefficient GCD */

                    vars = flint_malloc(deg * sizeof(slong));
                    for (k = 0; k < deg; k++)
                        vars[k] = i + k;

                    binom = 1;
                    for (k = 1; k <= deg; k++)
                    {
                        binom = binom * (deg - k + 1) / k;
                        /*  alternative way to define a cutoff:
                        if (binom > (ulong) vieta_limit)
                            break;
                        */

                        fmpz_mpoly_init(r, CA_FIELD_MCTX(K, ctx));
                        fmpz_mpoly_symmetric_gens(r, k, vars, deg, CA_FIELD_MCTX(K, ctx));

                        /* todo: cancel gcd between coefficients */
                        fmpz_mpoly_scalar_mul_fmpz(r, r, QQBAR_COEFFS(CA_EXT_QQBAR(x)) + deg, CA_FIELD_MCTX(K, ctx));
                        if (k % 2 == 1)
                            fmpz_mpoly_add_fmpz(r, r, QQBAR_COEFFS(CA_EXT_QQBAR(x)) + deg - k, CA_FIELD_MCTX(K, ctx));
                        else
                            fmpz_mpoly_sub_fmpz(r, r, QQBAR_COEFFS(CA_EXT_QQBAR(x)) + deg - k, CA_FIELD_MCTX(K, ctx));

                        _ca_field_ideal_insert_clear_mpoly(K, r, CA_FIELD_MCTX(K, ctx), ctx);
                    }

                    flint_free(vars);
                }

                i += num_conjugates;
            }
        }

        for (i = 0; i < len; i++)
        {
            slong a, b;
            ca_ext_struct * x = CA_FIELD_EXT_ELEM(K, i);

            /* Absolute annihilating polynomial for number field */
            if (CA_EXT_IS_QQBAR(x))
            {
                fmpz_mpoly_t poly;
                fmpz_mpoly_init(poly, CA_FIELD_MCTX(K, ctx));
                fmpz_mpoly_set_gen_fmpz_poly(poly, i, QQBAR_POLY(CA_EXT_QQBAR(x)), CA_FIELD_MCTX(K, ctx));
                _ca_field_ideal_insert_clear_mpoly(K, poly, CA_FIELD_MCTX(K, ctx), ctx);
                continue;
            }

            /* x = sqrt(t) -> x^2 - t = 0 */
            /* x = t^(a/b) -> x^b - t^a = 0 */
            /* A relation will only be added if t can be expressed in the
               present field. */
            if (ext_as_pow_pq(&a, &b, x, ctx))
            {
                ca_srcptr t;
                ca_field_struct * L;   /* Field of t */
                slong L_len;
                slong * tgen_map;
                slong j, k;
                int success;

                t = CA_EXT_FUNC_ARGS(x);

                if (CA_IS_SPECIAL(t))
                    flint_throw(FLINT_ERROR, "(%s)\n", __func__);

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
                    /* u^b - (p/q)^a  -->  q^a u^b - p^a */
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
                    fmpz_mpoly_pow_ui(u2, u2, b, CA_FIELD_MCTX(K, ctx));

                    if (a < 0)
                    {
                        fmpz_mpoly_swap(p, q, CA_FIELD_MCTX(K, ctx));
                        a = -a;
                    }

                    if (a != 1)
                        fmpz_mpoly_pow_ui(q, q, a, CA_FIELD_MCTX(K, ctx));

                    fmpz_mpoly_mul(u2, u2, q, CA_FIELD_MCTX(K, ctx));

                    if (a != 1)
                        fmpz_mpoly_pow_ui(p, p, a, CA_FIELD_MCTX(K, ctx));

                    fmpz_mpoly_sub(u2, u2, p, CA_FIELD_MCTX(K, ctx));

                    _ca_field_ideal_insert_clear_mpoly(K, u2, CA_FIELD_MCTX(K, ctx), ctx);
                    /* u2 does not need to be cleared */

                    fmpz_mpoly_clear(p, CA_FIELD_MCTX(K, ctx));
                    fmpz_mpoly_clear(q, CA_FIELD_MCTX(K, ctx));
                }

                flint_free(tgen_map);
            }
        }
    }

    ca_field_build_ideal_logs(K, ctx);
    ca_field_build_ideal_multiplicative(K, ctx);
    /* ca_field_build_ideal_sin_cos(K, ctx); */
    ca_field_build_ideal_erf(K, ctx);
    ca_field_build_ideal_gamma(K, ctx);

    if (ctx->options[CA_OPT_USE_GROEBNER])
    {
        slong i;
        int want_groebner;

        want_groebner = 1;
        for (i = 0; i < CA_FIELD_IDEAL_LENGTH(K); i++)
        {
            if (fmpz_mpoly_total_degree_si(CA_FIELD_IDEAL_ELEM(K, i), CA_FIELD_MCTX(K, ctx)) > ctx->options[CA_OPT_QQBAR_DEG_LIMIT])
            {
                want_groebner = 0;
                break;
            }
        }

        if (want_groebner && CA_FIELD_IDEAL(K)->length > 0)
        {
            if (ctx->options[CA_OPT_VERBOSE])
            {
                flint_printf("Attempt to compute Groebner basis for:\n    ");
                fmpz_mpoly_vec_print(CA_FIELD_IDEAL(K), CA_FIELD_MCTX(K, ctx)); flint_printf("\n\n");
            }

            /* Consider autoreduction as a preprocessing step. This is not necessarily a win,
               because the vector is usually nearly autoreduced by construction. */
            if (0)
            {
                slong before, after;
                before = CA_FIELD_IDEAL_LENGTH(K);
                fmpz_mpoly_vec_autoreduction(CA_FIELD_IDEAL(K), CA_FIELD_IDEAL(K), CA_FIELD_MCTX(K, ctx));
                after = CA_FIELD_IDEAL_LENGTH(K);
                if (before > 4)
                    flint_printf("before, after: %wd %wd\n", before, after);
            }

            if (fmpz_mpoly_buchberger_naive_with_limits(CA_FIELD_IDEAL(K), CA_FIELD_IDEAL(K),
                ctx->options[CA_OPT_GROEBNER_LENGTH_LIMIT],
                ctx->options[CA_OPT_GROEBNER_POLY_LENGTH_LIMIT],
                ctx->options[CA_OPT_GROEBNER_POLY_BITS_LIMIT],
                CA_FIELD_MCTX(K, ctx)))
            {
                fmpz_mpoly_vec_autoreduction_groebner(CA_FIELD_IDEAL(K), CA_FIELD_IDEAL(K), CA_FIELD_MCTX(K, ctx));

                if (ctx->options[CA_OPT_VERBOSE])
                {
                    flint_printf("Computed Groebner basis:\n    "); fmpz_mpoly_vec_print(CA_FIELD_IDEAL(K), CA_FIELD_MCTX(K, ctx)); flint_printf("\n\n");
                }
            }
            else
            {
                if (ctx->options[CA_OPT_VERBOSE])
                {
                    flint_printf("WARNING: Failed to compute a Groebner basis\n");
                    flint_printf("Current ideal:\n    "); fmpz_mpoly_vec_print(CA_FIELD_IDEAL(K), CA_FIELD_MCTX(K, ctx)); flint_printf("\n\n");
                }
            }

        }
    }
}
