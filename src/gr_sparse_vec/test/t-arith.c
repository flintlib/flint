/*
    Copyright (C) 2024 Kartik Venkatram and Alden Walker

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_sparse_vec.h"
#include "fmpz.h"
#include "fmpq.h"

#define CHECK_TEST(x, name) { if (GR_SUCCESS != (x)) { flint_printf("FAIL %s\n", (name)); flint_abort(); } }
int


test_add_sub_mul(flint_rand_t state, gr_ctx_t ctx)
{
    slong i, k, l, m;
    slong N = 20;
    slong n_tests = 20;
    int status = GR_SUCCESS;
    gr_sparse_vec_t vec, vec2, vec_other;
    gr_ptr dvec, dvec2, dvec3, dvec_other;
    gr_ptr temp;
    gr_ctx_t ctx_other;
    truth_t eq;

    gr_ctx_init_fmpz(ctx_other);
    GR_TMP_INIT(temp, ctx);
    GR_TMP_INIT_VEC(dvec, N, ctx);
    GR_TMP_INIT_VEC(dvec2, N, ctx);
    GR_TMP_INIT_VEC(dvec3, N, ctx);
    GR_TMP_INIT_VEC(dvec_other, N, ctx_other);
    gr_sparse_vec_init(vec, N, ctx);
    gr_sparse_vec_init(vec2, N, ctx);
    gr_sparse_vec_init(vec_other, N, ctx_other);

    for (i = 0; i < 18 * n_tests; i++)
    {
        status = GR_SUCCESS;
        //flint_printf("%d\n", i);
        k = i % 2; // Into sparse or dense
        l = (i / 2) % 3; // Add, subtract, or mul
        m = (i / 6) % 3; // For into sparse, also check other ctx usage
        if (k > 0 && m > 0) continue;

        status |= gr_sparse_vec_randtest(vec, 10, 0, state, ctx);
        if (k == 0)
        {
            status |= gr_vec_set_sparse_vec(dvec, vec, ctx);
            if (m == 0)
            {
                status |= gr_sparse_vec_randtest(vec2, 10, 0, state, ctx);
                status |= gr_vec_set_sparse_vec(dvec2, vec2, ctx);
                if (l == 0)
                {
                    status |= gr_sparse_vec_add(vec, vec, vec2, ctx);
                    status |= _gr_vec_add(dvec, dvec, dvec2, N, ctx);
                }
                else if (l == 1)
                {
                    status |= gr_sparse_vec_sub(vec, vec, vec2, ctx);
                    status |= _gr_vec_sub(dvec, dvec, dvec2, N, ctx);
                }
                else
                {
                    status |= gr_sparse_vec_mul(vec, vec, vec2, ctx);
                    status |= _gr_vec_mul(dvec, dvec, dvec2, N, ctx);
                }
            }
            else
            {
                status |= gr_sparse_vec_randtest(vec_other, 10, 0, state, ctx_other);
                status |= gr_vec_set_sparse_vec(dvec_other, vec_other, ctx_other);
                if (m == 1)
                {
                    if (l == 0)
                    {
                        status |= gr_sparse_vec_add_other(vec, vec, vec_other, ctx_other, ctx);
                        status |= _gr_vec_add_other(dvec, dvec, dvec_other, ctx_other, N, ctx);
                    }
                    else if (l == 1)
                    {
                        status |= gr_sparse_vec_sub_other(vec, vec, vec_other, ctx_other, ctx);
                        status |= _gr_vec_sub_other(dvec, dvec, dvec_other, ctx_other, N, ctx);
                    }
                    else
                    {
                        status |= gr_sparse_vec_mul_other(vec, vec, vec_other, ctx_other, ctx);
                        status |= _gr_vec_mul_other(dvec, dvec, dvec_other, ctx_other, N, ctx);
                    }
                }
                else
                {
                    if (l == 0)
                    {
                        status |= gr_other_add_sparse_vec(vec, vec_other, ctx_other, vec, ctx);
                        status |= _gr_other_add_vec(dvec, dvec_other, ctx_other, dvec, N, ctx);
                    }
                    else if (l == 1)
                    {
                        status |= gr_other_sub_sparse_vec(vec, vec_other, ctx_other, vec, ctx);
                        status |= _gr_other_sub_vec(dvec, dvec_other, ctx_other, dvec, N, ctx);
                    }
                    else
                    {
                        status |= gr_other_mul_sparse_vec(vec, vec_other, ctx_other, vec, ctx);
                        status |= _gr_other_mul_vec(dvec, dvec_other, ctx_other, dvec, N, ctx);
                    }
                }
            }
        }
        else
        {
            status |= _gr_vec_randtest(dvec, state, N, ctx);
            status |= _gr_vec_set(dvec2, dvec, N, ctx);
            status |= gr_vec_set_sparse_vec(dvec3, vec, ctx);

            if (l == 0) 
            {
                status |= gr_vec_add_sparse_vec(dvec, dvec, vec, ctx);
                status |= _gr_vec_add(dvec2, dvec2, dvec3, N, ctx);
            }
            else if (l == 1)
            {
                status |= gr_vec_sub_sparse_vec(dvec, dvec, vec, ctx);
                status |= _gr_vec_sub(dvec2, dvec2, dvec3, N, ctx);
            }
            else
            {
                continue; // TODO
            }
        }

        if (status == GR_UNABLE)
            continue;

        if (k == 0)
            status |= gr_vec_set_sparse_vec(dvec2, vec, ctx);
        eq = _gr_vec_equal(dvec, dvec2, N, ctx);
        if (eq == T_FALSE || status != GR_SUCCESS)
        {
            flint_printf(
                "\ni = %d, k = %d, l = %d, m = %d, equal = %d, status = %d\n",
                i, k, l, m, eq, status
            );
            gr_ctx_println(ctx);
            flint_printf("dvec = "); _gr_vec_print(dvec, N, ctx); flint_printf("\n");
            flint_printf("dvec2 = "); _gr_vec_print(dvec2, N, ctx); flint_printf("\n");
            return GR_TEST_FAIL;
        }
    }

    gr_sparse_vec_clear(vec, ctx);
    gr_sparse_vec_clear(vec2, ctx);
    gr_sparse_vec_clear(vec_other, ctx_other);
    GR_TMP_CLEAR_VEC(dvec, N, ctx);
    GR_TMP_CLEAR_VEC(dvec2, N, ctx);
    GR_TMP_CLEAR_VEC(dvec3, N, ctx);
    GR_TMP_CLEAR_VEC(dvec_other, N, ctx_other);
    GR_TMP_CLEAR(temp, ctx);
    gr_ctx_clear(ctx_other);
    return status;
}

#define TEST_ACCUM_MUL_SCALAR(STATUS, K, L, TYPE, VEC, VEC2, DVEC, DVEC2, DVEC3, C, CTX) { \
    if (K == 0) \
    { \
        if (L == 0) \
        { \
            STATUS |= gr_sparse_vec_addmul_##TYPE(VEC, VEC2, C, CTX); \
            STATUS |= _gr_vec_addmul_##TYPE(DVEC, DVEC2, N, C, CTX); \
        } \
        else \
        { \
            STATUS |= gr_sparse_vec_submul_##TYPE(VEC, VEC2, C, CTX); \
            STATUS |= _gr_vec_submul_##TYPE(DVEC, DVEC2, N, C, CTX); \
        } \
    } \
    else \
    { \
        if (L == 0) \
        { \
            status |= gr_vec_addmul_sparse_vec_##TYPE(DVEC, VEC, C, CTX); \
            status |= _gr_vec_addmul_##TYPE(DVEC2, DVEC3, N, C, CTX); \
        } \
        else \
        { \
            status |= gr_vec_submul_sparse_vec_##TYPE(DVEC, VEC, C, CTX); \
            status |= _gr_vec_submul_##TYPE(DVEC2, DVEC3, N, C, CTX); \
        } \
    } \
}

int test_accum_mul_scalar(flint_rand_t state, gr_ctx_t ctx)
{
    slong i, j, k, l, c;
    slong N = 20;
    slong n_tests = 20;
    int status = GR_SUCCESS;
    gr_sparse_vec_t vec, vec2, vec_other;
    gr_ptr dvec, dvec2, dvec3;
    gr_ptr temp;
    gr_ctx_t ctx_other;
    truth_t eq;

    gr_ctx_init_fmpz(ctx_other);
    GR_TMP_INIT(temp, ctx);
    GR_TMP_INIT_VEC(dvec, N, ctx);
    GR_TMP_INIT_VEC(dvec2, N, ctx);
    GR_TMP_INIT_VEC(dvec3, N, ctx);
    gr_sparse_vec_init(vec, N, ctx);
    gr_sparse_vec_init(vec2, N, ctx);
    gr_sparse_vec_init(vec_other, N, ctx_other);

    for (i = 0; i < 8 * n_tests; i++)
    {
        status = GR_SUCCESS;
        //flint_printf("%d\n", i);
        j = i % 2; // Scaler = ctx or si
        k = (i / 2) % 2; // Into sparse or dense
        l = (i / 4) % 2; // Add or subtract

        status |= gr_sparse_vec_randtest(vec, 10, 0, state, ctx);
        if (k == 0)
        {
            status |= gr_vec_set_sparse_vec(dvec, vec, ctx);
            status |= gr_sparse_vec_randtest(vec2, 10, 0, state, ctx);
            status |= gr_vec_set_sparse_vec(dvec2, vec2, ctx);
        }
        else
        {
            status |= _gr_vec_randtest(dvec, state, N, ctx);
            status |= _gr_vec_set(dvec2, dvec, N, ctx);
            status |= gr_vec_set_sparse_vec(dvec3, vec, ctx);
        }
        if (j == 0)
        {
            status |= gr_randtest_not_zero(temp, state, ctx);
            TEST_ACCUM_MUL_SCALAR(status, k, l, scalar, vec, vec2, dvec, dvec2, dvec3, temp, ctx)
            break;
        }
        else
        {
            c = n_randint(state, 0);
            TEST_ACCUM_MUL_SCALAR(status, k, l, scalar_si, vec, vec2, dvec, dvec2, dvec3, c, ctx)
            break;
        }
        if (status == GR_UNABLE)
            continue;

        if (k == 0)
            status |= gr_vec_set_sparse_vec(dvec2, vec, ctx);
        eq = _gr_vec_equal(dvec, dvec2, N, ctx);
        if (eq == T_FALSE || status != GR_SUCCESS)
        {
            flint_printf(
                "\ni = %d, j = %d, k = %d, l = %d, equal = %d, status = %d\n",
                i, j, k, l, eq, status
            );
            gr_ctx_println(ctx);
            flint_printf("dvec = "); _gr_vec_print(dvec, N, ctx); flint_printf("\n");
            flint_printf("dvec2 = "); _gr_vec_print(dvec2, N, ctx); flint_printf("\n");
            return GR_TEST_FAIL;
        }
    }

    gr_sparse_vec_clear(vec, ctx);
    gr_sparse_vec_clear(vec2, ctx);
    gr_sparse_vec_clear(vec_other, ctx_other);
    GR_TMP_CLEAR_VEC(dvec, N, ctx);
    GR_TMP_CLEAR_VEC(dvec2, N, ctx);
    GR_TMP_CLEAR_VEC(dvec3, N, ctx);
    GR_TMP_CLEAR(temp, ctx);
    gr_ctx_clear(ctx_other);
    return status;
}


#define TEST_MUL_SCALAR(STATUS, K, TYPE, VEC, VEC2, DVEC, DVEC2, C, CTX) { \
    if (K == 1) \
    { \
        STATUS |= gr_sparse_vec_div_##TYPE(VEC2, VEC, C, CTX); \
        STATUS |= _gr_vec_div_##TYPE(DVEC2, DVEC, N, C, CTX); \
    } \
    else \
    { \
        STATUS |= gr_sparse_vec_mul_##TYPE(VEC2, VEC, C, CTX); \
        STATUS |= _gr_vec_mul_##TYPE(DVEC2, DVEC, N, C, CTX); \
        if (K == 2) \
        { \
            STATUS |= gr_sparse_vec_divexact_##TYPE(VEC2, VEC2, C, CTX); \
            STATUS |= _gr_vec_divexact_##TYPE(DVEC2, DVEC2, N, C, CTX); \
        } \
    } \
}

int test_mul_div_scalar(flint_rand_t state, gr_ctx_t ctx)
{
    slong i, j, k, c;
    ulong uc;
    slong N = 100;
    slong n_tests = 20;
    int status = GR_SUCCESS;
    gr_sparse_vec_t vec, vec2;
    gr_ptr dvec, dvec2;
    gr_ptr temp;
    fmpz_t zc;
    fmpq_t qc;
    truth_t eq;

    GR_TMP_INIT(temp, ctx);
    GR_TMP_INIT_VEC(dvec, N, ctx);
    GR_TMP_INIT_VEC(dvec2, N, ctx);
    fmpz_init(zc);
    fmpq_init(qc);
    gr_sparse_vec_init(vec, N, ctx);
    gr_sparse_vec_init(vec2, N, ctx);

    for (i = 0; i < 18 * n_tests; i++)
    {
        j = i % 6; // Which type of scalar
        k = (i / 6) % 3; // Mul, div, or mul + divexact
        if ((j == 4 || k == 1) && gr_ctx_is_field(ctx) != T_TRUE)
            continue;
        if (k == 2 && (gr_ctx_is_exact(ctx) != T_TRUE || gr_ctx_is_integral_domain(ctx) != T_TRUE))
            continue;
        status |= gr_sparse_vec_randtest(vec, 10, 0, state, ctx);
        status |= gr_vec_set_sparse_vec(dvec, vec, ctx);

        switch(j)
        {
        case 0: 
            //flint_printf("Testing scalar\n");
            status |= gr_randtest_not_zero(temp, state, ctx);
            TEST_MUL_SCALAR(status, k, scalar, vec, vec2, dvec, dvec2, temp, ctx)
            break;
        case 1:
            //flint_printf("Testing scalar_si\n");
            c = n_randint(state, 0);
            TEST_MUL_SCALAR(status, k, scalar_si, vec, vec2, dvec, dvec2, c, ctx)
            break;
        case 2:
            //flint_printf("Testing scalar_ui\n");
            uc = n_randint(state, 0);
            TEST_MUL_SCALAR(status, k, scalar_ui, vec, vec2, dvec, dvec2, uc, ctx)
            break;
        case 3:
            //flint_printf("Testing scalar_fmpz\n");
            fmpz_randtest_not_zero(zc, state, 32);
            TEST_MUL_SCALAR(status, k, scalar_fmpz, vec, vec2, dvec, dvec2, zc, ctx)
            break;
        case 4:
            //flint_printf("Testing scalar_fmpq\n");
            fmpq_randtest_not_zero(qc, state, 32);
            TEST_MUL_SCALAR(status, k, scalar_fmpq, vec, vec2, dvec, dvec2, qc, ctx)
            break;
        case 5:
            //flint_printf("Testing scalar_2exp_si\n");
            c = n_randint(state, 32) + 1;
            if (k == 1)
            {
                status |= gr_sparse_vec_mul_scalar_2exp_si(vec2, vec, -c, ctx);
                status |= _gr_vec_mul_scalar_2exp_si(dvec2, dvec, N, -c, ctx);
            }
            else
            {
                status |= gr_sparse_vec_mul_scalar_2exp_si(vec2, vec, c, ctx);
                status |= _gr_vec_mul_scalar_2exp_si(dvec2, dvec, N, c, ctx);
                if (k == 2)
                {
                    status |= gr_sparse_vec_mul_scalar_2exp_si(vec2, vec2, -c, ctx);
                    status |= _gr_vec_mul_scalar_2exp_si(dvec2, dvec2, N, -c, ctx);
                }
            }
            break;
        }
        // If any operation not allowed, just skip test
        if (status == GR_UNABLE || status == GR_DOMAIN) // TODO: FIXME
        {
            status = GR_SUCCESS;
            continue;
        }
        //gr_sparse_vec_print_nz(vec, ctx); flint_printf("\n");
        //gr_sparse_vec_print_nz(vec2, ctx); flint_printf("\n");
 
        status |= gr_vec_set_sparse_vec(dvec, vec2, ctx);
        eq = _gr_vec_equal(dvec, dvec2, N, ctx);
        if (eq == T_FALSE || status != GR_SUCCESS)
        {
            flint_printf(
                "j = %d, k = %d, equal = %d, status = %d\n",
                j, k, eq, status
            );
            gr_ctx_println(ctx);
            _gr_vec_print(dvec, N, ctx); flint_printf("\n");
            _gr_vec_print(dvec2, N, ctx); flint_printf("\n");
            status |= gr_sparse_vec_set_vec(vec, dvec2, N, ctx);
            gr_sparse_vec_print_nz(vec, ctx); flint_printf("\n");
            gr_sparse_vec_print_nz(vec2, ctx); flint_printf("\n");
            return GR_TEST_FAIL;
        }
    }

    gr_sparse_vec_clear(vec, ctx);
    gr_sparse_vec_clear(vec2, ctx);
    GR_TMP_CLEAR_VEC(dvec, N, ctx);
    GR_TMP_CLEAR_VEC(dvec2, N, ctx);
    fmpz_clear(zc);
    fmpq_clear(qc);
    GR_TMP_CLEAR(temp, ctx);
    return status;
}

TEST_FUNCTION_START(gr_sparse_vec_arith, state)
{   
    slong i;
    gr_ctx_t ctx;
    for (i = 0; i < 16; ++i)
    {
        while (1)
        {
            //gr_ctx_init_fmpz(ctx);
            gr_ctx_init_random(ctx, state);
            if (gr_ctx_is_zero_ring(ctx) == T_TRUE)
                gr_ctx_clear(ctx);
            else
                break;
        }
        CHECK_TEST(test_add_sub_mul(state, ctx), "Addition, subtraction, and multiplication into sparse and dense");
        CHECK_TEST(test_accum_mul_scalar(state, ctx), "Scalar addmul and submul into sparse and dense");
        CHECK_TEST(test_mul_div_scalar(state, ctx), "Scalar multiplication and division");
        gr_ctx_clear(ctx);
    }
    TEST_FUNCTION_END(state);
}
