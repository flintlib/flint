/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "gr_special.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"

/* The examples of examples/elementary.c (Calcium), in the lazy field. */

static gr_ctx_t K;
static gr_ptr I, PI;

#define OP(e) GR_MUST_SUCCEED(e)
static void set_qqi(gr_ptr r, slong a, slong b, slong c, slong d) { gr_ptr t = gr_heap_init(K); OP(gr_set_si(r, a, K)); OP(gr_div_si(r, r, b, K)); OP(gr_mul_si(t, I, c, K)); OP(gr_div_si(t, t, d, K)); OP(gr_add(r, r, t, K)); gr_heap_clear(t, K); }
static void sqrtui(gr_ptr r, ulong n) { OP(gr_set_ui(r, n, K)); OP(gr_sqrt(r, r, K)); }
static void logsqr_term(gr_ptr x, slong a, slong b, slong c, slong d, slong div)
{ gr_ptr y = gr_heap_init(K), z = gr_heap_init(K); set_qqi(y, a, b, c, d); OP(gr_log(y, y, K)); OP(gr_sqr(y, y, K)); OP(gr_mul(z, PI, I, K)); OP(gr_mul(y, y, z, K)); OP(gr_div_si(y, y, div, K)); OP(gr_add(x, x, y, K)); gr_heap_clear(y, K); gr_heap_clear(z, K); }
static void logpi2_term(gr_ptr x, slong a, slong b, slong c, slong d, slong div)
{ gr_ptr y = gr_heap_init(K), z = gr_heap_init(K); set_qqi(y, a, b, c, d); OP(gr_log(y, y, K)); OP(gr_sqr(z, PI, K)); OP(gr_mul(y, y, z, K)); OP(gr_div_si(y, y, div, K)); OP(gr_add(x, x, y, K)); gr_heap_clear(y, K); gr_heap_clear(z, K); }

static void check(const char * label, gr_ptr x, truth_t expected)
{
    truth_t z = gr_is_zero(x, K);
    if (z != expected)
    {
        flint_printf("FAIL: %s: got %s\n", label, z == T_TRUE ? "zero" : z == T_FALSE ? "nonzero" : "unknown");
        gr_println(x, K);
        flint_abort();
    }
}

TEST_FUNCTION_START(gr_tower_elementary, state)
{
    gr_ctx_t QQ;
    gr_ptr x, y, z;
    fmpq_t q;

    gr_ctx_init_fmpq(QQ);
    gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
    x = gr_heap_init(K); y = gr_heap_init(K); z = gr_heap_init(K); I = gr_heap_init(K); PI = gr_heap_init(K);
    fmpq_init(q);
    OP(gr_set_si(I, -1, K)); OP(gr_sqrt(I, I, K)); OP(gr_pi(PI, K));

#define T0
#define T1(label) check(label, x, T_TRUE);
#define TF(label) check(label, x, T_FALSE);

    T0 OP(gr_mul(x, PI, I, K)); OP(gr_exp(x, x, K)); OP(gr_add_ui(x, x, 1, K)); T1("Exp(Pi*I) + 1")
    T0 OP(gr_set_si(x, -1, K)); OP(gr_log(x, x, K)); OP(gr_mul(y, PI, I, K)); OP(gr_div(x, x, y, K)); OP(gr_sub_ui(x, x, 1, K)); T1("Log(-1) / (Pi*I) - 1")
    T0 OP(gr_neg(x, I, K)); OP(gr_log(x, x, K)); OP(gr_mul(y, PI, I, K)); OP(gr_div(x, x, y, K)); fmpq_set_si(q, 1, 2); OP(gr_add_fmpq(x, x, q, K)); T1("Log(-I) / (Pi*I) + 1/2")
    T0 OP(gr_set_ui(x, 10, K)); OP(gr_pow_ui(x, x, 123, K)); OP(gr_inv(x, x, K)); OP(gr_log(x, x, K)); OP(gr_set_ui(y, 100, K)); OP(gr_log(y, y, K)); OP(gr_div(x, x, y, K)); fmpq_set_si(q, 123, 2); OP(gr_add_fmpq(x, x, q, K)); T1("Log(1 / 10^123) / Log(100) + 123/2")
    T0 sqrtui(x, 2); OP(gr_add_ui(x, x, 1, K)); OP(gr_log(x, x, K)); sqrtui(y, 2); OP(gr_mul_ui(y, y, 2, K)); OP(gr_add_ui(y, y, 3, K)); OP(gr_log(y, y, K)); OP(gr_div(x, x, y, K)); fmpq_set_si(q, 1, 2); OP(gr_sub_fmpq(x, x, q, K)); T1("Log(1 + Sqrt(2)) / Log(3 + 2*Sqrt(2)) - 1/2")
    T0 sqrtui(x, 2); sqrtui(y, 3); sqrtui(z, 6); OP(gr_mul(x, x, y, K)); OP(gr_sub(x, x, z, K)); T1("Sqrt(2)*Sqrt(3) - Sqrt(6)")
    T0 sqrtui(x, 2); OP(gr_add_ui(x, x, 1, K)); OP(gr_exp(x, x, K)); sqrtui(y, 2); OP(gr_neg(y, y, K)); OP(gr_add_ui(y, y, 1, K)); OP(gr_exp(y, y, K)); OP(gr_one(z, K)); OP(gr_exp(z, z, K)); OP(gr_sqr(z, z, K)); OP(gr_mul(x, x, y, K)); OP(gr_div(x, x, z, K)); OP(gr_sub_ui(x, x, 1, K)); T1("Exp(1+Sqrt(2)) * Exp(1-Sqrt(2)) / (Exp(1)^2) - 1")
    T0 OP(gr_log(x, I, K)); OP(gr_mul(x, x, I, K)); OP(gr_exp(x, x, K)); OP(gr_div_ui(y, PI, 2, K)); OP(gr_neg(y, y, K)); OP(gr_exp(y, y, K)); OP(gr_sub(x, x, y, K)); T1("I^I - Exp(-Pi/2)")
    T0 sqrtui(x, 3); OP(gr_exp(x, x, K)); OP(gr_sqr(x, x, K)); sqrtui(y, 12); OP(gr_exp(y, y, K)); OP(gr_sub(x, x, y, K)); T1("Exp(Sqrt(3))^2 - Exp(Sqrt(12))")
    T0 OP(gr_mul(x, PI, I, K)); OP(gr_log(x, x, K)); OP(gr_mul_ui(x, x, 2, K)); OP(gr_sqrt(y, PI, K)); OP(gr_log(y, y, K)); OP(gr_mul_ui(y, y, 4, K)); OP(gr_sub(x, x, y, K)); OP(gr_mul(y, PI, I, K)); OP(gr_sub(x, x, y, K)); T1("2*Log(Pi*I) - 4*Log(Sqrt(Pi)) - Pi*I")
    T0 OP(gr_zero(x, K)); logsqr_term(x, 2,3,-2,3, -8); logsqr_term(x, 2,3,2,3, 8); logpi2_term(x, -1,1,-1,1, 12); logpi2_term(x, -1,1,1,1, 12); logpi2_term(x, 1,3,-1,3, 12); logpi2_term(x, 1,3,1,3, 12);
       OP(gr_set_ui(y, 18, K)); OP(gr_log(y, y, K)); OP(gr_sqr(z, PI, K)); OP(gr_mul(y, y, z, K)); OP(gr_div_si(y, y, -48, K)); OP(gr_sub(x, x, y, K)); T1("BBK2014 example 1")
    T0 sqrtui(x, 6); OP(gr_mul_ui(x, x, 2, K)); OP(gr_add_ui(x, x, 5, K)); OP(gr_sqrt(x, x, K)); sqrtui(y, 2); OP(gr_sub(x, x, y, K)); sqrtui(y, 3); OP(gr_sub(x, x, y, K)); T1("Sqrt(5 + 2*Sqrt(6)) - Sqrt(2) - Sqrt(3)")
    T0 OP(gr_sqrt(x, I, K)); OP(gr_add_ui(y, I, 1, K)); sqrtui(z, 2); OP(gr_div(y, y, z, K)); OP(gr_sub(x, x, y, K)); T1("Sqrt(I) - (1+I)/Sqrt(2)")
    T0 sqrtui(y, 163); OP(gr_mul(x, PI, y, K)); OP(gr_exp(x, x, K)); OP(gr_set_ui(y, 640320, K)); OP(gr_pow_ui(y, y, 3, K)); OP(gr_add_ui(y, y, 744, K)); OP(gr_sub(x, x, y, K)); TF("Exp(Pi*Sqrt(163)) - (640320^3 + 744)")


    fmpq_clear(q);
    gr_heap_clear(x, K); gr_heap_clear(y, K); gr_heap_clear(z, K); gr_heap_clear(I, K); gr_heap_clear(PI, K);
    gr_ctx_clear(K);
    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
