/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "fexpr.h"
#include "test_helpers.h"
#include "fmpq.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "fmpz_vec.h"
#include "gr_special.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"

/*
    The real and algebraic subfields (operations leaving the subfield
    fail with GR_DOMAIN, the values in the subfield are computed), the
    global generator names and gens(), and the symbolic expressions.
*/

#define CHECK_DOMAIN(expr, what) \
    do { if ((expr) != GR_DOMAIN) { flint_printf("FAIL: %s should be GR_DOMAIN\n", what); flint_abort(); } } while (0)

#define CHECK_OK(expr, what) \
    do { if ((expr) != GR_SUCCESS) { flint_printf("FAIL: %s failed\n", what); flint_abort(); } } while (0)

TEST_FUNCTION_START(gr_tower_subfields, state)
{
    gr_ctx_t QQ, R, A, C;
    gr_ptr x, y, z;
    fmpq_t q;

    gr_ctx_init_fmpq(QQ);
    gr_ctx_init_tower_lazy(R, QQ, GR_TOWER_MERGE_EXPRESS | GR_TOWER_LAZY_REAL);
    gr_ctx_init_tower_lazy(A, QQ, GR_TOWER_MERGE_EXPRESS | GR_TOWER_LAZY_ALGEBRAIC);
    gr_ctx_init_tower_lazy(C, QQ, GR_TOWER_MERGE_EXPRESS);
    fmpq_init(q);

    /* the real field */
    x = gr_heap_init(R);
    y = gr_heap_init(R);
    z = gr_heap_init(R);

    CHECK_DOMAIN(gr_i(x, R), "i");
    CHECK_OK(gr_set_si(x, -2, R), "set");
    CHECK_DOMAIN(gr_sqrt(y, x, R), "sqrt(-2)");
    CHECK_DOMAIN(gr_log(y, x, R), "log(-2)");
    fmpq_set_si(q, 1, 3);
    CHECK_DOMAIN(gr_pow_fmpq(y, x, q, R), "(-2)^(1/3)");
    CHECK_OK(gr_set_si(x, 2, R), "set");
    CHECK_DOMAIN(gr_asin(y, x, R), "asin(2)");
    CHECK_DOMAIN(gr_acos(y, x, R), "acos(2)");
    CHECK_DOMAIN(gr_atanh(y, x, R), "atanh(2)");
    CHECK_OK(gr_set_si(x, 0, R), "set");
    CHECK_DOMAIN(gr_acosh(y, x, R), "acosh(0)");

    /* values in the real field: acos(1/2) = pi/3, sqrt(2)^sqrt(2), log(2) */
    fmpq_set_si(q, 1, 2);
    CHECK_OK(gr_set_fmpq(x, q, R), "set");
    CHECK_OK(gr_acos(y, x, R), "acos(1/2)");
    CHECK_OK(gr_pi(z, R), "pi");
    CHECK_OK(gr_div_ui(z, z, 3, R), "pi/3");
    if (gr_equal(y, z, R) != T_TRUE)
    {
        flint_printf("FAIL: acos(1/2) != pi/3: "); gr_println(y, R);
        flint_abort();
    }
    CHECK_OK(gr_set_si(x, 2, R), "set");
    CHECK_OK(gr_sqrt(x, x, R), "sqrt 2");
    CHECK_OK(gr_pow(y, x, x, R), "sqrt2^sqrt2");
    CHECK_OK(gr_log(y, y, R), "log");
    CHECK_OK(gr_log(z, x, R), "log sqrt2");
    CHECK_OK(gr_mul(z, z, x, R), "sqrt2 log sqrt2");
    if (gr_equal(y, z, R) != T_TRUE)
    {
        flint_printf("FAIL: log(sqrt2^sqrt2) != sqrt2 log(sqrt2)\n");
        flint_abort();
    }
    if (gr_tower_lazy_is_real(y, R) != T_TRUE)
    {
        flint_printf("FAIL: is_real in the real field\n");
        flint_abort();
    }

    gr_heap_clear(x, R);
    gr_heap_clear(y, R);
    gr_heap_clear(z, R);

    /* the algebraic field */
    x = gr_heap_init(A);
    y = gr_heap_init(A);
    z = gr_heap_init(A);

    CHECK_DOMAIN(gr_pi(x, A), "pi");
    CHECK_OK(gr_set_si(x, 1, A), "set");
    CHECK_DOMAIN(gr_exp(y, x, A), "exp(1)");
    CHECK_DOMAIN(gr_sin(y, x, A), "sin(1)");
    CHECK_DOMAIN(gr_atan(y, x, A), "atan(1)");
    CHECK_OK(gr_set_si(x, 2, A), "set");
    CHECK_DOMAIN(gr_log(y, x, A), "log(2)");
    CHECK_OK(gr_set_si(x, 0, A), "set");
    CHECK_OK(gr_exp(y, x, A), "exp(0)");
    CHECK_OK(gr_cos(z, x, A), "cos(0)");
    if (gr_is_one(y, A) != T_TRUE || gr_is_one(z, A) != T_TRUE)
    {
        flint_printf("FAIL: exp(0), cos(0)\n");
        flint_abort();
    }
    CHECK_OK(gr_set_si(x, 2, A), "set");
    CHECK_OK(gr_sqrt(x, x, A), "sqrt 2");
    CHECK_DOMAIN(gr_pow(y, x, x, A), "sqrt2^sqrt2");
    CHECK_OK(gr_set_si(x, -2, A), "set");
    CHECK_OK(gr_sqrt(y, x, A), "sqrt(-2)");
    CHECK_OK(gr_sqr(y, y, A), "sqr");
    if (gr_equal(x, y, A) != T_TRUE)
    {
        flint_printf("FAIL: sqrt(-2)^2\n");
        flint_abort();
    }

    gr_heap_clear(x, A);
    gr_heap_clear(y, A);
    gr_heap_clear(z, A);

    /* global names and gens() in the complex field */
    {
        gr_vec_t gens;
        char * s;
        slong i;
        fexpr_t e;

        x = gr_heap_init(C);
        y = gr_heap_init(C);

        CHECK_OK(gr_set_si(x, 2, C), "set");
        CHECK_OK(gr_sqrt(x, x, C), "sqrt 2");           /* a1 */
        CHECK_OK(gr_set_si(y, 3, C), "set");
        CHECK_OK(gr_sqrt(y, y, C), "sqrt 3");           /* a2 */
        CHECK_OK(gr_exp(y, y, C), "exp(sqrt 3)");       /* t1 */
        CHECK_OK(gr_add(x, x, y, C), "add");
        CHECK_OK(gr_pi(y, C), "pi");
        CHECK_OK(gr_add(x, x, y, C), "add");

        CHECK_OK(gr_get_str(&s, x, C), "get_str");
        if (strcmp(s, "pi+t1+a1 {a1 = sqrt(2); a2 = sqrt(3); t1 = exp(a2)}") != 0)
        {
            flint_printf("FAIL: printing: %s\n", s);
            flint_abort();
        }
        flint_free(s);

        gr_vec_init(gens, 0, C);
        CHECK_OK(gr_gens(gens, C), "gens");
        if (gens->length != 4)
        {
            flint_printf("FAIL: gens: %wd\n", gens->length);
            flint_abort();
        }
        for (i = 0; i < 4; i++)
        {
            const char * names[] = { "a1 {a1 = sqrt(2)}", "a2 {a2 = sqrt(3)}", "t1 {a2 = sqrt(3); t1 = exp(a2)}", "pi" };
            CHECK_OK(gr_get_str(&s, gr_vec_entry_ptr(gens, i, C), C), "get_str");
            if (strcmp(s, names[i]) != 0)
            {
                flint_printf("FAIL: gens[%wd] = %s\n", i, s);
                flint_abort();
            }
            flint_free(s);
        }
        gr_vec_clear(gens, C);

        gr_tower_lazy_ctx_set_print(C, GR_TOWER_PRINT_NUMERIC, 6);
        CHECK_OK(gr_get_str(&s, x, C), "get_str");
        if (strcmp(s, "10.2080") != 0)
        {
            flint_printf("FAIL: numeric printing: %s\n", s);
            flint_abort();
        }
        flint_free(s);
        gr_tower_lazy_ctx_set_print(C, GR_TOWER_PRINT_SYMBOLIC | GR_TOWER_PRINT_DEFS, 6);

        fexpr_init(e);
        /* (gr_get_fexpr is defined only when fexpr.h precedes gr.h) */
        CHECK_OK(((int (*)(fexpr_t, gr_srcptr, gr_ctx_t)) C->methods[GR_METHOD_GET_FEXPR])(e, x, C), "get_fexpr");
        s = fexpr_get_str(e);
        if (strcmp(s, "Where(Add(Pi, t_1, a_1), Def(a_1, Sqrt(2)), Def(a_2, Sqrt(3)), Def(t_1, Exp(a_2)))") != 0)
        {
            flint_printf("FAIL: fexpr: %s\n", s);
            flint_abort();
        }
        flint_free(s);
        fexpr_clear(e);

        /* the printed form is read back, in this and in another context */
        {
            gr_ctx_t C2;
            gr_ptr x2, y2;
            gr_ctx_init_tower_lazy(C2, QQ, GR_TOWER_MERGE_EXPRESS);
            x2 = gr_heap_init(C2);
            CHECK_OK(gr_get_str(&s, x, C), "get_str");
            CHECK_OK(gr_set_str(y, s, C), "set_str");
            if (gr_equal(x, y, C) != T_TRUE)
            {
                flint_printf("FAIL: round trip: %s\n", s);
                flint_abort();
            }
            CHECK_OK(gr_set_str(x2, s, C2), "set_str in another context");
            CHECK_OK(gr_set_other(y, x2, C2, C), "set_other from another context");
            if (gr_equal(x, y, C) != T_TRUE)
            {
                flint_printf("FAIL: round trip through another context: %s\n", s);
                flint_abort();
            }
            flint_free(s);

            /* (elements of C2 for the computations in C2) */
            y2 = gr_heap_init(C2);
            /* a root of a polynomial, identified numerically */
            CHECK_OK(gr_set_str(x2, "a1 {a1 = root(-1 - a1 + a1^5, 0.181232 - 1.08395*i)}", C2), "root");
            CHECK_OK(gr_pow_ui(y2, x2, 5, C2), "pow");
            CHECK_OK(gr_sub(y2, y2, x2, C2), "sub");
            CHECK_OK(gr_sub_ui(y2, y2, 1, C2), "sub");
            if (gr_is_zero(y2, C2) != T_TRUE)
            {
                flint_printf("FAIL: root of x^5 - x - 1\n");
                flint_abort();
            }
            CHECK_OK(gr_set_ui(y2, 5, C2), "set");
            CHECK_OK(gr_tower_lazy_root_ui(y2, y2, 3, C2), "root_ui");
            CHECK_OK(gr_set_str(x2, "sqrt(2)", C2), "sqrt(2)");
            CHECK_OK(gr_add(y2, y2, x2, C2), "add");
            CHECK_OK(gr_set_str(x2, "a1 + a2 {a1 = root(5, 3); a2 = sqrt(2)}", C2), "defs");
            if (gr_equal(x2, y2, C2) != T_TRUE)
            {
                flint_printf("FAIL: root(5, 3)\n");
                flint_abort();
            }

            gr_heap_clear(y2, C2);
            gr_heap_clear(x2, C2);
            gr_ctx_clear(C2);
        }

        gr_heap_clear(x, C);
        gr_heap_clear(y, C);
    }

    /* views sharing the state of a complex field */
    {
        gr_ctx_t CV, RV, AV, PR;
        gr_ptr u, v;
        char * str;

        gr_ctx_init_tower_lazy(CV, QQ, GR_TOWER_MERGE_EXPRESS);
        gr_ctx_init_tower_lazy_view(RV, CV, GR_TOWER_LAZY_REAL);
        gr_ctx_init_tower_lazy_view(AV, RV, GR_TOWER_LAZY_ALGEBRAIC);
        if (!gr_tower_lazy_ctx_same_state(CV, AV) || gr_tower_lazy_ctx_field_flags(AV) != (GR_TOWER_LAZY_REAL | GR_TOWER_LAZY_ALGEBRAIC))
        {
            flint_printf("FAIL: views\n");
            flint_abort();
        }

        u = gr_heap_init(CV);
        v = gr_heap_init(CV);

        /* conversions: real elements pass, others do not */
        CHECK_OK(gr_set_ui(u, 2, CV), "set");
        CHECK_OK(gr_sqrt(u, u, CV), "sqrt");
        CHECK_OK(gr_set_other(v, u, CV, AV), "set_other to the view");
        if (gr_equal(u, v, CV) != T_TRUE)
        {
            flint_printf("FAIL: view conversion\n");
            flint_abort();
        }
        CHECK_OK(gr_pi(u, CV), "pi");
        CHECK_DOMAIN(gr_i(v, RV), "i in the real view");
        if (gr_set_other(v, u, CV, AV) == GR_SUCCESS)
        {
            flint_printf("FAIL: pi in the real algebraic view\n");
            flint_abort();
        }
        CHECK_OK(gr_i(u, CV), "i");
        CHECK_DOMAIN(gr_set_other(v, u, CV, RV), "i to the real view");

        /* parsing: auxiliary complex generators, but a real value */
        CHECK_DOMAIN(gr_set_str(v, "i", RV), "parse i");
        CHECK_DOMAIN(gr_set_str(v, "a1 {a1 = exp(2*pi*i/8)}", RV), "parse exp(2 pi i/8)");
        CHECK_OK(gr_set_str(v, "a1^2+a1^6 {a1 = exp(2*pi*i/8)}", RV), "parse 0");
        if (gr_is_zero(v, RV) != T_TRUE)
        {
            flint_printf("FAIL: parse 0\n");
            flint_abort();
        }

        /* gens(): the definitions in the field */
        {
            gr_vec_t gens;
            slong i;
            gr_vec_init(gens, 0, RV);
            CHECK_OK(gr_gens(gens, RV), "gens");
            for (i = 0; i < gens->length; i++)
            {
                if (gr_tower_lazy_is_real(gr_vec_entry_ptr(gens, i, RV), CV) != T_TRUE)
                {
                    flint_printf("FAIL: a generator of the real view is not real\n");
                    flint_abort();
                }
            }
            gr_vec_clear(gens, RV);
        }

        /* real forms: cos(pi/5) = (7 - tan(pi/5)^2)/8 (the tangent normal
           form); the real factors of x^4 + 1 are x^2 -/+ sqrt(2) x + 1
           (tan(pi/8) = sqrt(2) - 1) */
        CHECK_OK(gr_pi(u, RV), "pi");
        CHECK_OK(gr_div_ui(u, u, 5, RV), "pi/5");
        CHECK_OK(gr_cos(u, u, RV), "cos(pi/5)");
        CHECK_OK(gr_get_str(&str, u, RV), "get_str");
        if (strstr(str, "tan(pi/5)") == NULL || strstr(str, "exp(") != NULL)
        {
            flint_printf("FAIL: cos(pi/5) = %s\n", str);
            flint_abort();
        }
        flint_free(str);

        gr_ctx_init_gr_poly(PR, RV);
        {
            gr_poly_t f, c;
            gr_vec_t fac;
            fmpz_vec_t e;
            slong i;
            gr_poly_init(f, RV);
            gr_poly_init(c, RV);
            gr_vec_init(fac, 0, PR);
            fmpz_vec_init(e, 0);
            CHECK_OK(gr_poly_set_coeff_si(f, 4, 1, RV), "coeff");
            CHECK_OK(gr_poly_set_coeff_si(f, 0, 1, RV), "coeff");
            CHECK_OK(gr_factor(c, fac, e, f, 0, PR), "factor x^4 + 1");
            if (fac->length != 2)
            {
                flint_printf("FAIL: x^4 + 1 over the real view\n");
                flint_abort();
            }
            for (i = 0; i < 2; i++)
            {
                CHECK_OK(gr_get_str(&str, gr_poly_coeff_srcptr(gr_vec_entry_ptr(fac, i, PR), 1, RV), RV), "get_str");
                if (strstr(str, "sqrt(2)") == NULL || strstr(str, "exp(") != NULL)
                {
                    flint_printf("FAIL: real factor coefficient %s\n", str);
                    flint_abort();
                }
                flint_free(str);
            }
            /* x^5 - 3: the real quadratic factors in terms of tan(pi/5)
               (the real generator is not re-expressed through the fifth
               roots of unity when the towers are merged) */
            CHECK_OK(gr_poly_zero(f, RV), "zero");
            CHECK_OK(gr_poly_set_coeff_si(f, 5, 1, RV), "coeff");
            CHECK_OK(gr_poly_set_coeff_si(f, 0, -3, RV), "coeff");
            CHECK_OK(gr_factor(c, fac, e, f, 0, PR), "factor x^5 - 3");
            if (fac->length != 3)
            {
                flint_printf("FAIL: x^5 - 3 over the real view\n");
                flint_abort();
            }
            for (i = 0; i < 3; i++)
            {
                const gr_poly_struct * h = gr_vec_entry_ptr(fac, i, PR);
                if (h->length != 3)
                    continue;
                CHECK_OK(gr_get_str(&str, gr_poly_coeff_srcptr(h, 1, RV), RV), "get_str");
                if (strstr(str, "tan(pi/5)") == NULL || strstr(str, "exp(") != NULL)
                {
                    flint_printf("FAIL: real factor coefficient of x^5 - 3: %s\n", str);
                    flint_abort();
                }
                flint_free(str);
            }

            gr_vec_clear(fac, PR);
            fmpz_vec_clear(e);
            gr_poly_clear(f, RV);
            gr_poly_clear(c, RV);
        }
        gr_ctx_clear(PR);

        /* random real combinations of roots of unity (computed in the
           complex field) enter the real view in real terms, and keep
           their values */
        {
            slong it;
            gr_ptr zz, acc, back;
            zz = gr_heap_init(CV); acc = gr_heap_init(CV); back = gr_heap_init(CV);
            for (it = 0; it < 10 * flint_test_multiplier(); it++)
            {
                const slong orders[] = { 5, 8, 12, 10 };
                slong nterms = 1 + n_randint(state, 3), t2;
                CHECK_OK(gr_zero(acc, CV), "zero");
                for (t2 = 0; t2 < nterms; t2++)
                {
                    slong n = orders[n_randint(state, 4)], kk = 1 + n_randint(state, n - 1);
                    fmpq_t r;
                    fmpq_init(r);
                    /* 2 cos(2 pi k / n) = z + 1/z, z = exp(2 pi i k / n) */
                    fmpq_set_si(r, 2 * kk, n);
                    CHECK_OK(gr_pi(zz, CV), "pi");
                    CHECK_OK(gr_mul_fmpq(zz, zz, r, CV), "mul");
                    CHECK_OK(gr_i(back, CV), "i");
                    CHECK_OK(gr_mul(zz, zz, back, CV), "mul");
                    CHECK_OK(gr_exp(zz, zz, CV), "exp");
                    CHECK_OK(gr_inv(back, zz, CV), "inv");
                    CHECK_OK(gr_add(zz, zz, back, CV), "add");
                    CHECK_OK(gr_mul_si(zz, zz, (slong) n_randint(state, 7) - 3, CV), "mul");
                    CHECK_OK(gr_add(acc, acc, zz, CV), "add");
                    fmpq_clear(r);
                }
                CHECK_OK(gr_set_other(v, acc, CV, RV), "to the real view");
                CHECK_OK(gr_get_str(&str, v, RV), "get_str");
                {
                    /* (no i other than in pi, as in tan(pi/15)) */
                    char * p;
                    for (p = str; *p; p++)
                        if (p[0] == 'p' && p[1] == 'i')
                            p[0] = p[1] = 'P';
                }
                if (strstr(str, "exp(") != NULL || strstr(str, "i") != NULL)
                {
                    flint_printf("FAIL: real form of a combination of roots of unity: %s\n", str);
                    flint_abort();
                }
                flint_free(str);
                if (gr_equal(v, acc, CV) != T_TRUE)
                {
                    flint_printf("FAIL: real form changes the value\n");
                    flint_abort();
                }
            }
            gr_heap_clear(zz, CV); gr_heap_clear(acc, CV); gr_heap_clear(back, CV);
        }

        gr_heap_clear(u, CV);
        gr_heap_clear(v, CV);

        /* (the parent may be cleared before its views) */
        gr_ctx_clear(CV);
        gr_ctx_clear(AV);
        gr_ctx_clear(RV);
    }

    fmpq_clear(q);
    gr_ctx_clear(R);
    gr_ctx_clear(A);
    gr_ctx_clear(C);
    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
