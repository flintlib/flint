/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Coverage of the public functions of fixed towers which the other tests
   do not exercise: predicates and inversion at a level, coordinates,
   multiplication matrices, square roots at a level, expression from an
   enclosure, refinement of a dynamic step, roots of unity and roots of
   integers, binomial irreducibility, polynomial roots, tangents and
   special functions, options tables, maps and absorption, and the flat
   representation helpers. */

#include "test_helpers.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "acb.h"
#include "qqbar.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_mat.h"
#include "gr_tower.h"

static truth_t
_equal_si(gr_srcptr x, slong c, gr_ctx_t ctx)
{
    gr_ptr t;
    truth_t res;
    GR_TMP_INIT(t, ctx);
    GR_MUST_SUCCEED(gr_set_si(t, c, ctx));
    res = gr_equal(x, t, ctx);
    GR_TMP_CLEAR(t, ctx);
    return res;
}

#define CHECK(cond, msg) do { if (!(cond)) { flint_printf("FAIL: %s\n", msg); flint_abort(); } } while (0)

static void
_api_adjoin_sqrt_ui(gr_tower_t T, ulong n, const char * name)
{
    gr_ctx_struct * top = gr_tower_field(T);
    gr_ptr x;
    GR_TMP_INIT(x, top);
    GR_MUST_SUCCEED(gr_set_ui(x, n, top));
    GR_MUST_SUCCEED(gr_tower_adjoin_sqrt(T, x, name));
    GR_TMP_CLEAR(x, top);
}

/* the generator a_k as an element of F_k */
static void
_gen_at(gr_ptr res, slong k, gr_tower_t T)
{
    GR_MUST_SUCCEED(gr_gen(res, gr_tower_field_at(T, k)));
}

/* the generator a_k as an element of the top field F_n */
static void
_gen_top(gr_ptr res, slong k, gr_tower_t T)
{
    gr_ctx_struct * Fk = gr_tower_field_at(T, k);
    gr_ptr t;
    GR_TMP_INIT(t, Fk);
    GR_MUST_SUCCEED(gr_gen(t, Fk));
    GR_MUST_SUCCEED(gr_tower_promote(res, t, k, T->length, T));
    GR_TMP_CLEAR(t, Fk);
}

TEST_FUNCTION_START(gr_tower_api, state)
{
    gr_ctx_t QQ;

    gr_ctx_init_fmpq(QQ);

    /* Q(sqrt2, sqrt3): predicates at a level, coordinates, matrices,
       square roots, expression, roots of polynomials */
    {
        gr_tower_t T;
        gr_ctx_struct * F1, * F2;
        gr_ptr a, b, x, y, z, c;
        gr_mat_t M;
        gr_poly_t p, q;
        acb_t w;
        slong D;
        int status;

        gr_tower_init(T, QQ);
        _api_adjoin_sqrt_ui(T, 2, "a");
        _api_adjoin_sqrt_ui(T, 3, "b");
        F1 = gr_tower_field_at(T, 1);
        F2 = gr_tower_field_at(T, 2);
        D = gr_tower_degree(T);
        CHECK(D == 4 && gr_tower_degree_at(T, 1) == 2, "degrees");

        GR_TMP_INIT(a, F1);
        GR_TMP_INIT5(b, x, y, z, c, F2);
        _gen_at(a, 1, T);
        _gen_at(b, 2, T);

        /* predicates and inversion at level 1 */
        {
            gr_ptr x1, y1, z1;
            GR_TMP_INIT3(x1, y1, z1, F1);
            CHECK(gr_tower_is_zero_at(a, 1, T) == T_FALSE, "is_zero_at");
            GR_MUST_SUCCEED(gr_mul(x1, a, a, F1));
            GR_MUST_SUCCEED(gr_sub_ui(x1, x1, 2, F1));
            CHECK(gr_tower_is_zero_at(x1, 1, T) == T_TRUE, "a^2 - 2 = 0 at level 1");
            GR_MUST_SUCCEED(gr_add_ui(x1, a, 1, F1));
            GR_MUST_SUCCEED(gr_tower_inv_at(y1, x1, 1, T));
            GR_MUST_SUCCEED(gr_mul(y1, y1, x1, F1));
            CHECK(gr_is_one(y1, F1) == T_TRUE, "inv_at");
            GR_MUST_SUCCEED(gr_tower_div_at(y1, a, x1, 1, T));
            GR_MUST_SUCCEED(gr_mul(y1, y1, x1, F1));
            CHECK(gr_tower_equal_at(y1, a, 1, T) == T_TRUE, "div_at");

            /* the square root of (1 + a)^2 = 3 + 2a at level 1 (k < n) */
            GR_MUST_SUCCEED(gr_mul(z1, x1, x1, F1));
            status = gr_tower_sqrt_at(y1, z1, 1, T);
            CHECK(status == GR_SUCCESS && gr_tower_equal_at(y1, x1, 1, T) == T_TRUE, "sqrt_at level 1");
            GR_MUST_SUCCEED(gr_set_ui(z1, 5, F1));
            CHECK(gr_tower_sqrt_at(y1, z1, 1, T) == GR_DOMAIN, "sqrt_at of a non-square");
            GR_TMP_CLEAR3(x1, y1, z1, F1);
        }

        /* coordinates over Q: x = 1 + 2a + 3b + 4ab at level 2 */
        GR_MUST_SUCCEED(gr_tower_promote(x, a, 1, 2, T));
        GR_MUST_SUCCEED(gr_mul(y, x, b, F2));
        GR_MUST_SUCCEED(gr_mul_ui(y, y, 4, F2));
        GR_MUST_SUCCEED(gr_mul_ui(x, x, 2, F2));
        GR_MUST_SUCCEED(gr_add(x, x, y, F2));
        GR_MUST_SUCCEED(gr_mul_ui(y, b, 3, F2));
        GR_MUST_SUCCEED(gr_add(x, x, y, F2));
        GR_MUST_SUCCEED(gr_add_ui(x, x, 1, F2));
        {
            gr_ptr coeffs;
            coeffs = gr_heap_init_vec(D, QQ);
            GR_MUST_SUCCEED(gr_tower_get_coeffs_at(coeffs, x, 2, T));
            CHECK(_equal_si(GR_ENTRY(coeffs, 0, QQ->sizeof_elem), 1, QQ) == T_TRUE, "coeffs[0]");
            CHECK(_equal_si(GR_ENTRY(coeffs, 1, QQ->sizeof_elem), 2, QQ) == T_TRUE, "coeffs[1]");
            CHECK(_equal_si(GR_ENTRY(coeffs, 2, QQ->sizeof_elem), 3, QQ) == T_TRUE, "coeffs[2]");
            CHECK(_equal_si(GR_ENTRY(coeffs, 3, QQ->sizeof_elem), 4, QQ) == T_TRUE, "coeffs[3]");
            GR_MUST_SUCCEED(gr_tower_set_coeffs_at(y, coeffs, 2, T));
            CHECK(gr_equal(x, y, F2) == T_TRUE, "set_coeffs_at round trip");
            gr_heap_clear_vec(coeffs, D, QQ);
        }

        /* the multiplication matrix of x: its characteristic polynomial is
           the characteristic polynomial of x, its determinant the norm */
        gr_mat_init(M, D, D, QQ);
        gr_poly_init(p, QQ);
        gr_poly_init(q, QQ);
        GR_MUST_SUCCEED(gr_tower_multiplication_matrix(M, x, T));
        GR_MUST_SUCCEED(gr_mat_charpoly(p, M, QQ));
        GR_MUST_SUCCEED(gr_tower_charpoly(q, x, T));
        CHECK(gr_poly_equal(p, q, QQ) == T_TRUE, "multiplication matrix charpoly");
        {
            /* the determinant is the norm: the constant term of the
               characteristic polynomial up to the sign (-1)^D */
            gr_ptr det;
            GR_TMP_INIT(det, QQ);
            GR_MUST_SUCCEED(gr_mat_det(det, M, QQ));
            CHECK(gr_equal(det, gr_poly_coeff_srcptr(q, 0, QQ), QQ) == T_TRUE, "determinant is the norm");
            GR_TMP_CLEAR(det, QQ);
        }
        gr_mat_clear(M, QQ);
        gr_poly_clear(p, QQ);
        gr_poly_clear(q, QQ);

        /* expression from an enclosure: the root -a of X^2 - 2 over F_2 */
        gr_poly_init(p, F2);
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(p, 2, 1, F2));
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(p, 0, -2, F2));
        acb_init(w);
        acb_set_d(w, -1.41421356);
        arb_add_error_2exp_si(acb_realref(w), -20);
        status = gr_tower_express_limit(y, p, w, 0, T);
        GR_MUST_SUCCEED(gr_tower_promote(x, a, 1, 2, T));
        GR_MUST_SUCCEED(gr_neg(x, x, F2));
        CHECK(status == GR_SUCCESS && gr_equal(x, y, F2) == T_TRUE, "express -sqrt(2) from an enclosure");
        status = gr_tower_express(y, p, w, T);
        CHECK(status == GR_SUCCESS && gr_equal(x, y, F2) == T_TRUE, "express (default limit)");
        gr_poly_clear(p, F2);

        /* roots of (X - a)^2 (X - b) over F_2 with multiplicities */
        {
            gr_vec_t roots;
            fmpz_vec_t mult;
            gr_poly_t lin;
            gr_ctx_t pctx;
            slong i, found_a = 0, found_b = 0;

            gr_ctx_init_gr_poly(pctx, F2);
            gr_poly_init(p, F2);
            gr_poly_init(lin, F2);
            gr_vec_init(roots, 0, F2);
            fmpz_vec_init(mult, 0);

            GR_MUST_SUCCEED(gr_poly_set_coeff_si(lin, 1, 1, F2));
            GR_MUST_SUCCEED(gr_neg(x, x, F2));   /* x = a in F_2 */
            GR_MUST_SUCCEED(gr_neg(z, x, F2));
            GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(lin, 0, z, F2));
            GR_MUST_SUCCEED(gr_poly_mul(p, lin, lin, F2));
            GR_MUST_SUCCEED(gr_neg(z, b, F2));
            GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(lin, 0, z, F2));
            GR_MUST_SUCCEED(gr_poly_mul(p, p, lin, F2));

            status = gr_tower_poly_roots(roots, mult, p, 2, T);
            CHECK(status == GR_SUCCESS && roots->length == 2, "poly_roots");
            for (i = 0; i < roots->length; i++)
            {
                gr_srcptr r = gr_vec_entry_srcptr(roots, i, F2);
                if (gr_equal(r, x, F2) == T_TRUE)
                    found_a = fmpz_get_si(mult->entries + i);
                if (gr_equal(r, b, F2) == T_TRUE)
                    found_b = fmpz_get_si(mult->entries + i);
            }
            CHECK(found_a == 2 && found_b == 1, "poly_roots multiplicities");

            gr_vec_clear(roots, F2);
            fmpz_vec_clear(mult);
            gr_poly_clear(p, F2);
            gr_poly_clear(lin, F2);
            gr_ctx_clear(pctx);
        }

        acb_clear(w);
        GR_TMP_CLEAR(a, F1);
        GR_TMP_CLEAR5(b, x, y, z, c, F2);
        gr_tower_clear(T);
    }

    /* roots of unity, roots of integers, binomial irreducibility, options */
    {
        gr_tower_t T;
        gr_ctx_struct * top;
        gr_ptr x, y;
        fmpz_t p;
        slong options[GR_TOWER_OPT_NUM_OPTIONS];
        int status;

        gr_tower_init(T, QQ);
        GR_MUST_SUCCEED(gr_tower_adjoin_root_of_unity(T, 8, "z8"));
        CHECK(gr_tower_degree(T) == 4, "degree of Q(zeta_8)");
        fmpz_init_set_ui(p, 5);
        GR_MUST_SUCCEED(gr_tower_adjoin_root_fmpz(T, p, 3, "c5"));
        CHECK(gr_tower_degree(T) == 12, "degree of Q(zeta_8, 5^(1/3))");
        top = gr_tower_field(T);
        GR_TMP_INIT2(x, y, top);

        /* zeta_8^4 = -1 and c5^3 = 5 */
        _gen_top(x, 1, T);
        GR_MUST_SUCCEED(gr_pow_ui(x, x, 4, top));
        CHECK(_equal_si(x, -1, top) == T_TRUE, "zeta_8^4");
        _gen_top(y, 2, T);
        GR_MUST_SUCCEED(gr_pow_ui(y, y, 3, top));
        CHECK(_equal_si(y, 5, top) == T_TRUE, "c5^3");

        /* a modular proof that X^3 - 7 is irreducible over the field, and
           that X^2 + 1 (= X^2 - zeta_8^4) is not (zeta_8^2 is a root) */
        {
            gr_tower_flat_struct * F = gr_tower_flat(T);
            fmpz_mpoly_q_t a;
            fmpz_mpoly_q_init(a, F->mctx);
            fmpz_mpoly_q_set_si(a, 7, F->mctx);
            status = gr_tower_binomial_irreducible_modular(T, a, F->mctx, 3, 20);
            CHECK(status == 1, "X^3 - 7 irreducible (modular)");
            fmpz_mpoly_q_set_si(a, -1, F->mctx);
            status = gr_tower_binomial_irreducible_modular(T, a, F->mctx, 2, 20);
            CHECK(status == 0, "X^2 + 1 reducible over Q(zeta_8)");
            fmpz_mpoly_q_clear(a, F->mctx);
        }

        /* an options table of the tower */
        gr_tower_options_init(options);
        options[GR_TOWER_OPT_PREC_LIMIT] = 1234;
        gr_tower_set_options(T, options);
        CHECK(GR_TOWER_OPTION(T, GR_TOWER_OPT_PREC_LIMIT) == 1234, "set_options");
        gr_tower_set_options(T, NULL);
        CHECK(GR_TOWER_OPTION(T, GR_TOWER_OPT_PREC_LIMIT) == gr_tower_option_default(GR_TOWER_OPT_PREC_LIMIT), "default options");

        fmpz_clear(p);
        GR_TMP_CLEAR2(x, y, top);
        gr_tower_clear(T);
    }

    /* refinement of a dynamic step: the root sqrt(2) of (X^2 - 2)(X^2 - 3) */
    {
        gr_tower_t T;
        gr_poly_t m, g;
        acb_t z;
        int status;

        gr_tower_init(T, QQ);
        gr_poly_init(m, QQ);
        gr_poly_init(g, QQ);
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(g, 2, 1, QQ));
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(g, 0, -2, QQ));
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, 2, 1, QQ));
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, 0, -3, QQ));
        GR_MUST_SUCCEED(gr_poly_mul(m, m, g, QQ));
        acb_init(z);
        acb_set_d(z, 1.41421356);
        arb_add_error_2exp_si(acb_realref(z), -20);
        GR_MUST_SUCCEED(gr_tower_adjoin_algebraic(T, m, z, GR_TOWER_STATUS_DYNAMIC, "r"));
        CHECK(gr_tower_step_degree(T, 1) == 4, "dynamic step of degree 4");
        status = gr_tower_refine_step(T, 1, g);
        CHECK(status == GR_SUCCESS && gr_tower_step_degree(T, 1) == 2, "refine_step to X^2 - 2");
        CHECK(gr_tower_refine(T) == 0, "nothing more to refine");
        acb_clear(z);
        gr_poly_clear(m, QQ);
        gr_poly_clear(g, QQ);
        gr_tower_clear(T);
    }

    /* tan, atan and Gamma at points of a fixed tower, checked numerically
       (adjoining a transcendental generator rebuilds the nested contexts:
       the argument is created afresh each time) */
    {
        gr_tower_t T;
        gr_ctx_struct * top;
        gr_ptr u;
        acb_t z, v;
        fmpq_t q;
        slong prec = 128;

        gr_tower_init(T, QQ);
        acb_init(z);
        acb_init(v);
        fmpq_init(q);

        top = gr_tower_field(T);
        GR_TMP_INIT(u, top);
        GR_MUST_SUCCEED(gr_set_ui(u, 1, top));
        GR_MUST_SUCCEED(gr_tower_adjoin_tan(T, u, "t"));
        GR_TMP_CLEAR(u, top);
        GR_MUST_SUCCEED(gr_tower_trans_get_acb(z, T, 1, prec));
        acb_one(v);
        acb_tan(v, v, prec);
        CHECK(acb_overlaps(z, v) && acb_rel_accuracy_bits(z) > 60, "tan(1)");

        top = gr_tower_field(T);
        GR_TMP_INIT(u, top);
        GR_MUST_SUCCEED(gr_set_ui(u, 2, top));
        GR_MUST_SUCCEED(gr_tower_adjoin_atan(T, u, "s"));
        /* the argument remains an element of the superseded top field
           (a rational function field in tan(1)), which the tower keeps */
        GR_MUST_SUCCEED(gr_add(u, u, u, top));
        CHECK(gr_is_zero(u, top) == T_FALSE, "argument usable after the rebuild");
        GR_TMP_CLEAR(u, top);
        GR_MUST_SUCCEED(gr_tower_trans_get_acb(z, T, 2, prec));
        acb_set_ui(v, 2);
        acb_atan(v, v, prec);
        CHECK(acb_overlaps(z, v) && acb_rel_accuracy_bits(z) > 60, "atan(2)");

        top = gr_tower_field(T);
        GR_TMP_INIT(u, top);
        fmpq_set_si(q, 1, 3);
        GR_MUST_SUCCEED(gr_set_fmpq(u, q, top));
        GR_MUST_SUCCEED(gr_tower_adjoin_special(T, GR_TOWER_GAMMA, 0, u, "g"));
        GR_TMP_CLEAR(u, top);
        GR_MUST_SUCCEED(gr_tower_trans_get_acb(z, T, 3, prec));
        acb_set_fmpq(v, q, prec);
        acb_gamma(v, v, prec);
        CHECK(acb_overlaps(z, v) && acb_rel_accuracy_bits(z) > 60, "Gamma(1/3)");

        CHECK(gr_tower_num_trans_si(T) == 3, "three transcendental generators");

        fmpq_clear(q);
        acb_clear(z);
        acb_clear(v);
        gr_tower_clear(T);
    }

    /* maps: the inclusion Q(sqrt2) -> Q(sqrt2, sqrt3), absorption of
       Q(sqrt3, sqrt6) into Q(sqrt2) (sqrt6 is expressed, sqrt3 adjoined),
       and the application to polynomials */
    {
        gr_tower_t U, B;
        gr_tower_map_t map;
        gr_ctx_struct * Ut, * Bt;
        gr_ptr x, y, s3, s6;
        gr_poly_t f, g;
        int status;

        gr_tower_init(U, QQ);
        gr_tower_init(B, QQ);
        _api_adjoin_sqrt_ui(U, 2, "a");
        _api_adjoin_sqrt_ui(B, 3, "b");
        _api_adjoin_sqrt_ui(B, 6, "c");

        gr_tower_map_init(map, B, U);
        status = gr_tower_absorb(U, map, B, GR_TOWER_MERGE_EXPRESS);
        CHECK(status == GR_SUCCESS, "absorb");
        CHECK(gr_tower_degree(U) == 4, "Q(sqrt2, sqrt3, sqrt6) has degree 4");

        Ut = gr_tower_field(U);
        Bt = gr_tower_field(B);
        GR_TMP_INIT2(x, y, Ut);
        GR_TMP_INIT2(s3, s6, Bt);
        _gen_top(s3, 1, B);
        _gen_top(s6, 2, B);

        GR_MUST_SUCCEED(gr_tower_map_apply(x, s3, map));
        GR_MUST_SUCCEED(gr_mul(y, x, x, Ut));
        CHECK(_equal_si(y, 3, Ut) == T_TRUE, "image of sqrt3 squares to 3");
        GR_MUST_SUCCEED(gr_tower_map_apply(y, s6, map));
        /* sqrt6 = sqrt2 sqrt3 up to sign: the image times sqrt2 is +-sqrt3's image */
        {
            gr_ptr a;
            GR_TMP_INIT(a, Ut);
            _gen_top(a, 1, U);
            GR_MUST_SUCCEED(gr_mul(y, y, a, Ut));
            GR_MUST_SUCCEED(gr_div(y, y, x, Ut));
            GR_MUST_SUCCEED(gr_mul(y, y, y, Ut));
            CHECK(_equal_si(y, 4, Ut) == T_TRUE, "sqrt6 sqrt2 / sqrt3 = +-2");
            GR_TMP_CLEAR(a, Ut);
        }

        /* a polynomial over B mapped to U: X^2 - sqrt3 maps to a polynomial
           with the image of sqrt3 as constant term */
        gr_poly_init(f, Bt);
        gr_poly_init(g, Ut);
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 2, 1, Bt));
        GR_MUST_SUCCEED(gr_neg(s3, s3, Bt));
        GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(f, 0, s3, Bt));
        GR_MUST_SUCCEED(gr_tower_map_apply_poly(g, f, map));
        GR_MUST_SUCCEED(gr_neg(x, x, Ut));
        CHECK(g->length == 3 && gr_equal(gr_poly_coeff_srcptr(g, 0, Ut), x, Ut) == T_TRUE, "map_apply_poly");
        gr_poly_clear(f, Bt);
        gr_poly_clear(g, Ut);

        GR_TMP_CLEAR2(x, y, Ut);
        GR_TMP_CLEAR2(s3, s6, Bt);
        gr_tower_map_clear(map);
        gr_tower_clear(U);
        gr_tower_clear(B);

        /* the inclusion map of a prefix */
        gr_tower_init(U, QQ);
        gr_tower_init(B, QQ);
        _api_adjoin_sqrt_ui(B, 2, "a");
        _api_adjoin_sqrt_ui(U, 2, "a");
        _api_adjoin_sqrt_ui(U, 3, "b");
        gr_tower_map_init(map, B, U);
        GR_MUST_SUCCEED(gr_tower_map_set_inclusion(map));
        Ut = gr_tower_field(U);
        Bt = gr_tower_field(B);
        GR_TMP_INIT2(x, y, Ut);
        GR_TMP_INIT(s3, Bt);
        _gen_at(s3, 1, B);
        GR_MUST_SUCCEED(gr_add_ui(s3, s3, 1, Bt));
        GR_MUST_SUCCEED(gr_tower_map_apply(x, s3, map));
        _gen_top(y, 1, U);
        GR_MUST_SUCCEED(gr_add_ui(y, y, 1, Ut));
        CHECK(gr_equal(x, y, Ut) == T_TRUE, "inclusion map");
        GR_TMP_CLEAR2(x, y, Ut);
        GR_TMP_CLEAR(s3, Bt);
        gr_tower_map_clear(map);
        gr_tower_clear(U);
        gr_tower_clear(B);
    }

    /* the flat representation: levels, conversion across contexts, and
       composition */
    {
        gr_tower_t T;
        gr_tower_flat_struct * F;
        fmpz_mpoly_ctx_struct * old;
        fmpz_mpoly_q_t x, y, z;
        gr_ptr a;
        gr_ctx_struct * F1;

        gr_tower_init(T, QQ);
        _api_adjoin_sqrt_ui(T, 2, "a");
        F = gr_tower_flat(T);
        F1 = gr_tower_field_at(T, 1);
        GR_TMP_INIT(a, F1);
        _gen_at(a, 1, T);
        GR_MUST_SUCCEED(gr_add_ui(a, a, 1, F1));

        fmpz_mpoly_q_init(x, F->mctx);
        GR_MUST_SUCCEED(gr_tower_flat_set_nested_at(x, a, 1, F));
        CHECK(gr_tower_flat_level(x, F) == 1 && gr_tower_flat_alg_level(x, F) == 1, "flat_level");
        old = F->mctx;

        /* growing the tower beyond the capacity of its context makes a
           new one; the old elements are converted */
        _api_adjoin_sqrt_ui(T, 3, "b");
        _api_adjoin_sqrt_ui(T, 5, "c");
        _api_adjoin_sqrt_ui(T, 7, "d");
        _api_adjoin_sqrt_ui(T, 11, "e");
        GR_MUST_SUCCEED(gr_tower_flat_ensure(F));
        fmpz_mpoly_q_init(y, F->mctx);
        gr_tower_flat_convert(y, x, old, F);
        CHECK(gr_tower_flat_level(y, F) == 1, "flat_convert keeps the level");

        /* composition: substitute the generator of x's context by y - 1 (so
           that 1 + a maps to y) and check the square: (1 + a)^2 = 3 + 2a */
        {
            fmpz_mpoly_q_struct ** imgs;
            slong i, nv = old->minfo->nvars;
            fmpz_mpoly_q_t t;
            imgs = flint_calloc(nv, sizeof(fmpz_mpoly_q_struct *));
            fmpz_mpoly_q_init(t, F->mctx);
            fmpz_mpoly_q_sub_si(t, y, 1, F->mctx);
            for (i = 0; i < nv; i++)
                imgs[i] = t;
            fmpz_mpoly_q_init(z, F->mctx);
            GR_MUST_SUCCEED(gr_tower_flat_compose(z, x, old, imgs, F));
            CHECK(fmpz_mpoly_q_equal(z, y, F->mctx), "flat_compose");
            fmpz_mpoly_q_mul(z, z, z, F->mctx);
            GR_MUST_SUCCEED(gr_tower_flat_reduce(z, F));
            fmpz_mpoly_q_mul_si(y, y, 2, F->mctx);
            fmpz_mpoly_q_add_si(y, y, 1, F->mctx);   /* 2(1 + a) + 1 = 3 + 2a */
            CHECK(fmpz_mpoly_q_equal(z, y, F->mctx), "flat_reduce of a square");
            fmpz_mpoly_q_clear(t, F->mctx);
            fmpz_mpoly_q_clear(z, F->mctx);
            flint_free(imgs);
        }

        fmpz_mpoly_q_clear(x, old);
        fmpz_mpoly_q_clear(y, F->mctx);
        GR_TMP_CLEAR(a, F1);
        gr_tower_clear(T);
    }

    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
