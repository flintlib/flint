/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "test_helpers.h"
#include "fmpq.h"
#include "acb.h"
#include "acb_hypgeom.h"
#include "fmpq_vec.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_special.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"

/* pFq(a; b; z) numerically for real rational parameters and argument
   (z = 1 for p = q + 1 by direct summation with a tail estimate: a
   numerical reference, not an enclosure; the parameters are chosen with
   Re(sum b - sum a) >= 2) */
static int
_hyp_reference(acb_t res, acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, int flags, slong prec)
{
    if (p == 0 && q == 1)
        acb_hypgeom_0f1(res, b, z, 0, prec);
    else if (p == 1 && q == 1)
        acb_hypgeom_m(res, a, b, z, 0, prec);
    else if (p == 2 && q == 1)
        acb_hypgeom_2f1(res, a, a + 1, b, z, flags, prec);
    else if (p == q + 1 && acb_is_one(z))
    {
        /* sum of N terms; the tail ~ t_N N / (s - 1) */
        double s = 0.0, t = 1.0, tail;
        slong n, i, N = 200000;
        double sd = 0.0;
        for (i = 0; i < q; i++)
            sd += arf_get_d(arb_midref(acb_realref(b + i)), ARF_RND_NEAR);
        for (i = 0; i < p; i++)
            sd -= arf_get_d(arb_midref(acb_realref(a + i)), ARF_RND_NEAR);
        for (n = 0; n < N; n++)
        {
            s += t;
            for (i = 0; i < p; i++)
                t *= arf_get_d(arb_midref(acb_realref(a + i)), ARF_RND_NEAR) + n;
            for (i = 0; i < q; i++)
                t /= arf_get_d(arb_midref(acb_realref(b + i)), ARF_RND_NEAR) + n;
            t /= (n + 1);
        }
        tail = t * N / (sd);
        acb_set_d(res, s + tail);
        mag_set_d(arb_radref(acb_realref(res)), 1e-8 * (1.0 + fabs(s)));
    }
    else
        acb_hypgeom_pfq(res, a, p, b, q, z, 0, prec);
    return acb_is_finite(res);
}

/* the flags of acb_hypgeom_2f1 for the rational parameters a0, a1, b
   (the integer differences, so that the reference is finite there) */
static int
_hyp_2f1_flags(const fmpq_t a0, const fmpq_t a1, const fmpq_t b)
{
    fmpq_t d;
    int flags = 0;
    fmpq_init(d);
    fmpq_sub(d, a0, a1);
    if (fmpz_is_one(fmpq_denref(d))) flags |= ACB_HYPGEOM_2F1_AB;
    fmpq_sub(d, a0, b);
    if (fmpz_is_one(fmpq_denref(d))) flags |= ACB_HYPGEOM_2F1_AC;
    fmpq_sub(d, a1, b);
    if (fmpz_is_one(fmpq_denref(d))) flags |= ACB_HYPGEOM_2F1_BC;
    fmpq_add(d, a0, a1);
    fmpq_sub(d, d, b);
    if (fmpz_is_one(fmpq_denref(d))) flags |= ACB_HYPGEOM_2F1_ABC;
    fmpq_clear(d);
    return flags;
}

static void
_rand_fmpq(fmpq_t x, flint_rand_t state, slong num_bound, slong den_bound)
{
    fmpz_set_si(fmpq_numref(x), (slong) n_randint(state, 2 * num_bound + 1) - num_bound);
    fmpz_set_si(fmpq_denref(x), 1 + n_randint(state, den_bound));
    fmpq_canonicalise(x);
}

/* whether x (an element of K) evaluates to an enclosure overlapping y */
static int
_hyp_check_overlap(gr_srcptr x, const acb_t y, gr_ctx_t K)
{
    acb_t z;
    int ok;
    acb_init(z);
    ok = (gr_tower_lazy_get_acb(z, x, 128, K) == GR_SUCCESS) && acb_overlaps(z, y);
    acb_clear(z);
    return ok;
}

/*
    Table entry i at the symbol values sym (those used) and z = c, against
    the numerical value of its left-hand side: returns 1 if checked, 0 if
    the entry does not apply (or the reference is not available there);
    aborts on a mismatch.
*/
static int
_hyp_table_check(slong i, slong p, slong q, const int * used, const fmpq * sym, const fmpq_t c, gr_ctx_t K)
{
    gr_srcptr vals[26];
    gr_ptr symv, params, value;
    acb_ptr av;
    fmpq * pq;
    acb_t zv, ref;
    int status, cond, ok = 1, checked = 0;
    slong k;

    symv = gr_heap_init_vec(26, K);
    params = gr_heap_init_vec(FLINT_MAX(p + q, 1), K);
    value = gr_heap_init(K);
    av = _acb_vec_init(FLINT_MAX(p + q, 1));
    pq = _fmpq_vec_init(FLINT_MAX(p + q, 3));
    acb_init(zv);
    acb_init(ref);

    for (k = 0; k < 26; k++)
        vals[k] = NULL;
    for (k = 0; k < 25; k++)
    {
        if (!used[k])
            continue;
        GR_MUST_SUCCEED(gr_set_fmpq(GR_ENTRY(symv, k, K->sizeof_elem), sym + k, K));
        vals[k] = GR_ENTRY(symv, k, K->sizeof_elem);
    }
    GR_MUST_SUCCEED(gr_set_fmpq(GR_ENTRY(symv, 25, K->sizeof_elem), c, K));
    vals[25] = GR_ENTRY(symv, 25, K->sizeof_elem);
    arb_set_fmpq(acb_realref(zv), c, 256);

    status = _gr_tower_hypgeom_table_eval(params, &cond, value, i, vals, K);
    if (status != GR_SUCCESS || cond != 1)
        ok = 0;

    /* rational parameters, not nonpositive integers; for the sums
       at z = 1, sum b - sum a >= 2 */
    for (k = 0; k < p + q && ok; k++)
    {
        if (gr_get_fmpq(pq + k, GR_ENTRY(params, k, K->sizeof_elem), K) != GR_SUCCESS ||
            (fmpz_is_one(fmpq_denref(pq + k)) && fmpz_sgn(fmpq_numref(pq + k)) <= 0))
            ok = 0;
        else
            arb_set_fmpq(acb_realref(av + k), pq + k, 256);
    }
    if (ok && p == q + 1 && fmpq_is_one(c))
    {
        double sd = 0.0;
        for (k = 0; k < p + q; k++)
            sd += ((k < p) ? -1 : 1) * arf_get_d(arb_midref(acb_realref(av + k)), ARF_RND_NEAR);
        if (sd < 2.0)
            ok = 0;
    }

    if (ok && _hyp_reference(ref, av, p, av + p, q, zv,
            (p == 2 && q == 1) ? _hyp_2f1_flags(pq + 0, pq + 1, pq + 2) : 0, 256))
    {
        if (!_hyp_check_overlap(value, ref, K))
        {
            flint_printf("FAIL: table entry %wd: %s\n", i, _gr_tower_hypgeom_table_entry(i));
            flint_printf("params = "); _gr_vec_print(params, p + q, K); flint_printf("\n");
            flint_printf("z = "); fmpq_print(c); flint_printf("\n");
            flint_printf("value = "); gr_println(value, K);
            flint_printf("reference = "); acb_printn(ref, 30, 0); flint_printf("\n");
            flint_abort();
        }
        checked = 1;
    }

    acb_clear(zv);
    acb_clear(ref);
    _acb_vec_clear(av, FLINT_MAX(p + q, 1));
    _fmpq_vec_clear(pq, FLINT_MAX(p + q, 3));
    gr_heap_clear(value, K);
    gr_heap_clear_vec(params, FLINT_MAX(p + q, 1), K);
    gr_heap_clear_vec(symv, 26, K);
    return checked;
}

TEST_FUNCTION_START(gr_tower_hypgeom, state)
{
    gr_ctx_t QQ, K;
    slong iter, i, n;

    gr_ctx_init_fmpq(QQ);
    gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);

    /* every entry of the table, at random symbols, against a numerical
       evaluation of its left-hand side */
    n = _gr_tower_hypgeom_table_length();
    for (i = 0; i < n; i++)
    {
        slong p, q, k, tries, checked = 0;
        int used[25], zfree;
        fmpq_t z0;

        fmpq_init(z0);
        if (_gr_tower_hypgeom_table_shape(&p, &q, used, &zfree, z0, i) != GR_SUCCESS)
        {
            flint_printf("FAIL: invalid table entry %wd: %s\n", i, _gr_tower_hypgeom_table_entry(i));
            flint_abort();
        }

        for (tries = 0; tries < 40 && checked < 2; tries++)
        {
            fmpq sym[25];
            fmpq_t c;

            fmpq_init(c);
            for (k = 0; k < 25; k++)
            {
                fmpq_init(sym + k);
                if (used[k])
                    _rand_fmpq(sym + k, state, 7, 4);
            }
            if (zfree)
            {
                do { _rand_fmpq(c, state, 4, 9); } while (fmpq_is_zero(c));
                fmpq_div_fmpz(c, c, (fmpz[1]) {2});
            }
            else
                fmpq_set(c, z0);

            checked += _hyp_table_check(i, p, q, used, sym, c, K);

            for (k = 0; k < 25; k++)
                fmpq_clear(sym + k);
            fmpq_clear(c);
        }

        /* symbols at half-integers, quarters and negative values, where
           parameters of the functions on the right become nonpositive
           integers (2F1(1/4, 3/4; -1/2; 3/8) by entry "a, a+1/2 | c":
           2F1(1/2, -1; -2; w) on the right, whose value there is the
           limit in c, not the terminating series) */
        for (tries = 0; tries < 12; tries++)
        {
            static const slong vn[] = { -5, -3, -1, -1, 1, 1, 3, 3, 1, -3 };
            static const slong vd[] = { 2, 2, 2, 4, 4, 2, 4, 2, 3, 4 };
            static const slong zn[] = { 3, -1, 1, -3 };
            static const slong zd[] = { 8, 3, 5, 4 };
            fmpq sym[25];
            fmpq_t c;
            slong j;

            fmpq_init(c);
            for (k = 0; k < 25; k++)
            {
                fmpq_init(sym + k);
                if (used[k])
                {
                    j = n_randint(state, 10);
                    fmpq_set_si(sym + k, vn[j], vd[j]);
                }
            }
            if (zfree)
            {
                j = n_randint(state, 4);
                fmpq_set_si(c, zn[j], zd[j]);
            }
            else
                fmpq_set(c, z0);

            (void) _hyp_table_check(i, p, q, used, sym, c, K);

            for (k = 0; k < 25; k++)
                fmpq_clear(sym + k);
            fmpq_clear(c);
        }

        if (checked == 0)
        {
            flint_printf("FAIL: table entry %wd not checked: %s\n", i, _gr_tower_hypgeom_table_entry(i));
            flint_abort();
        }

        fmpq_clear(z0);
    }

    /* random 0F1, 1F1, 2F1 with rational parameters against acb (the
       contiguity module, the cells, the transformations, the reductions) */
    for (iter = 0; iter < 60 * flint_test_multiplier(); iter++)
    {
        slong p, q, k;
        fmpq_t c;
        gr_ptr a, b, z, y;
        acb_ptr av;
        acb_t zv, ref;
        fmpq pq[3];
        int status, flags;

        p = n_randint(state, 3);
        q = 1;
        fmpq_init(c);
        a = gr_heap_init_vec(3, K);
        b = gr_heap_init(K);
        z = gr_heap_init(K);
        y = gr_heap_init(K);
        av = _acb_vec_init(4);
        for (k = 0; k < 3; k++)
            fmpq_init(pq + k);
        acb_init(zv);
        acb_init(ref);

        for (k = 0; k < p + q; k++)
        {
            /* small denominators: closed forms and generators */
            fmpz_set_si(fmpq_numref(c), (slong) n_randint(state, 13) - 6);
            fmpz_set_si(fmpq_denref(c), 1 + n_randint(state, 4));
            fmpq_canonicalise(c);
            if (k == p + q - 1 && fmpz_is_one(fmpq_denref(c)) && fmpz_sgn(fmpq_numref(c)) <= 0)
                fmpq_set_si(c, 1, 2);
            arb_set_fmpq(acb_realref(av + k), c, 256);
            fmpq_set(pq + k, c);
            GR_MUST_SUCCEED(gr_set_fmpq((k < p) ? GR_ENTRY(a, k, K->sizeof_elem) : b, c, K));
        }
        flags = (p == 2) ? _hyp_2f1_flags(pq + 0, pq + 1, pq + 2) : 0;

        do { _rand_fmpq(c, state, 6, 7); } while (fmpq_is_zero(c));
        if (p == 2)
            fmpq_div_fmpz(c, c, (fmpz[1]) {3});   /* (mostly inside the unit disk, some outside) */
        arb_set_fmpq(acb_realref(zv), c, 256);
        GR_MUST_SUCCEED(gr_set_fmpq(z, c, K));

        /* a complex argument a third of the time (the branches of the
           factors (1 - z)^a, z^a of the transformations, the cells) */
        if (n_randint(state, 3) == 0)
        {
            gr_ptr t;
            GR_TMP_INIT(t, K);
            do { _rand_fmpq(c, state, 4, 5); } while (fmpq_is_zero(c));
            if (p == 2)
                fmpq_div_fmpz(c, c, (fmpz[1]) {3});
            arb_set_fmpq(acb_imagref(zv), c, 256);
            GR_MUST_SUCCEED(gr_i(t, K));
            GR_MUST_SUCCEED(gr_mul_fmpq(t, t, c, K));
            GR_MUST_SUCCEED(gr_add(z, z, t, K));
            GR_TMP_CLEAR(t, K);
            fmpq_zero(c);   /* (not on the cut) */
        }

        if (p == 0)
            status = gr_hypgeom_0f1(y, b, z, 0, K);
        else if (p == 1)
            status = gr_hypgeom_1f1(y, a, b, z, 0, K);
        else
            status = gr_hypgeom_2f1(y, a, GR_ENTRY(a, 1, K->sizeof_elem), b, z, 0, K);

        _hyp_reference(ref, av, p, av + p, q, zv, flags, 256);

        if (status == GR_SUCCESS)
        {
            if (!_hyp_check_overlap(y, ref, K))
            {
                flint_printf("FAIL: value (p = %wd)\n", p);
                flint_printf("a = "); _gr_vec_print(a, p, K); flint_printf("\n");
                flint_printf("b = "); gr_println(b, K);
                flint_printf("z = "); gr_println(z, K);
                flint_printf("y = "); gr_println(y, K);
                flint_printf("ref = "); acb_printn(ref, 30, 0); flint_printf("\n");
                flint_abort();
            }
        }
        else if (acb_is_finite(ref) && !(p == 2 && fmpq_cmp_ui(c, 1) >= 0))
        {
            /* (the cut z >= 1 of 2F1 and poles aside, every value is
               expected) */
            flint_printf("FAIL: status %d (p = %wd)\n", status, p);
            flint_printf("a = "); _gr_vec_print(a, p, K); flint_printf("\n");
            flint_printf("b = "); gr_println(b, K);
            flint_printf("z = "); gr_println(z, K);
            flint_abort();
        }

        fmpq_clear(c);
        gr_heap_clear_vec(a, 3, K);
        gr_heap_clear(b, K);
        gr_heap_clear(z, K);
        gr_heap_clear(y, K);
        _acb_vec_clear(av, 4);
        for (k = 0; k < 3; k++)
            fmpq_clear(pq + k);
        acb_clear(zv);
        acb_clear(ref);
    }

    /* identities decided exactly */
    {
        gr_ptr a, b, c, z, u, v, w, t, s;
        GR_TMP_INIT5(a, b, c, z, u, K);
        GR_TMP_INIT4(v, w, t, s, K);

        /* Gauss's contiguous relation in c, with irrational parameters:
           c (c - 1) (z - 1) F(c - 1) + c (c - 1 - (2c - a - b - 1) z) F(c)
             + (c - a)(c - b) z F(c + 1) = 0 */
        GR_MUST_SUCCEED(gr_set_str(a, "sqrt(2)", K));
        GR_MUST_SUCCEED(gr_set_str(b, "1/3", K));
        GR_MUST_SUCCEED(gr_set_str(c, "7/5", K));
        GR_MUST_SUCCEED(gr_set_str(z, "1/4", K));

        GR_MUST_SUCCEED(gr_sub_ui(t, c, 1, K));
        GR_MUST_SUCCEED(gr_hypgeom_2f1(u, a, b, t, z, 0, K));
        GR_MUST_SUCCEED(gr_mul(u, u, c, K));
        GR_MUST_SUCCEED(gr_mul(u, u, t, K));
        GR_MUST_SUCCEED(gr_sub_ui(s, z, 1, K));
        GR_MUST_SUCCEED(gr_mul(u, u, s, K));

        GR_MUST_SUCCEED(gr_hypgeom_2f1(v, a, b, c, z, 0, K));
        GR_MUST_SUCCEED(gr_mul_ui(s, c, 2, K));
        GR_MUST_SUCCEED(gr_sub(s, s, a, K));
        GR_MUST_SUCCEED(gr_sub(s, s, b, K));
        GR_MUST_SUCCEED(gr_sub_ui(s, s, 1, K));
        GR_MUST_SUCCEED(gr_mul(s, s, z, K));
        GR_MUST_SUCCEED(gr_sub(s, t, s, K));
        GR_MUST_SUCCEED(gr_mul(s, s, c, K));
        GR_MUST_SUCCEED(gr_mul(v, v, s, K));

        GR_MUST_SUCCEED(gr_add_ui(t, c, 1, K));
        GR_MUST_SUCCEED(gr_hypgeom_2f1(w, a, b, t, z, 0, K));
        GR_MUST_SUCCEED(gr_sub(s, c, a, K));
        GR_MUST_SUCCEED(gr_mul(w, w, s, K));
        GR_MUST_SUCCEED(gr_sub(s, c, b, K));
        GR_MUST_SUCCEED(gr_mul(w, w, s, K));
        GR_MUST_SUCCEED(gr_mul(w, w, z, K));

        GR_MUST_SUCCEED(gr_add(u, u, v, K));
        GR_MUST_SUCCEED(gr_add(u, u, w, K));
        if (gr_is_zero(u, K) != T_TRUE)
        {
            flint_printf("FAIL: Gauss's contiguous relation\n");
            gr_println(u, K);
            flint_abort();
        }

        /* Euler's and Pfaff's transformations (canonical forms):
           F(a, b; c; z) = (1 - z)^(c - a - b) F(c - a, c - b; c; z)
                         = (1 - z)^(-a) F(a, c - b; c; z / (z - 1)) */
        GR_MUST_SUCCEED(gr_set_str(a, "1/3", K));
        GR_MUST_SUCCEED(gr_set_str(b, "2/7", K));
        GR_MUST_SUCCEED(gr_set_str(c, "4/5", K));
        GR_MUST_SUCCEED(gr_set_str(z, "-2/3", K));
        GR_MUST_SUCCEED(gr_hypgeom_2f1(u, a, b, c, z, 0, K));
        GR_MUST_SUCCEED(gr_sub(s, c, a, K));
        GR_MUST_SUCCEED(gr_sub(t, c, b, K));
        GR_MUST_SUCCEED(gr_hypgeom_2f1(v, s, t, c, z, 0, K));
        GR_MUST_SUCCEED(gr_sub(w, c, a, K));
        GR_MUST_SUCCEED(gr_sub(w, w, b, K));
        GR_MUST_SUCCEED(gr_sub_ui(s, z, 1, K));
        GR_MUST_SUCCEED(gr_neg(s, s, K));
        GR_MUST_SUCCEED(gr_pow(s, s, w, K));
        GR_MUST_SUCCEED(gr_mul(v, v, s, K));
        if (gr_equal(u, v, K) != T_TRUE)
        {
            flint_printf("FAIL: Euler's transformation\n");
            flint_abort();
        }
        GR_MUST_SUCCEED(gr_sub_ui(s, z, 1, K));
        GR_MUST_SUCCEED(gr_div(w, z, s, K));
        GR_MUST_SUCCEED(gr_hypgeom_2f1(v, a, t, c, w, 0, K));
        GR_MUST_SUCCEED(gr_neg(s, s, K));
        GR_MUST_SUCCEED(gr_neg(w, a, K));
        GR_MUST_SUCCEED(gr_pow(s, s, w, K));
        GR_MUST_SUCCEED(gr_mul(v, v, s, K));
        if (gr_equal(u, v, K) != T_TRUE)
        {
            flint_printf("FAIL: Pfaff's transformation\n");
            flint_abort();
        }

        /* Kummer's transformation 1F1(a; b; z) = exp(z) 1F1(b - a; b; -z) */
        GR_MUST_SUCCEED(gr_set_str(a, "1/3+sqrt(3)", K));
        GR_MUST_SUCCEED(gr_set_str(b, "5/7", K));
        GR_MUST_SUCCEED(gr_set_str(z, "2/5", K));
        GR_MUST_SUCCEED(gr_hypgeom_1f1(u, a, b, z, 0, K));
        GR_MUST_SUCCEED(gr_sub(s, b, a, K));
        GR_MUST_SUCCEED(gr_neg(t, z, K));
        GR_MUST_SUCCEED(gr_hypgeom_1f1(v, s, b, t, 0, K));
        GR_MUST_SUCCEED(gr_exp(s, z, K));
        GR_MUST_SUCCEED(gr_mul(v, v, s, K));
        if (gr_equal(u, v, K) != T_TRUE)
        {
            flint_printf("FAIL: Kummer's transformation\n");
            flint_abort();
        }

        /* symmetry in the parameters (irrational a: generators) */
        GR_MUST_SUCCEED(gr_set_str(a, "sqrt(2)", K));
        GR_MUST_SUCCEED(gr_set_str(b, "1/3", K));
        GR_MUST_SUCCEED(gr_set_str(c, "7/5", K));
        GR_MUST_SUCCEED(gr_set_str(z, "1/4 + i/5", K));
        GR_MUST_SUCCEED(gr_hypgeom_2f1(u, a, b, c, z, 0, K));
        GR_MUST_SUCCEED(gr_hypgeom_2f1(v, b, a, c, z, 0, K));
        if (gr_equal(u, v, K) != T_TRUE)
        {
            flint_printf("FAIL: 2F1(a, b; c; z) = 2F1(b, a; c; z)\n");
            flint_abort();
        }

        /* integer parameter differences at inexact parameters (the
           numerical values of the generators by the limits) */
        GR_MUST_SUCCEED(gr_set_str(a, "2/3", K));
        GR_MUST_SUCCEED(gr_set_str(b, "1/3", K));
        GR_MUST_SUCCEED(gr_set_str(c, "1", K));
        GR_MUST_SUCCEED(gr_set_str(z, "9/10", K));
        GR_MUST_SUCCEED(gr_hypgeom_2f1(u, a, b, c, z, 0, K));
        GR_MUST_SUCCEED(gr_set_str(z, "1 + i/2", K));
        GR_MUST_SUCCEED(gr_hypgeom_2f1(u, a, b, c, z, 0, K));
        {
            acb_t x, y, zz, cc, r;
            acb_init(x); acb_init(y); acb_init(zz); acb_init(cc); acb_init(r);
            arb_set_si(acb_realref(x), 2); arb_div_ui(acb_realref(x), acb_realref(x), 3, 256);
            arb_set_si(acb_realref(y), 1); arb_div_ui(acb_realref(y), acb_realref(y), 3, 256);
            acb_one(cc);
            acb_set_d_d(zz, 1.0, 0.5);
            acb_hypgeom_2f1(r, x, y, cc, zz, ACB_HYPGEOM_2F1_AB | ACB_HYPGEOM_2F1_ABC, 256);
            if (!_hyp_check_overlap(u, r, K))
            {
                flint_printf("FAIL: 2F1(2/3, 1/3; 1; 1 + i/2)\n");
                flint_abort();
            }
            acb_clear(x); acb_clear(y); acb_clear(zz); acb_clear(cc); acb_clear(r);
        }

        /* the restrictions of the views: 2F1(1, 1; 2; 2) = -log(-1)/2 is
           not real; 2F1(1, 1; 2; 1/2) = 2 log 2 not known to be algebraic */
        {
            gr_ctx_t KR, KA;
            gr_ptr x1, x2, w2, r;
            gr_ctx_init_tower_lazy(KR, QQ, GR_TOWER_MERGE_EXPRESS | GR_TOWER_LAZY_REAL);
            gr_ctx_init_tower_lazy(KA, QQ, GR_TOWER_MERGE_EXPRESS | GR_TOWER_LAZY_ALGEBRAIC);
            GR_TMP_INIT4(x1, x2, w2, r, KR);
            GR_MUST_SUCCEED(gr_one(x1, KR));
            GR_MUST_SUCCEED(gr_set_ui(x2, 2, KR));
            if (gr_hypgeom_2f1(r, x1, x1, x2, x2, 0, KR) != GR_DOMAIN)
            {
                flint_printf("FAIL: 2F1(1, 1; 2; 2) in the real view\n");
                flint_abort();
            }
            GR_TMP_CLEAR4(x1, x2, w2, r, KR);
            GR_TMP_INIT4(x1, x2, w2, r, KA);
            GR_MUST_SUCCEED(gr_one(x1, KA));
            GR_MUST_SUCCEED(gr_set_ui(x2, 2, KA));
            GR_MUST_SUCCEED(gr_set_str(w2, "1/2", KA));
            if (gr_hypgeom_2f1(r, x1, x1, x2, w2, 0, KA) == GR_SUCCESS)
            {
                flint_printf("FAIL: 2F1(1, 1; 2; 1/2) in the algebraic view\n");
                flint_abort();
            }
            GR_TMP_CLEAR4(x1, x2, w2, r, KA);
            gr_ctx_clear(KR);
            gr_ctx_clear(KA);
        }

        /* regularization */
        GR_MUST_SUCCEED(gr_set_str(a, "2/3", K));
        GR_MUST_SUCCEED(gr_set_str(b, "1/5", K));
        GR_MUST_SUCCEED(gr_set_str(c, "-2", K));
        GR_MUST_SUCCEED(gr_set_str(z, "1/7", K));
        if (gr_hypgeom_2f1(u, a, b, c, z, 0, K) != GR_DOMAIN)
        {
            flint_printf("FAIL: pole at c = -2\n");
            flint_abort();
        }
        GR_MUST_SUCCEED(gr_hypgeom_2f1(u, a, b, c, z, 1, K));
        {
            acb_t x, y, zz, cc;
            acb_init(x); acb_init(y); acb_init(zz); acb_init(cc);
            acb_set_d(x, 2.0 / 3);
            acb_set_d(y, 1.0 / 5);
            arb_set_si(acb_realref(cc), -2);
            arb_set_d(acb_realref(zz), 1.0 / 7);
            arb_set_si(acb_realref(x), 2); arb_div_ui(acb_realref(x), acb_realref(x), 3, 128);
            arb_set_si(acb_realref(y), 1); arb_div_ui(acb_realref(y), acb_realref(y), 5, 128);
            arb_set_si(acb_realref(zz), 1); arb_div_ui(acb_realref(zz), acb_realref(zz), 7, 128);
            acb_hypgeom_2f1(x, x, y, cc, zz, 1, 128);
            if (!_hyp_check_overlap(u, x, K))
            {
                flint_printf("FAIL: regularized 2F1 at c = -2\n");
                flint_abort();
            }
            acb_clear(x); acb_clear(y); acb_clear(zz); acb_clear(cc);
        }

        /* elliptic integrals through 2F1 */
        GR_MUST_SUCCEED(gr_set_str(z, "1/2", K));
        GR_MUST_SUCCEED(gr_elliptic_e(u, z, K));
        GR_MUST_SUCCEED(gr_set_str(v, "gamma(1/4)^2/(8*sqrt(pi)) + pi^(3/2)/gamma(1/4)^2", K));
        if (gr_equal(u, v, K) != T_TRUE)
        {
            flint_printf("FAIL: E(1/2)\n");
            gr_println(u, K);
            flint_abort();
        }

        GR_TMP_CLEAR5(a, b, c, z, u, K);
        GR_TMP_CLEAR4(v, w, t, s, K);
    }

    /* the links to the generators of the context: quadratic
       transformations of 2F1 (both orders of evaluation), Landen's
       transformation of K and E */
    {
        static const char * qt[][9] = {
            /* a, b, c, z; A, B, C, u; F(a, b; c; z) / F(A, B; C; u) */
            {"1/3", "1/5", "2/5", "1/7", "1/6", "2/3", "7/10", "1/169", "(13/14)^(-1/3)"},
            {"1/3", "1/5", "17/15", "1/7", "1/6", "2/3", "17/15", "7/16", "(8/7)^(-1/3)"},
            {"1/3", "1/5", "23/30", "1/7", "1/6", "1/10", "23/30", "24/49", "1"},
            {"1/3", "1/5", "17/15", "-2/3", "1/6", "2/3", "17/15", "-24", "3^(1/3)"},
        };
        slong k, order;

        for (k = 0; k < 4; k++)
        {
            for (order = 0; order < 2; order++)
            {
                gr_ctx_t K2;
                gr_ptr p[9], x, y;
                slong i;

                gr_ctx_init_tower_lazy(K2, QQ, GR_TOWER_MERGE_EXPRESS);
                for (i = 0; i < 9; i++)
                {
                    p[i] = gr_heap_init(K2);
                    GR_MUST_SUCCEED(gr_set_str(p[i], qt[k][i], K2));
                }
                x = gr_heap_init(K2);
                y = gr_heap_init(K2);
                if (order == 0)
                {
                    GR_MUST_SUCCEED(gr_hypgeom_2f1(y, p[4], p[5], p[6], p[7], 0, K2));
                    GR_MUST_SUCCEED(gr_hypgeom_2f1(x, p[0], p[1], p[2], p[3], 0, K2));
                }
                else
                {
                    GR_MUST_SUCCEED(gr_hypgeom_2f1(x, p[0], p[1], p[2], p[3], 0, K2));
                    GR_MUST_SUCCEED(gr_hypgeom_2f1(y, p[4], p[5], p[6], p[7], 0, K2));
                }
                GR_MUST_SUCCEED(gr_mul(y, y, p[8], K2));
                if (gr_equal(x, y, K2) != T_TRUE)
                {
                    flint_printf("FAIL: quadratic transformation %wd (order %wd)\n", k, order);
                    gr_println(x, K2);
                    gr_println(y, K2);
                    flint_abort();
                }
                for (i = 0; i < 9; i++)
                    gr_heap_clear(p[i], K2);
                gr_heap_clear(x, K2);
                gr_heap_clear(y, K2);
                gr_ctx_clear(K2);
            }
        }

        /* the (1/4, 3/4; 1) and (1/4, 1/4; 1) families at the images of
           K(1/7), and Landen's transformation (one and two steps, both
           orders) */
        for (order = 0; order < 2; order++)
        {
            gr_ctx_t K2;
            gr_ptr m, m1, m2, k1, k2, x, y, w;

            gr_ctx_init_tower_lazy(K2, QQ, GR_TOWER_MERGE_EXPRESS);
            GR_TMP_INIT4(m, m1, m2, k1, K2);
            GR_TMP_INIT4(k2, x, y, w, K2);

            if (order == 0)
            {
                GR_MUST_SUCCEED(gr_set_str(m, "1/7", K2));
                GR_MUST_SUCCEED(gr_elliptic_k(w, m, K2));
                GR_MUST_SUCCEED(gr_set_str(x, "1/4", K2));
                GR_MUST_SUCCEED(gr_set_str(y, "3/4", K2));
                GR_MUST_SUCCEED(gr_set_str(k1, "7/16", K2));
                GR_MUST_SUCCEED(gr_one(k2, K2));
                GR_MUST_SUCCEED(gr_hypgeom_2f1(x, x, y, k2, k1, 0, K2));
                GR_MUST_SUCCEED(gr_set_str(y, "sqrt(8/7)*2/pi", K2));
                GR_MUST_SUCCEED(gr_mul(y, y, w, K2));
                if (gr_equal(x, y, K2) != T_TRUE)
                {
                    flint_printf("FAIL: 2F1(1/4, 3/4; 1; 7/16)\n");
                    gr_println(x, K2);
                    flint_abort();
                }
                GR_MUST_SUCCEED(gr_set_str(x, "1/4", K2));
                GR_MUST_SUCCEED(gr_set_str(k1, "24/49", K2));
                GR_MUST_SUCCEED(gr_hypgeom_2f1(x, x, x, k2, k1, 0, K2));
                GR_MUST_SUCCEED(gr_set_str(y, "2/pi", K2));
                GR_MUST_SUCCEED(gr_mul(y, y, w, K2));
                if (gr_equal(x, y, K2) != T_TRUE)
                {
                    flint_printf("FAIL: 2F1(1/4, 1/4; 1; 24/49)\n");
                    gr_println(x, K2);
                    flint_abort();
                }
            }

            /* m = 1/3, m1, m2 down the Landen chain */
            GR_MUST_SUCCEED(gr_set_str(m, "1/3", K2));
            GR_MUST_SUCCEED(gr_set_str(k1, "(1-sqrt(2/3))/(1+sqrt(2/3))", K2));
            GR_MUST_SUCCEED(gr_sqr(m1, k1, K2));
            GR_MUST_SUCCEED(gr_sub_ui(k2, m1, 1, K2));
            GR_MUST_SUCCEED(gr_neg(k2, k2, K2));
            GR_MUST_SUCCEED(gr_sqrt(k2, k2, K2));
            GR_MUST_SUCCEED(gr_sub_ui(w, k2, 1, K2));
            GR_MUST_SUCCEED(gr_neg(w, w, K2));
            GR_MUST_SUCCEED(gr_add_ui(k2, k2, 1, K2));
            GR_MUST_SUCCEED(gr_div(k2, w, k2, K2));
            GR_MUST_SUCCEED(gr_sqr(m2, k2, K2));

            /* K(m) = (1 + k1) (1 + k2) K(m2) */
            if (order == 0)
            {
                GR_MUST_SUCCEED(gr_elliptic_k(y, m2, K2));
                GR_MUST_SUCCEED(gr_elliptic_k(x, m, K2));
            }
            else
            {
                GR_MUST_SUCCEED(gr_elliptic_k(x, m, K2));
                GR_MUST_SUCCEED(gr_elliptic_k(y, m2, K2));
            }
            GR_MUST_SUCCEED(gr_add_ui(w, k1, 1, K2));
            GR_MUST_SUCCEED(gr_mul(y, y, w, K2));
            GR_MUST_SUCCEED(gr_add_ui(w, k2, 1, K2));
            GR_MUST_SUCCEED(gr_mul(y, y, w, K2));
            if (gr_equal(x, y, K2) != T_TRUE)
            {
                flint_printf("FAIL: Landen (K, order %wd)\n", order);
                gr_println(x, K2);
                gr_println(y, K2);
                flint_abort();
            }

            /* E(m) = (1 + k') E(m1) - k' K(m), k' = sqrt(2/3) */
            if (order == 0)
            {
                GR_MUST_SUCCEED(gr_elliptic_e(y, m1, K2));
                GR_MUST_SUCCEED(gr_elliptic_e(x, m, K2));
            }
            else
            {
                GR_MUST_SUCCEED(gr_elliptic_e(x, m, K2));
                GR_MUST_SUCCEED(gr_elliptic_e(y, m1, K2));
            }
            GR_MUST_SUCCEED(gr_set_str(w, "1+sqrt(2/3)", K2));
            GR_MUST_SUCCEED(gr_mul(y, y, w, K2));
            GR_MUST_SUCCEED(gr_elliptic_k(w, m, K2));
            GR_MUST_SUCCEED(gr_set_str(k2, "sqrt(2/3)", K2));
            GR_MUST_SUCCEED(gr_mul(w, w, k2, K2));
            GR_MUST_SUCCEED(gr_sub(y, y, w, K2));
            if (gr_equal(x, y, K2) != T_TRUE)
            {
                flint_printf("FAIL: Landen (E, order %wd)\n", order);
                gr_println(x, K2);
                gr_println(y, K2);
                flint_abort();
            }

            GR_TMP_CLEAR4(m, m1, m2, k1, K2);
            GR_TMP_CLEAR4(k2, x, y, w, K2);
            gr_ctx_clear(K2);
        }
    }

    /* Landen's E relation with K(m) not yet adjoined when the two
       E values meet: E(m) = (1 + s) E(k1^2) - s K(m), s = sqrt(1 - m),
       k1 = (1 - s)/(1 + s), at a complex m */
    {
        gr_ctx_t K2;
        gr_ptr m, s, k1, x, y, w;

        gr_ctx_init_tower_lazy(K2, QQ, GR_TOWER_MERGE_EXPRESS);
        GR_TMP_INIT4(m, s, k1, x, K2);
        GR_TMP_INIT2(y, w, K2);

        /* an earlier computation in the same context (this changes the
           order in which the relations are found) */
        GR_MUST_SUCCEED(gr_set_str(s, "1/3", K2));
        GR_MUST_SUCCEED(gr_set_str(k1, "1/5", K2));
        GR_MUST_SUCCEED(gr_set_str(m, "3/7", K2));
        GR_MUST_SUCCEED(gr_set_str(w, "2*i", K2));
        GR_MUST_SUCCEED(gr_hypgeom_2f1(y, s, k1, m, w, 0, K2));

        GR_MUST_SUCCEED(gr_set_str(m, "2/7 + i/3", K2));
        GR_MUST_SUCCEED(gr_sub_ui(s, m, 1, K2));
        GR_MUST_SUCCEED(gr_neg(s, s, K2));
        GR_MUST_SUCCEED(gr_sqrt(s, s, K2));
        GR_MUST_SUCCEED(gr_sub_ui(w, s, 1, K2));
        GR_MUST_SUCCEED(gr_neg(w, w, K2));
        GR_MUST_SUCCEED(gr_add_ui(k1, s, 1, K2));
        GR_MUST_SUCCEED(gr_div(k1, w, k1, K2));
        GR_MUST_SUCCEED(gr_sqr(k1, k1, K2));

        GR_MUST_SUCCEED(gr_elliptic_e(x, m, K2));
        GR_MUST_SUCCEED(gr_elliptic_e(y, k1, K2));
        GR_MUST_SUCCEED(gr_add_ui(w, s, 1, K2));
        GR_MUST_SUCCEED(gr_mul(y, y, w, K2));
        GR_MUST_SUCCEED(gr_sub(x, x, y, K2));
        GR_MUST_SUCCEED(gr_elliptic_k(w, m, K2));
        GR_MUST_SUCCEED(gr_mul(w, w, s, K2));
        GR_MUST_SUCCEED(gr_add(x, x, w, K2));
        if (gr_is_zero(x, K2) != T_TRUE)
        {
            flint_printf("FAIL: Landen (E, complex m)\n");
            gr_println(x, K2);
            flint_abort();
        }

        GR_TMP_CLEAR4(m, s, k1, x, K2);
        GR_TMP_CLEAR2(y, w, K2);
        gr_ctx_clear(K2);
    }

    gr_ctx_clear(K);
    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
