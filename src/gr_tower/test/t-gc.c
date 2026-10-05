/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "gr_poly.h"
#include "gr_vec.h"
#include "gr_special.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"

/*
    Collection of unused towers in lazy fields: a pool of elements is
    repeatedly updated with random operations (square roots, polynomial
    roots, exp, log, arithmetic), the replaced elements being cleared.
    The elements of the pool keep their values (checked against
    enclosures recorded when they were computed, and by exact
    comparisons with recomputed values), their printed forms read back,
    and the number of towers stays bounded.
*/

#define POOL 6

TEST_FUNCTION_START(gr_tower_gc, state)
{
    gr_ctx_t QQ;
    slong iter;

    gr_ctx_init_fmpq(QQ);

    for (iter = 0; iter < 2 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t K;
        gr_ptr pool, t, u;
        acb_ptr enc;
        slong i, step, steps = 40, max_towers = 0;
        int flags = n_randint(state, 2) ? GR_TOWER_MERGE_EXPRESS : 0;

        if (n_randint(state, 3) == 0)
            flags |= GR_TOWER_LAZY_REAL;
        gr_ctx_init_tower_lazy(K, QQ, flags);

        pool = gr_heap_init_vec(POOL, K);
        enc = _acb_vec_init(POOL);
        GR_TMP_INIT2(t, u, K);

        for (i = 0; i < POOL; i++)
        {
            GR_MUST_SUCCEED(gr_set_si(GR_ENTRY(pool, i, K->sizeof_elem), 1 + i, K));
            GR_MUST_SUCCEED(gr_tower_lazy_get_acb(enc + i, GR_ENTRY(pool, i, K->sizeof_elem), 64, K));
        }

        for (step = 0; step < steps; step++)
        {
            slong a = n_randint(state, POOL), b = n_randint(state, POOL), c = n_randint(state, POOL);
            gr_ptr x = GR_ENTRY(pool, a, K->sizeof_elem);
            gr_ptr y = GR_ENTRY(pool, b, K->sizeof_elem);
            int status = GR_SUCCESS;
            slong op = n_randint(state, 7);

            switch (op)
            {
                case 0:
                    status = gr_add_si(t, x, 1 + n_randint(state, 3), K);
                    status |= gr_sqrt(t, t, K);
                    break;
                case 1:
                    status = gr_mul(t, x, y, K);
                    status |= gr_add_si(t, t, 1, K);
                    break;
                case 2:
                {
                    /* a root of x^3 - sqrt(k) x - (1 + small) (random
                       elements of the pool would make the splitting
                       towers expensive) */
                    gr_poly_t f;
                    gr_vec_t r;
                    fmpz_vec_t m;
                    gr_ctx_t P;
                    gr_ctx_init_gr_poly(P, K);
                    gr_poly_init(f, K);
                    gr_vec_init(r, 0, K);
                    fmpz_vec_init(m, 0);
                    status = gr_poly_set_coeff_si(f, 3, 1, K);
                    status |= gr_set_si(u, 2 + n_randint(state, 5), K);
                    status |= gr_sqrt(u, u, K);
                    status |= gr_neg(u, u, K);
                    status |= gr_poly_set_coeff_scalar(f, 1, u, K);
                    status |= gr_poly_set_coeff_si(f, 0, -1 - (slong) n_randint(state, 5), K);
                    if (status == GR_SUCCESS)
                        status = gr_poly_roots(r, m, f, 0, K);
                    if (status == GR_SUCCESS && r->length > 0)
                        status = gr_set(t, gr_vec_entry_srcptr(r, n_randint(state, r->length), K), K);
                    else
                        status = GR_UNABLE;
                    gr_poly_clear(f, K);
                    gr_vec_clear(r, K);
                    fmpz_vec_clear(m);
                    gr_ctx_clear(P);
                    break;
                }
                case 3:
                    status = gr_div_si(t, x, 4, K);
                    status |= gr_exp(t, t, K);
                    break;
                case 4:
                    status = gr_mul(t, x, x, K);
                    status |= gr_add_si(t, t, 1, K);
                    status |= gr_log(t, t, K);
                    break;
                case 5:
                    status = gr_sub(t, x, y, K);
                    status |= gr_add(t, t, y, K);   /* = x */
                    break;
                default:
                    status = gr_set_si(t, 2 + n_randint(state, 30), K);
                    status |= gr_sqrt(t, t, K);
                    status |= gr_add(t, t, y, K);
                    break;
            }

            if (status != GR_SUCCESS)
                continue;

            /* keep the values of the pool small */
            {
                acb_t z;
                acb_init(z);
                if (gr_tower_lazy_get_acb(z, t, 64, K) != GR_SUCCESS)
                    status = GR_UNABLE;
                else if (!acb_is_finite(z) || arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 20) > 0 ||
                         arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 20) > 0)
                    status = GR_UNABLE;
                if (status == GR_SUCCESS)
                {
                    gr_swap(GR_ENTRY(pool, c, K->sizeof_elem), t, K);
                    acb_set(enc + c, z);
                }
                acb_clear(z);
            }

            max_towers = FLINT_MAX(max_towers, gr_tower_lazy_ctx_num_towers(K));

            /* values of the pool */
            for (i = 0; i < POOL; i++)
            {
                acb_t z;
                acb_init(z);
                GR_MUST_SUCCEED(gr_tower_lazy_get_acb(z, GR_ENTRY(pool, i, K->sizeof_elem), 64, K));
                if (!acb_overlaps(z, enc + i))
                {
                    flint_printf("FAIL: value changed (iter %wd, step %wd, element %wd)\n", iter, step, i);
                    acb_printd(z, 20); flint_printf("\n");
                    acb_printd(enc + i, 20); flint_printf("\n");
                    gr_println(GR_ENTRY(pool, i, K->sizeof_elem), K);
                    flint_abort();
                }
                acb_clear(z);
            }
        }

        /* the printed forms read back (with their definitions) */
        for (i = 0; i < POOL; i++)
        {
            char * s;
            gr_ptr x = GR_ENTRY(pool, i, K->sizeof_elem);
            GR_MUST_SUCCEED(gr_get_str(&s, x, K));
            if (gr_set_str(t, s, K) == GR_SUCCESS && gr_equal(t, x, K) == T_FALSE)
            {
                flint_printf("FAIL: printed form (iter %wd)\n%s\n", iter, s);
                flint_abort();
            }
            flint_free(s);
        }

        /* the towers of the cleared elements are collected */
        GR_MUST_SUCCEED(gr_zero(t, K));
        GR_MUST_SUCCEED(gr_zero(u, K));
        GR_MUST_SUCCEED(_gr_vec_zero(pool, POOL, K));
        {
            slong n = gr_tower_lazy_ctx_num_towers(K);
            /* (the trivial tower, and those pinned by the registries of
               canonical definitions: roots of unity, radicals of
               integers, constants such as pi) */
            if (n > 40)
            {
                flint_printf("FAIL: %wd towers remain (max %wd during the run)\n", n, max_towers);
                gr_tower_lazy_ctx_stats(K);
                flint_abort();
            }
        }

        GR_TMP_CLEAR2(t, u, K);
        gr_heap_clear_vec(pool, POOL, K);
        _acb_vec_clear(enc, POOL);
        gr_ctx_clear(K);
    }

    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
