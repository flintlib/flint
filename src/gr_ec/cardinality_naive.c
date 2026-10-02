/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include "fmpz.h"
#include "gr.h"
#include "gr_ec.h"
#include "impl.h"

/* ------------------------------------------------------------------ */
/* walking over the elements of a finite field                        */
/* ------------------------------------------------------------------ */

/*
    F_q = F_p[g] with g the generator, so the elements are exactly the
    sums c_0 + c_1 g + ... + c_{d-1} g^{d-1} with c_i in [0, p). We walk
    over them with an odometer on the digits.
*/
typedef struct
{
    gr_ctx_struct * R;
    gr_ptr pows;            /* 1, g, ..., g^(d-1) */
    ulong * digits;
    slong d;
    ulong p;
}
field_iter_struct;

static void
field_iter_clear(field_iter_struct * it)
{
    if (it->pows != NULL)
        gr_heap_clear_vec(it->pows, it->d, it->R);
    flint_free(it->digits);
}

static int
field_iter_init(field_iter_struct * it, gr_ctx_t R)
{
    fmpz_t p;
    slong i, sz = R->sizeof_elem;
    int status = GR_SUCCESS;

    it->R = R;
    it->pows = NULL;
    it->digits = NULL;

    fmpz_init(p);

    if (gr_ctx_fq_prime(p, R) != GR_SUCCESS || !fmpz_abs_fits_ui(p)
            || gr_ctx_fq_degree(&it->d, R) != GR_SUCCESS || it->d < 1)
    {
        fmpz_clear(p);
        return GR_UNABLE;
    }

    it->p = fmpz_get_ui(p);
    fmpz_clear(p);

    it->digits = flint_calloc(it->d, sizeof(ulong));
    it->pows = gr_heap_init_vec(it->d, R);

    status |= gr_one(it->pows, R);

    if (it->d > 1)
    {
        gr_ptr g;
        GR_TMP_INIT(g, R);
        status |= gr_gen(g, R);

        for (i = 1; i < it->d; i++)
            status |= gr_mul(GR_ENTRY(it->pows, i, sz),
                        GR_ENTRY(it->pows, i - 1, sz), g, R);

        GR_TMP_CLEAR(g, R);
    }

    if (status != GR_SUCCESS)
        field_iter_clear(it);

    return status;
}

/* the element for the current digits */
static int
field_iter_get(gr_ptr res, const field_iter_struct * it)
{
    gr_ctx_struct * R = it->R;
    slong i, sz = R->sizeof_elem;
    int status = GR_SUCCESS;

    if (it->d == 1)
        return gr_set_ui(res, it->digits[0], R);

    status |= gr_zero(res, R);

    for (i = 0; i < it->d; i++)
        if (it->digits[i] != 0)
        {
            gr_ptr t;
            GR_TMP_INIT(t, R);
            status |= gr_mul_ui(t, GR_ENTRY(it->pows, i, sz), it->digits[i], R);
            status |= gr_add(res, res, t, R);
            GR_TMP_CLEAR(t, R);
        }

    return status;
}

/* returns 0 once it has wrapped around, that is, once every element is done */
static int
field_iter_next(field_iter_struct * it)
{
    slong i;

    for (i = 0; i < it->d; i++)
    {
        it->digits[i]++;

        if (it->digits[i] < it->p)
            return 1;

        it->digits[i] = 0;
    }

    return 0;
}

/* ------------------------------------------------------------------ */
/* naive counting                                                     */
/* ------------------------------------------------------------------ */

/* GR_EC_NAIVE_MAX_Q, the largest field this will walk, is in impl.h */

int
gr_ec_ctx_cardinality_naive(fmpz_t res, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    field_iter_struct it;
    fmpz_t q, p, count;
    gr_ptr x, b, c, t;
    slong d;
    int char_two, status = GR_SUCCESS;

    if (gr_ctx_is_field(R) != T_TRUE)
        return GR_DOMAIN;

    fmpz_init(q);
    fmpz_init(p);
    fmpz_init(count);

    if (gr_ctx_cardinality_fmpz(q, R) != GR_SUCCESS
            || gr_ctx_fq_prime(p, R) != GR_SUCCESS
            || !fmpz_fits_si(q) || fmpz_cmp_si(q, GR_EC_NAIVE_MAX_Q) > 0)
    {
        fmpz_clear(q); fmpz_clear(p); fmpz_clear(count);
        return GR_UNABLE;
    }

    char_two = (fmpz_cmp_ui(p, 2) == 0);

    if (gr_ctx_fq_degree(&d, R) != GR_SUCCESS)
    {
        fmpz_clear(q); fmpz_clear(p); fmpz_clear(count);
        return GR_UNABLE;
    }

    status = field_iter_init(&it, R);

    if (status != GR_SUCCESS)
    {
        fmpz_clear(q); fmpz_clear(p); fmpz_clear(count);
        return status;
    }

    GR_TMP_INIT4(x, b, c, t, R);

    fmpz_one(count);            /* the point at infinity */

    do
    {
        status |= field_iter_get(x, &it);

        /* b = a1 x + a3 */
        status |= gr_mul(b, GR_EC_A1(ctx), x, R);
        status |= gr_add(b, b, GR_EC_A3(ctx), R);

        /* c = x^3 + a2 x^2 + a4 x + a6, by Horner */
        status |= gr_add(c, x, GR_EC_A2(ctx), R);
        status |= gr_mul(c, c, x, R);
        status |= gr_add(c, c, GR_EC_A4(ctx), R);
        status |= gr_mul(c, c, x, R);
        status |= gr_add(c, c, GR_EC_A6(ctx), R);

        if (status != GR_SUCCESS)
            break;

        /* number of y with y^2 + b y - c = 0 */
        if (!char_two)
        {
            /* the discriminant of the quadratic, b^2 + 4c */
            status |= gr_sqr(t, b, R);
            status |= gr_mul_ui(c, c, 4, R);
            status |= gr_add(t, t, c, R);

            if (status != GR_SUCCESS)
                break;

            if (gr_is_zero(t, R) == T_TRUE)
                fmpz_add_ui(count, count, 1);
            else if (gr_is_square(t, R) == T_TRUE)
                fmpz_add_ui(count, count, 2);
        }
        else if (gr_is_zero(b, R) == T_TRUE)
        {
            /* y^2 = c, and squaring is a bijection in characteristic 2 */
            fmpz_add_ui(count, count, 1);
        }
        else
        {
            /*
                y = b z turns y^2 + b y = c into z^2 + z = c / b^2, which
                is soluble exactly when the absolute trace vanishes, and
                then has two solutions.
            */
            fmpz_t tr;

            status |= gr_sqr(t, b, R);
            status |= gr_div(t, c, t, R);

            fmpz_init(tr);

            if (d == 1)
                status |= gr_get_fmpz(tr, t, R);
            else
                status |= gr_fq_trace(tr, t, R);

            if (status == GR_SUCCESS && fmpz_is_even(tr))
                fmpz_add_ui(count, count, 2);

            fmpz_clear(tr);

            if (status != GR_SUCCESS)
                break;
        }
    }
    while (field_iter_next(&it));

    if (status == GR_SUCCESS)
        fmpz_set(res, count);

    GR_TMP_CLEAR4(x, b, c, t, R);
    field_iter_clear(&it);
    fmpz_clear(q);
    fmpz_clear(p);
    fmpz_clear(count);

    return status;
}
