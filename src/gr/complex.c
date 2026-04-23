/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr.h"
#include "gr_generic.h"
#include "gr_vec.h"
#include "gr_special.h"
#include "gr_generic.h"

typedef struct
{
    gr_ctx_struct * real_ctx;
    int flags;
} _gr_complex_ctx_struct;

typedef gr_ctx_struct gr_complex_ctx_struct;
typedef gr_complex_ctx_struct gr_complex_ctx_t[1];


#define GR_COMPLEX_CTX(ring_ctx) ((_gr_complex_ctx_struct *)((ring_ctx)))
#define GR_COMPLEX_REAL_CTX(ring_ctx) (GR_COMPLEX_CTX(ring_ctx)->real_ctx)
#define GR_COMPLEX_FLAGS(ring_ctx) (GR_COMPLEX_CTX(ring_ctx)->flags)

#define REAL_CTX GR_COMPLEX_REAL_CTX
#define RE(x, ctx) (x)
#define IM(x, ctx) GR_ENTRY((x), 1, GR_COMPLEX_REAL_CTX(ctx)->sizeof_elem)


static int
_gr_complex_ctx_write(gr_stream_t out, gr_complex_ctx_t ctx)
{
    gr_stream_write(out, "Complex algebra over ");
    gr_ctx_write(out, REAL_CTX(ctx));
    return GR_SUCCESS;
}

static void
_gr_complex_ctx_clear(gr_ctx_t ctx)
{
}

static truth_t _gr_complex_ctx_is_ring(gr_ctx_t ctx) { return gr_ctx_is_ring(REAL_CTX(ctx)); }
static truth_t _gr_complex_ctx_is_commutative_ring(gr_ctx_t ctx) { return gr_ctx_is_commutative_ring(REAL_CTX(ctx)); }
static truth_t _gr_complex_ctx_is_integral_domain(gr_ctx_t ctx) { return gr_ctx_is_integral_domain(REAL_CTX(ctx)); }
static truth_t _gr_complex_ctx_is_field(gr_ctx_t ctx) { return gr_ctx_is_field(REAL_CTX(ctx)); }
static truth_t _gr_complex_ctx_is_rational_vector_space(gr_ctx_t ctx) { return gr_ctx_is_rational_vector_space(REAL_CTX(ctx)); }
static truth_t _gr_complex_ctx_is_real_vector_space(gr_ctx_t ctx) { return gr_ctx_is_real_vector_space(REAL_CTX(ctx)); }
static truth_t _gr_complex_ctx_is_approx_commutative_ring(gr_ctx_t ctx) { return gr_ctx_is_approx_commutative_ring(REAL_CTX(ctx)); }

static truth_t _gr_complex_ctx_is_threadsafe(gr_ctx_t ctx) { return gr_ctx_is_threadsafe(REAL_CTX(ctx)); }
static truth_t _gr_complex_ctx_is_finite(gr_ctx_t ctx) { return gr_ctx_is_finite(REAL_CTX(ctx)); }
static truth_t _gr_complex_ctx_is_finite_characteristic(gr_ctx_t ctx) { return gr_ctx_is_finite_characteristic(REAL_CTX(ctx)); }
static truth_t _gr_complex_ctx_is_exact(gr_ctx_t ctx) { return gr_ctx_is_exact(REAL_CTX(ctx)); }

static int
_gr_complex_ctx_gen_name(char ** name, slong i, gr_ctx_t ctx)
{
    if (i != 0)
        return GR_DOMAIN;

    *name = flint_malloc(2);
    if (*name == NULL)
        return GR_UNABLE;
    (*name)[0] = 'I';
    (*name)[1] = '\0';
    return GR_SUCCESS;
}

static gr_ptr _gr_complex_ctx_base(gr_ctx_t ctx) { return REAL_CTX(ctx); }

static void
_gr_complex_init(gr_ptr x, gr_complex_ctx_t ctx)
{
    gr_ptr a = RE(x, ctx), b = IM(x, ctx);
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);

    gr_init(a, real_ctx);
    gr_init(b, real_ctx);
}

static void
_gr_complex_clear(gr_ptr x, gr_complex_ctx_t ctx)
{
    gr_ptr a = RE(x, ctx), b = IM(x, ctx);
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);

    gr_clear(a, real_ctx);
    gr_clear(b, real_ctx);
}

static int
_gr_complex_write(gr_stream_t out, gr_srcptr x, gr_complex_ctx_t ctx)
{
    gr_srcptr a = RE(x, ctx), b = IM(x, ctx);
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_stream_write(out, "(");
    status |= gr_write(out, a, real_ctx);
    status |= gr_stream_write(out, ") + (");
    status |= gr_write(out, b, real_ctx);
    status |= gr_stream_write(out, ") * I");

    return status;
}

static void
_gr_complex_set_shallow(gr_ptr x, gr_srcptr y, gr_complex_ctx_t ctx)
{
    gr_ptr a = RE(x, ctx), b = IM(x, ctx);
    gr_srcptr c = RE(y, ctx), d = IM(y, ctx);
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);

    gr_set_shallow(a, c, real_ctx);
    gr_set_shallow(b, d, real_ctx);
}

static void
_gr_complex_swap(gr_ptr x, gr_ptr y, gr_complex_ctx_t ctx)
{
    gr_ptr a = RE(x, ctx), b = IM(x, ctx);
    gr_ptr c = RE(y, ctx), d = IM(y, ctx);
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);

    gr_swap(a, c, real_ctx);
    gr_swap(b, d, real_ctx);
}

static int
_gr_complex_randtest(gr_ptr x, flint_rand_t state, gr_complex_ctx_t ctx)
{
    gr_ptr a = RE(x, ctx), b = IM(x, ctx);
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_randtest(a, state, real_ctx);
    status |= n_randint(state, 2) ? gr_zero(b, real_ctx) : gr_randtest(b, state, real_ctx);

    return status;
}

static int
_gr_complex_zero(gr_ptr x, gr_complex_ctx_t ctx)
{
    gr_ptr a = RE(x, ctx), b = IM(x, ctx);
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_zero(a, real_ctx);
    status |= gr_zero(b, real_ctx);
    return status;
}

static int
_gr_complex_one(gr_ptr x, gr_complex_ctx_t ctx)
{
    gr_ptr a = RE(x, ctx), b = IM(x, ctx);
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_one(a, real_ctx);
    status |= gr_zero(b, real_ctx);
    return status;
}

static int
_gr_complex_neg_one(gr_ptr x, gr_complex_ctx_t ctx)
{
    gr_ptr a = RE(x, ctx), b = IM(x, ctx);
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_neg_one(a, real_ctx);
    status |= gr_zero(b, real_ctx);
    return status;
}

static int
_gr_complex_set(gr_ptr res, gr_srcptr x, gr_complex_ctx_t ctx)
{
    gr_ptr a = RE(res, ctx), b = IM(res, ctx);
    gr_srcptr c = RE(x, ctx), d = IM(x, ctx);
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set(a, c, real_ctx);
    status |= gr_set(b, d, real_ctx);
    return status;
}

static int
_gr_complex_set_si(gr_ptr res, slong x, gr_complex_ctx_t ctx)
{
    gr_ptr a = RE(res, ctx), b = IM(res, ctx);
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set_si(a, x, real_ctx);
    status |= gr_zero(b, real_ctx);
    return status;
}

static int
_gr_complex_set_ui(gr_ptr res, ulong x, gr_complex_ctx_t ctx)
{
    gr_ptr a = RE(res, ctx), b = IM(res, ctx);
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set_ui(a, x, real_ctx);
    status |= gr_zero(b, real_ctx);
    return status;
}

static int
_gr_complex_set_fmpz(gr_ptr res, const fmpz_t x, gr_complex_ctx_t ctx)
{
    gr_ptr a = RE(res, ctx), b = IM(res, ctx);
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set_fmpz(a, x, real_ctx);
    status |= gr_zero(b, real_ctx);
    return status;
}

static int
_gr_complex_set_other(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_complex_ctx_t ctx)
{
    gr_ptr a = RE(res, ctx), b = IM(res, ctx);
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);
    int status = GR_SUCCESS;

    if (x_ctx == real_ctx)
    {
        status |= gr_set(a, x, real_ctx);
        status |= gr_zero(b, real_ctx);
    }
    else if (x_ctx->which_ring == GR_CTX_GR_COMPLEX)
    {
        gr_srcptr c = RE(x, x_ctx), d = IM(x, x_ctx);
        gr_ctx_struct * x_real_ctx = REAL_CTX(x_ctx);

        status |= gr_set_other(a, c, x_real_ctx, real_ctx);
        status |= gr_set_other(b, d, x_real_ctx, real_ctx);
    }
    else
    {
        status = gr_set_other(a, x, x_ctx, real_ctx);

        if (status == GR_SUCCESS)
        {
            status = gr_zero(b, real_ctx);
        }
        else
        {
            gr_ptr re, im;

            GR_TMP_INIT2(re, im, x_ctx);

            status = gr_re(re, x, x_ctx);
            status |= gr_im(im, x, x_ctx);

            if (status == GR_SUCCESS)
            {
                status = gr_set_other(a, re, x_ctx, real_ctx);
                status |= gr_set_other(b, im, x_ctx, real_ctx);
            }

            GR_TMP_CLEAR2(re, im, x_ctx);
        }
    }

    return status;
}

static int
_gr_complex_i(gr_ptr x, gr_complex_ctx_t ctx)
{
    gr_ptr a = RE(x, ctx), b = IM(x, ctx);
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_zero(a, real_ctx);
    status |= gr_one(b, real_ctx);
    return status;
}

static int
_gr_complex_gens_recursive(gr_vec_t vec, gr_ctx_t ctx)
{
    int status;
    gr_vec_t vec1;
    slong i, n;

    /* Get generators of the real ring */
    gr_vec_init(vec1, 0, REAL_CTX(ctx));
    status = gr_gens_recursive(vec1, REAL_CTX(ctx));
    n = vec1->length;

    gr_vec_set_length(vec, n + 1, ctx);

    /* Promote to complex */
    for (i = 0; i < n; i++)
    {
        gr_ptr x = gr_vec_entry_ptr(vec, i, ctx);
        gr_srcptr y = gr_vec_entry_srcptr(vec1, i, REAL_CTX(ctx));

        status |= gr_set(RE(x, ctx), y, REAL_CTX(ctx));
    }

    status |= _gr_complex_i(gr_vec_entry_ptr(vec, n, ctx), ctx);

    gr_vec_clear(vec1, REAL_CTX(ctx));

    return status;
}

static truth_t
_gr_complex_equal(gr_srcptr x, gr_srcptr y, gr_complex_ctx_t ctx)
{
    gr_srcptr a, b, c, d;
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);

    a = RE(x, ctx);
    b = IM(x, ctx);
    c = RE(y, ctx);
    d = IM(y, ctx);

    truth_t eq1, eq2;

    eq1 = gr_equal(a, c, real_ctx);
    if (eq1 == T_FALSE)
        return eq1;

    eq2 = gr_equal(b, d, real_ctx);

    return truth_and(eq1, eq2);
}

static truth_t
_gr_complex_is_zero(gr_srcptr x, gr_complex_ctx_t ctx)
{
    gr_srcptr a, b;
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);

    a = RE(x, ctx);
    b = IM(x, ctx);

    truth_t eq1, eq2;
    eq1 = gr_is_zero(a, real_ctx);
    if (eq1 == T_FALSE)
        return eq1;

    eq2 = gr_is_zero(b, real_ctx);
    return truth_and(eq1, eq2);
}

static truth_t
_gr_complex_is_one(gr_srcptr x, gr_complex_ctx_t ctx)
{
    gr_srcptr a, b;
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);

    a = RE(x, ctx);
    b = IM(x, ctx);

    truth_t eq1, eq2;
    eq1 = gr_is_one(a, real_ctx);
    if (eq1 == T_FALSE)
        return eq1;

    eq2 = gr_is_zero(b, real_ctx);
    return truth_and(eq1, eq2);
}

static truth_t
_gr_complex_is_neg_one(gr_srcptr x, gr_complex_ctx_t ctx)
{
    gr_srcptr a, b;
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);

    a = RE(x, ctx);
    b = IM(x, ctx);

    truth_t eq1, eq2;
    eq1 = gr_is_neg_one(a, real_ctx);
    if (eq1 == T_FALSE)
        return eq1;

    eq2 = gr_is_zero(b, real_ctx);
    return truth_and(eq1, eq2);
}


static int
_gr_complex_neg(gr_ptr res, gr_srcptr x, gr_complex_ctx_t ctx)
{
    gr_srcptr a = RE(x, ctx), b = IM(x, ctx);
    gr_ptr c = RE(res, ctx), d = IM(res, ctx);
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_neg(c, a, real_ctx);
    status |= gr_neg(d, b, real_ctx);
    return status;
}

static int
_gr_complex_add(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_complex_ctx_t ctx)
{
    gr_srcptr a = RE(x, ctx), b = IM(x, ctx);
    gr_srcptr c = RE(y, ctx), d = IM(y, ctx);
    gr_ptr e = RE(res, ctx), f = IM(res, ctx);
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_add(e, a, c, real_ctx);
    status |= gr_add(f, b, d, real_ctx);

    return status;
}

static int
_gr_complex_sub(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_complex_ctx_t ctx)
{
    gr_srcptr a = RE(x, ctx), b = IM(x, ctx);
    gr_srcptr c = RE(y, ctx), d = IM(y, ctx);
    gr_ptr e = RE(res, ctx), f = IM(res, ctx);
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_sub(e, a, c, real_ctx);
    status |= gr_sub(f, b, d, real_ctx);

    return status;
}

static int
_gr_complex_mul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_complex_ctx_t ctx)
{
    gr_srcptr a = RE(x, ctx), b = IM(x, ctx);
    gr_srcptr c = RE(y, ctx), d = IM(y, ctx);
    gr_ptr e = RE(res, ctx), f = IM(res, ctx);
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);
    int status = GR_SUCCESS;

    /* todo: squaring if x == y */

    if (res == x || res == y)
    {
        gr_ptr t, u, v, w;
        GR_TMP_INIT4(t, u, v, w, real_ctx);

        status |= gr_mul(t, a, c, real_ctx);
        status |= gr_mul(u, b, d, real_ctx);
        status |= gr_mul(v, a, d, real_ctx);
        status |= gr_mul(w, b, c, real_ctx);

        status |= gr_sub(e, t, u, real_ctx);
        status |= gr_add(f, v, w, real_ctx);

        GR_TMP_CLEAR4(t, u, v, w, real_ctx);
    }
    else
    {
        gr_ptr t, u;
        GR_TMP_INIT2(t, u, real_ctx);

        status |= gr_mul(t, a, c, real_ctx);
        status |= gr_mul(u, b, d, real_ctx);
        status |= gr_sub(e, t, u, real_ctx);

        status |= gr_mul(t, a, d, real_ctx);
        status |= gr_mul(u, b, c, real_ctx);
        status |= gr_add(f, t, u, real_ctx);

        GR_TMP_CLEAR2(t, u, real_ctx);
    }

    return status;
}

static int
_gr_complex_inv(gr_ptr res, gr_srcptr x, gr_complex_ctx_t ctx)
{
    gr_srcptr a = RE(x, ctx), b = IM(x, ctx);
    gr_ptr c = RE(res, ctx), d = IM(res, ctx);
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);
    int status = GR_SUCCESS;

    gr_ptr t, u;
    GR_TMP_INIT2(t, u, real_ctx);

    status |= gr_sqr(t, a, real_ctx);
    status |= gr_sqr(u, b, real_ctx);
    status |= gr_add(t, t, u, real_ctx);
    status |= _gr_vec_div_scalar(c, x, 2, t, real_ctx);
    status |= gr_neg(d, d, real_ctx);

    GR_TMP_CLEAR2(t, u, real_ctx);

    return status;
}

static int
_gr_complex_div(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_complex_ctx_t ctx)
{
    gr_srcptr a = RE(x, ctx), b = IM(x, ctx);
    gr_srcptr c = RE(y, ctx), d = IM(y, ctx);
    gr_ptr e = RE(res, ctx), f = IM(res, ctx);
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);
    int status = GR_SUCCESS;

    gr_ptr s, t, u, v, w;
    GR_TMP_INIT5(s, t, u, v, w, real_ctx);

    status |= gr_sqr(s, c, real_ctx);
    status |= gr_sqr(u, d, real_ctx);
    status |= gr_add(s, s, u, real_ctx);

    status |= gr_mul(t, a, c, real_ctx);
    status |= gr_mul(u, b, d, real_ctx);
    status |= gr_mul(v, a, d, real_ctx);
    status |= gr_mul(w, b, c, real_ctx);

    status |= gr_add(e, t, u, real_ctx);
    status |= gr_sub(f, w, v, real_ctx);

    status |= _gr_vec_div_scalar(e, e, 2, s, real_ctx);

    GR_TMP_CLEAR5(s, t, u, v, w, real_ctx);

    return status;
}

static int
_gr_complex_abs(gr_ptr res, gr_srcptr x, gr_complex_ctx_t ctx)
{
    gr_srcptr a = RE(x, ctx), b = IM(x, ctx);
    gr_ptr c = RE(res, ctx), d = IM(res, ctx);
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);
    int status = GR_SUCCESS;

    if (gr_is_zero(b, ctx) == T_TRUE)
    {
        status |= gr_abs(c, a, real_ctx);
    }
    else if (gr_is_zero(a, ctx) == T_TRUE)
    {
        status |= gr_abs(c, b, real_ctx);
    }
    else
    {
        status |= gr_sqr(c, a, real_ctx);
        status |= gr_addmul(c, b, b, real_ctx);
        status |= gr_sqrt(c, c, real_ctx);
    }
    
    status |= gr_zero(d, real_ctx);

    return status;
}

static int
_gr_complex_re(gr_ptr res, gr_srcptr x, gr_complex_ctx_t ctx)
{
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set(RE(res, ctx), RE(x, ctx), real_ctx);
    status |= gr_zero(IM(res, ctx), real_ctx);
    return status;
}

static int
_gr_complex_im(gr_ptr res, gr_srcptr x, gr_complex_ctx_t ctx)
{
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set(RE(res, ctx), IM(x, ctx), real_ctx);
    status |= gr_zero(IM(res, ctx), real_ctx);
    return status;
}

static int
_gr_complex_conj(gr_ptr res, gr_srcptr x, gr_complex_ctx_t ctx)
{
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set(RE(res, ctx), RE(x, ctx), real_ctx);
    status |= gr_neg(IM(res, ctx), IM(x, ctx), real_ctx);
    return status;
}

static int
_gr_complex_pi(gr_ptr x, gr_complex_ctx_t ctx)
{
    gr_ptr a = RE(x, ctx), b = IM(x, ctx);
    gr_ctx_struct * real_ctx = REAL_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_pi(a, real_ctx);
    status |= gr_zero(b, real_ctx);
    return status;
}


int _gr_complex_methods_initialized = 0;

gr_static_method_table _gr_complex_methods;

gr_method_tab_input _gr_complex_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_complex_ctx_write},
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_complex_ctx_clear},

    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) _gr_complex_ctx_is_ring},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr)  _gr_complex_ctx_is_commutative_ring},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr)  _gr_complex_ctx_is_integral_domain},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr)  _gr_complex_ctx_is_field},
    {GR_METHOD_CTX_IS_RATIONAL_VECTOR_SPACE, (gr_funcptr) _gr_complex_ctx_is_rational_vector_space},
    {GR_METHOD_CTX_IS_REAL_VECTOR_SPACE, (gr_funcptr) _gr_complex_ctx_is_real_vector_space},
    {GR_METHOD_CTX_IS_COMPLEX_VECTOR_SPACE, (gr_funcptr) _gr_complex_ctx_is_real_vector_space},
    {GR_METHOD_CTX_IS_THREADSAFE,       (gr_funcptr)  _gr_complex_ctx_is_threadsafe},
    {GR_METHOD_CTX_IS_FINITE,           (gr_funcptr) _gr_complex_ctx_is_finite},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,    (gr_funcptr) _gr_complex_ctx_is_finite_characteristic},
    {GR_METHOD_CTX_IS_APPROX_COMMUTATIVE_RING, (gr_funcptr)  _gr_complex_ctx_is_approx_commutative_ring},

    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) _gr_complex_ctx_is_exact},
    {GR_METHOD_CTX_NGENS,       (gr_funcptr) gr_generic_ctx_ngens_1},
    {GR_METHOD_CTX_GEN_NAME,    (gr_funcptr) _gr_complex_ctx_gen_name},
    {GR_METHOD_CTX_BASE,                (gr_funcptr) _gr_complex_ctx_base},
    {GR_METHOD_INIT,            (gr_funcptr) _gr_complex_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_complex_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_complex_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _gr_complex_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_complex_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_complex_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_complex_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_complex_one},
    {GR_METHOD_NEG_ONE,         (gr_funcptr) _gr_complex_neg_one},
    {GR_METHOD_GEN,             (gr_funcptr) _gr_complex_i},
    {GR_METHOD_GENS,           (gr_funcptr) gr_generic_gens_single},
    {GR_METHOD_GENS_RECURSIVE,  (gr_funcptr) _gr_complex_gens_recursive},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_complex_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_complex_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_complex_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_complex_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_complex_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_complex_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_complex_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_complex_set_fmpz},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_complex_set_other},
    {GR_METHOD_SET_STR,     (gr_funcptr) gr_generic_set_str_balance_additions},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_complex_neg},

    {GR_METHOD_ADD,             (gr_funcptr) _gr_complex_add},
/*
    {GR_METHOD_ADD_UI,          (gr_funcptr) _gr_complex_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) _gr_complex_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) _gr_complex_add_fmpz},
    {GR_METHOD_ADD_FMPQ,        (gr_funcptr) _gr_complex_add_fmpq},
*/
    {GR_METHOD_SUB,             (gr_funcptr) _gr_complex_sub},
/*
    {GR_METHOD_SUB_UI,          (gr_funcptr) _gr_complex_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) _gr_complex_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) _gr_complex_sub_fmpz},
    {GR_METHOD_SUB_FMPQ,        (gr_funcptr) _gr_complex_sub_fmpq},
*/
    {GR_METHOD_MUL,             (gr_funcptr) _gr_complex_mul},
/*
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_complex_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_complex_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_complex_mul_fmpz},
    {GR_METHOD_MUL_FMPQ,        (gr_funcptr) _gr_complex_mul_fmpq},

    {GR_METHOD_MUL_TWO,         (gr_funcptr) _gr_complex_mul_two},
    {GR_METHOD_SQR,             (gr_funcptr) _gr_complex_sqr},
*/
    {GR_METHOD_DIV,             (gr_funcptr) _gr_complex_div},
/*
    {GR_METHOD_DIV_UI,          (gr_funcptr) _gr_complex_div_ui},
    {GR_METHOD_DIV_SI,          (gr_funcptr) _gr_complex_div_si},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) _gr_complex_div_fmpz},
    {GR_METHOD_DIV_FMPQ,        (gr_funcptr) _gr_complex_div_fmpq},
*/

/*
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_complex_is_invertible},
*/
    {GR_METHOD_INV,             (gr_funcptr) _gr_complex_inv},
/*
    {GR_METHOD_SQRT,             (gr_funcptr) _gr_complex_sqrt},

    {GR_METHOD_NUMERATOR,       (gr_funcptr) _gr_complex_numerator},
    {GR_METHOD_DENOMINATOR,     (gr_funcptr) _gr_complex_denominator},
*/

    {GR_METHOD_ABS,               (gr_funcptr) _gr_complex_abs},
    {GR_METHOD_RE,               (gr_funcptr) _gr_complex_re},
    {GR_METHOD_IM,               (gr_funcptr) _gr_complex_im},
    {GR_METHOD_CONJ,             (gr_funcptr) _gr_complex_conj},

    {GR_METHOD_I,               (gr_funcptr) _gr_complex_i},
    {GR_METHOD_PI,              (gr_funcptr) _gr_complex_pi},

    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_gr_complex(gr_ctx_t ctx, gr_ctx_t real_ctx)
{
    ctx->which_ring = GR_CTX_GR_COMPLEX;
    ctx->sizeof_elem = 2 * real_ctx->sizeof_elem;
    ctx->size_limit = WORD_MAX;

    GR_COMPLEX_REAL_CTX(ctx) = real_ctx;
    GR_COMPLEX_FLAGS(ctx) = 0;

    ctx->methods = _gr_complex_methods;

    if (!_gr_complex_methods_initialized)
    {
        gr_method_tab_init(_gr_complex_methods, _gr_complex_methods_input);
        _gr_complex_methods_initialized = 1;
    }
}

