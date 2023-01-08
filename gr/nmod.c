/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"

#define NMOD_CTX_REF(ring_ctx) (((nmod_t *)((ring_ctx))))
#define NMOD_CTX(ring_ctx) (*NMOD_CTX_REF(ring_ctx))

void
_gr_nmod_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Integers mod ");
    gr_stream_write_ui(out, NMOD_CTX(ctx).n);
    gr_stream_write(out, " (_gr_nmod)");
}

/* todo: n_is_prime is fast, but this should still be cached
   or use a fixed table lookup */
truth_t
_gr_nmod_ctx_is_field(const gr_ctx_t ctx)
{
    return n_is_prime(NMOD_CTX(ctx).n) ? T_TRUE : T_FALSE;
}

void
_gr_nmod_init(ulong * x, const gr_ctx_t ctx)
{
    x[0] = 0;
}

void
_gr_nmod_clear(ulong * x, const gr_ctx_t ctx)
{
}

void
_gr_nmod_swap(ulong * x, ulong * y, const gr_ctx_t ctx)
{
    ulong t;
    t = *x;
    *x = *y;
    *y = t;
}

int
_gr_nmod_randtest(ulong * res, flint_rand_t state, const gr_ctx_t ctx)
{
    res[0] = n_randtest(state) % NMOD_CTX(ctx).n;
    return GR_SUCCESS;
}

int
_gr_nmod_write(gr_stream_t out, const ulong * x, const gr_ctx_t ctx)
{
    gr_stream_write_si(out, x[0]);
    return GR_SUCCESS;
}

int
_gr_nmod_zero(ulong * x, const gr_ctx_t ctx)
{
    x[0] = 0;
    return GR_SUCCESS;
}

int
_gr_nmod_one(ulong * x, const gr_ctx_t ctx)
{
    x[0] = (NMOD_CTX(ctx).n != 1);
    return GR_SUCCESS;
}

int
_gr_nmod_set_si(ulong * res, slong v, const gr_ctx_t ctx)
{
    ulong t;
    nmod_t mod = NMOD_CTX(ctx);
    t = FLINT_ABS(v);
    NMOD_RED(t, t, mod);
    if (v < 0)
        t = nmod_neg(t, mod);
    res[0] = t;
    return GR_SUCCESS;
}

int
_gr_nmod_set_ui(ulong * res, ulong v, const gr_ctx_t ctx)
{
    ulong t;
    nmod_t mod = NMOD_CTX(ctx);
    NMOD_RED(t, v, mod);
    res[0] = t;
    return GR_SUCCESS;
}

int
_gr_nmod_set_fmpz(ulong * res, const fmpz_t v, const gr_ctx_t ctx)
{
    nmod_t mod = NMOD_CTX(ctx);
    res[0] = fmpz_get_nmod(v, mod);
    return GR_SUCCESS;
}

#include "flint/fmpz_mod.h"

/* todo: public interface */

typedef struct
{
    fmpz_mod_ctx_struct ctx;
    truth_t is_prime;
}
fmpz_mod_ctx_extended_struct;

#define FMPZ_MOD_CTX(ring_ctx) (&(((fmpz_mod_ctx_extended_struct *)(GR_CTX_DATA_AS_PTR(ring_ctx)))->ctx))

int
_gr_nmod_set_other(ulong * res, gr_ptr v, gr_ctx_t v_ctx, const gr_ctx_t ctx)
{
    if (v_ctx->which_ring == GR_CTX_NMOD)
    {
        if (NMOD_CTX(ctx).n != NMOD_CTX(v_ctx).n)
            return GR_DOMAIN;

        *res = *((ulong *) v);
        return GR_SUCCESS;
    }

    if (v_ctx->which_ring == GR_CTX_FMPZ_MOD)
    {
        if (!fmpz_equal_ui(FMPZ_MOD_CTX(v_ctx)->n, NMOD_CTX(ctx).n))
            return GR_DOMAIN;

        res[0] = fmpz_get_ui(v);
        return GR_SUCCESS;
    }

    return GR_UNABLE;
}

truth_t
_gr_nmod_is_zero(const ulong * x, const gr_ctx_t ctx)
{
    return (x[0] == 0) ? T_TRUE : T_FALSE;
}

truth_t
_gr_nmod_is_one(const ulong * x, const gr_ctx_t ctx)
{
    return (x[0] == (NMOD_CTX(ctx).n != 1)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_nmod_is_neg_one(const ulong * x, const gr_ctx_t ctx)
{
    return (x[0] == NMOD_CTX(ctx).n - 1) ? T_TRUE : T_FALSE;
}

truth_t
_gr_nmod_equal(const ulong * x, const ulong * y, const gr_ctx_t ctx)
{
    return (x[0] == y[0]) ? T_TRUE : T_FALSE;
}

int
_gr_nmod_set(ulong * res, const ulong * x, const gr_ctx_t ctx)
{
    res[0] = x[0];
    return GR_SUCCESS;
}

int
_gr_nmod_neg(ulong * res, const ulong * x, const gr_ctx_t ctx)
{
    res[0] = nmod_neg(x[0], NMOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nmod_add(ulong * res, const ulong * x, const ulong * y, const gr_ctx_t ctx)
{
    res[0] = nmod_add(x[0], y[0], NMOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nmod_add_si(ulong * res, const ulong * x, slong y, const gr_ctx_t ctx)
{
    ulong t;
    _gr_nmod_set_si(&t, y, ctx);
    return _gr_nmod_add(res, x, &t, ctx);
}

int
_gr_nmod_sub(ulong * res, const ulong * x, const ulong * y, const gr_ctx_t ctx)
{
    res[0] = nmod_sub(x[0], y[0], NMOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nmod_mul(ulong * res, const ulong * x, const ulong * y, const gr_ctx_t ctx)
{
    res[0] = nmod_mul(x[0], y[0], NMOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nmod_mul_si(ulong * res, const ulong * x, slong y, const gr_ctx_t ctx)
{
    ulong t;
    _gr_nmod_set_si(&t, y, ctx);
    return _gr_nmod_mul(res, x, &t, ctx);
}

int
_gr_nmod_addmul(ulong * res, const ulong * x, const ulong * y, const gr_ctx_t ctx)
{
    ulong r = res[0];
    NMOD_ADDMUL(r, x[0], y[0], NMOD_CTX(ctx));
    res[0] = r;
    return GR_SUCCESS;
}

int
_gr_nmod_submul(ulong * res, const ulong * x, const ulong * y, const gr_ctx_t ctx)
{
    ulong r = res[0];
    ulong t = nmod_neg(y[0], NMOD_CTX(ctx));
    NMOD_ADDMUL(r, x[0], t, NMOD_CTX(ctx));
    res[0] = r;
    return GR_SUCCESS;
}

int
_gr_nmod_mul_two(ulong * res, const ulong * x, const gr_ctx_t ctx)
{
    return _gr_nmod_add(res, x, x, ctx);
}

int
_gr_nmod_sqr(ulong * res, const ulong * x, const gr_ctx_t ctx)
{
    return _gr_nmod_mul(res, x, x, ctx);
}

int
_gr_nmod_inv(ulong * res, const ulong * x, const gr_ctx_t ctx)
{
    ulong r, g;

    /* todo: also handle -1 fast? */
    if (x[0] == 1)
    {
        res[0] = x[0];
        return GR_SUCCESS;
    }

    g = n_gcdinv(&r, x[0], NMOD_CTX(ctx).n);
    if (g == 1)
    {
        res[0] = r;
        return GR_SUCCESS;
    }
    else
    {
        res[0] = 0;
        return GR_DOMAIN;
    }
}

int
_gr_nmod_div(ulong * res, const ulong * x, const ulong * y, const gr_ctx_t ctx)
{
    ulong t;
    int status;

    status = _gr_nmod_inv(&t, y, ctx);
    _gr_nmod_mul(res, x, &t, ctx);
    return status;
}

int
_gr_nmod_div_si(ulong * res, const ulong * x, slong y, const gr_ctx_t ctx)
{
    ulong t;
    _gr_nmod_set_si(&t, y, ctx);
    return _gr_nmod_div(res, x, &t, ctx);
}

truth_t
_gr_nmod_is_invertible(const ulong * x, const gr_ctx_t ctx)
{
    ulong r, g;
    g = n_gcdinv(&r, x[0], NMOD_CTX(ctx).n);
    return (g == 1) ? T_TRUE : T_FALSE;
}

int
__gr_nmod_vec_dot(ulong * res, const ulong * initial, int subtract, const ulong * vec1, const ulong * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    ulong s;
    int nlimbs;
    nmod_t mod;

    if (len <= 0)
    {
        if (initial == NULL)
            _gr_nmod_zero(res, ctx);
        else
            _gr_nmod_set(res, initial, ctx);
        return GR_SUCCESS;
    }

    mod = NMOD_CTX(ctx);

    nlimbs = _nmod_vec_dot_bound_limbs(len, mod);

    NMOD_VEC_DOT(s, i, len, vec1[i], vec2[i], mod, nlimbs);

    if (initial == NULL)
    {
        if (subtract)
            s = nmod_neg(s, mod);
    }
    else
    {
        if (subtract)
            s = nmod_sub(initial[0], s, mod);
        else
            s = nmod_add(initial[0], s, mod);
    }

    *res = s;

    return GR_SUCCESS;
}

int
__gr_nmod_vec_dot_rev(ulong * res, const ulong * initial, int subtract, const ulong * vec1, const ulong * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    ulong s;
    int nlimbs;
    nmod_t mod;

    if (len <= 0)
    {
        if (initial == NULL)
            _gr_nmod_zero(res, ctx);
        else
            _gr_nmod_set(res, initial, ctx);
        return GR_SUCCESS;
    }

    mod = NMOD_CTX(ctx);

    nlimbs = _nmod_vec_dot_bound_limbs(len, mod);

    NMOD_VEC_DOT(s, i, len, vec1[i], vec2[len - 1 - i], mod, nlimbs);

    if (initial == NULL)
    {
        if (subtract)
            s = nmod_neg(s, mod);
    }
    else
    {
        if (subtract)
            s = nmod_sub(initial[0], s, mod);
        else
            s = nmod_add(initial[0], s, mod);
    }

    *res = s;

    return GR_SUCCESS;
}

int
_gr_nmod_roots_gr_poly(gr_vec_t roots, gr_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx)
{
    if (poly->length == 0)
        return GR_DOMAIN;

    {
        gr_poly_t z_poly;
        gr_vec_t z_roots;
        gr_ctx_t z_ctx;
        slong i;
        int status = GR_SUCCESS;

        fmpz_t t;
        fmpz_init(t);

        fmpz_set_ui(t, NMOD_CTX(ctx).n);

        gr_ctx_init_fmpz_mod(z_ctx, t);
        gr_poly_init(z_poly, z_ctx);
        gr_vec_init(z_roots, 0, z_ctx);

        status |= gr_poly_set_gr_poly_other(z_poly, poly, ctx, z_ctx);
        status |= gr_poly_roots(z_roots, mult, z_poly, flags, z_ctx);

        if (status == GR_SUCCESS)
        {
            gr_vec_set_length(roots, z_roots->length, ctx);
            for (i = 0; i < z_roots->length; i++)
                status |= gr_set_other(gr_vec_entry_ptr(roots, i, ctx), gr_vec_entry_ptr(z_roots, i, z_ctx), z_ctx, ctx);
        }

        gr_poly_clear(z_poly, z_ctx);
        gr_vec_clear(z_roots, z_ctx);
        gr_ctx_clear(z_ctx);
        fmpz_clear(t);

        return status;
    }

}

int __gr_nmod_methods_initialized = 0;

gr_static_method_table __gr_nmod_methods;

gr_method_tab_input __gr_nmod_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_nmod_ctx_write},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) _gr_nmod_ctx_is_field},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) _gr_nmod_ctx_is_field},
    {GR_METHOD_CTX_IS_FINITE,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_INIT,            (gr_funcptr) _gr_nmod_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_nmod_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_nmod_swap},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_nmod_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_nmod_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_nmod_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_nmod_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_nmod_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_nmod_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_nmod_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_nmod_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_nmod_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_nmod_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_nmod_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_nmod_set_fmpz},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_nmod_set_other},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_nmod_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_nmod_add},
    {GR_METHOD_ADD_SI,          (gr_funcptr) _gr_nmod_add_si},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_nmod_sub},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_nmod_mul},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_nmod_mul_si},
    {GR_METHOD_ADDMUL,          (gr_funcptr) _gr_nmod_addmul},
    {GR_METHOD_SUBMUL,          (gr_funcptr) _gr_nmod_submul},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) _gr_nmod_mul_two},
    {GR_METHOD_SQR,             (gr_funcptr) _gr_nmod_sqr},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_nmod_div},
    {GR_METHOD_DIV_SI,          (gr_funcptr) _gr_nmod_div_si},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_nmod_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) _gr_nmod_inv},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) __gr_nmod_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) __gr_nmod_vec_dot_rev},
    {GR_METHOD_POLY_ROOTS,      (gr_funcptr) _gr_nmod_roots_gr_poly},
    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_nmod(gr_ctx_t ctx, ulong n)
{
    ctx->which_ring = GR_CTX_NMOD;
    ctx->sizeof_elem = sizeof(ulong);
    ctx->size_limit = WORD_MAX;

    nmod_init(NMOD_CTX_REF(ctx), n);

    ctx->methods = __gr_nmod_methods;

    if (!__gr_nmod_methods_initialized)
    {
        gr_method_tab_init(__gr_nmod_methods, __gr_nmod_methods_input);
        __gr_nmod_methods_initialized = 1;
    }
}
