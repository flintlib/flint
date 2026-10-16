/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Quotient ring F[x]/(m) of a polynomial ring over an arbitrary
    commutative ring F by a monic polynomial m.

    Elements are represented as polynomials over F (gr_poly_struct)
    reduced modulo m, i.e. of degree < deg(m). Elements of larger
    degree may temporarily exist after the modulus has been refined
    (see below); all operations reduce such operands on the fly, so
    that stale elements remain valid representatives.

    Reduction, multiplication and powers use the precomputed modulus
    (gr_poly_preinv_t), whose representation is chosen by the base
    ring (GR_METHOD_POLY_PREINV_SET).

    Dynamic evaluation (D5): with gr_ctx_set_is_pretend_field, the ring
    pretends to be a field even if m is only conjecturally irreducible
    (gr_ctx_is_field remains T_UNKNOWN unless m has been declared
    irreducible with gr_ctx_set_is_field). An inversion which encounters
    a zero divisor then fails with GR_UNABLE, and the nontrivial factor
    of m that was found is recorded in the context's zero divisor table
    (see also gr_ctx_recover_zero_divisor). The caller may then replace
    the modulus by a factor (gr_poly_quotient_ctx_refine); since the
    new modulus divides the old one, existing elements remain valid
    representatives and need no conversion.
*/

#include <string.h>
#include "fmpq.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_mat.h"
#include "gr_generic.h"

#if FLINT_USES_PTHREAD
#include <pthread.h>
#endif

/*
    The modulus together with precomputed data for reduction modulo it
    (gr_poly_preinv_t: the representation is chosen by the base ring's
    GR_METHOD_POLY_PREINV_SET method). A refinement replaces the whole
    object; superseded objects are kept alive until the context is
    cleared, since concurrent operations may still be using them.
*/
typedef struct
{
    gr_poly_struct modulus;
    gr_poly_preinv_struct preinv;
}
_gr_poly_quotient_mod_struct;

typedef struct
{
    gr_ctx_struct * base;
    _gr_poly_quotient_mod_struct * mod;   /* current (monic) modulus */
    truth_t is_field;                 /* T_TRUE: m is known irreducible */
    truth_t is_pretend_field;         /* T_TRUE: pretend to be a field (dynamic evaluation) */
    char * var;
    ulong version;                    /* incremented on each refinement */
    gr_poly_struct * zero_divisors;   /* nontrivial factors of the modulus found by inv */
    slong num_zero_divisors;
    slong alloc_zero_divisors;
    _gr_poly_quotient_mod_struct ** old_moduli;   /* superseded moduli (kept alive until ctx_clear) */
    slong num_old_moduli;
    slong alloc_old_moduli;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
}
gr_poly_quotient_ctx_struct;

static const char * default_var = "x";

#define QCTX(ctx) ((gr_poly_quotient_ctx_struct *) GR_CTX_DATA_AS_PTR(ctx))
#define BASE(ctx) (QCTX(ctx)->base)
#define MODULUS(ctx) (&QCTX(ctx)->mod->modulus)
#define PREINV(ctx) (&QCTX(ctx)->mod->preinv)
#define DEG(ctx) (QCTX(ctx)->mod->modulus.length - 1)

static _gr_poly_quotient_mod_struct *
_mod_new(const gr_poly_t m, gr_ctx_t base)
{
    _gr_poly_quotient_mod_struct * M = flint_malloc(sizeof(_gr_poly_quotient_mod_struct));

    gr_poly_init(&M->modulus, base);
    GR_MUST_SUCCEED(gr_poly_set(&M->modulus, m, base));

    gr_poly_preinv_init(&M->preinv, base);
    if (gr_poly_preinv_set(&M->preinv, &M->modulus, base) != GR_SUCCESS)
    {
        /* (e.g. a Newton inverse which cannot be computed in an inexact
           ring): use ordinary division by the monic modulus */
        gr_poly_preinv_clear(&M->preinv, base);
        gr_poly_preinv_init(&M->preinv, base);
        GR_MUST_SUCCEED(_gr_poly_preinv_set_plain(&M->preinv, M->modulus.coeffs, M->modulus.length, base));
    }

    return M;
}

static void
_mod_clear(_gr_poly_quotient_mod_struct * M, gr_ctx_t base)
{
    gr_poly_preinv_clear(&M->preinv, base);
    gr_poly_clear(&M->modulus, base);
    flint_free(M);
}

static void
_qctx_lock(gr_ctx_t ctx)
{
#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&QCTX(ctx)->mutex);
#endif
}

static void
_qctx_unlock(gr_ctx_t ctx)
{
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&QCTX(ctx)->mutex);
#endif
}

/* -------------------------------------------------------------------- */
/* context                                                               */
/* -------------------------------------------------------------------- */

static int
_gr_poly_quotient_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    status |= gr_stream_write(out, "Quotient ring (");
    status |= gr_ctx_write(out, BASE(ctx));
    status |= gr_stream_write(out, ")[");
    status |= gr_stream_write(out, QCTX(ctx)->var);
    status |= gr_stream_write(out, "] / (");
    status |= gr_poly_write(out, MODULUS(ctx), QCTX(ctx)->var, BASE(ctx));
    status |= gr_stream_write(out, ")");
    return status;
}

static void
_gr_poly_quotient_ctx_clear(gr_ctx_t ctx)
{
    gr_poly_quotient_ctx_struct * q = QCTX(ctx);
    slong i;

    _mod_clear(q->mod, q->base);

    for (i = 0; i < q->num_old_moduli; i++)
        _mod_clear(q->old_moduli[i], q->base);
    flint_free(q->old_moduli);

    for (i = 0; i < q->num_zero_divisors; i++)
        gr_poly_clear(q->zero_divisors + i, q->base);
    flint_free(q->zero_divisors);

    if (q->var != default_var)
        flint_free(q->var);

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&q->mutex);
#endif

    flint_free(q);
}

static truth_t _gr_poly_quotient_ctx_is_ring(gr_ctx_t ctx) { return gr_ctx_is_commutative_ring(BASE(ctx)); }
static truth_t _gr_poly_quotient_ctx_is_commutative_ring(gr_ctx_t ctx) { return gr_ctx_is_commutative_ring(BASE(ctx)); }

static truth_t
_gr_poly_quotient_ctx_is_field(gr_ctx_t ctx)
{
    if (QCTX(ctx)->is_field == T_TRUE)
        return gr_ctx_is_field(BASE(ctx)) == T_TRUE ? T_TRUE : T_UNKNOWN;
    return QCTX(ctx)->is_field;
}

static truth_t
_gr_poly_quotient_ctx_is_integral_domain(gr_ctx_t ctx)
{
    return _gr_poly_quotient_ctx_is_field(ctx);
}

static int
_gr_poly_quotient_ctx_set_is_field(gr_ctx_t ctx, truth_t is_field)
{
    QCTX(ctx)->is_field = is_field;
    return GR_SUCCESS;
}

/* A field, or pretending to be one: the modulus is irreducible or
   assumed to be, over a base ring which is (or pretends to be) a field. */
static truth_t
_gr_poly_quotient_ctx_is_pretend_field(gr_ctx_t ctx)
{
    if (QCTX(ctx)->is_field != T_TRUE && QCTX(ctx)->is_pretend_field != T_TRUE)
        return T_FALSE;

    return (gr_ctx_is_pretend_field(BASE(ctx)) == T_TRUE) ? T_TRUE : T_FALSE;
}

static void _gr_poly_quotient_clear_zero_divisors(gr_ctx_t ctx);

static int
_gr_poly_quotient_ctx_set_is_pretend_field(gr_ctx_t ctx, truth_t flag)
{
    if (flag == T_TRUE)
    {
        QCTX(ctx)->is_pretend_field = T_TRUE;
    }
    else
    {
        QCTX(ctx)->is_pretend_field = T_FALSE;
        _gr_poly_quotient_clear_zero_divisors(ctx);
    }

    return GR_SUCCESS;
}

static truth_t
_gr_poly_quotient_ctx_is_rational_vector_space(gr_ctx_t ctx)
{
    return gr_ctx_is_rational_vector_space(BASE(ctx));
}

static truth_t
_gr_poly_quotient_ctx_is_real_vector_space(gr_ctx_t ctx)
{
    return gr_ctx_is_real_vector_space(BASE(ctx));
}

static truth_t
_gr_poly_quotient_ctx_is_complex_vector_space(gr_ctx_t ctx)
{
    return gr_ctx_is_complex_vector_space(BASE(ctx));
}

static truth_t _gr_poly_quotient_ctx_is_finite(gr_ctx_t ctx) { return gr_ctx_is_finite(BASE(ctx)); }
static truth_t _gr_poly_quotient_ctx_is_finite_characteristic(gr_ctx_t ctx) { return gr_ctx_is_finite_characteristic(BASE(ctx)); }
static truth_t _gr_poly_quotient_ctx_is_exact(gr_ctx_t ctx) { return gr_ctx_is_exact(BASE(ctx)); }
static truth_t _gr_poly_quotient_ctx_is_canonical(gr_ctx_t ctx) { return gr_ctx_is_canonical(BASE(ctx)); }

static truth_t
_gr_poly_quotient_ctx_is_threadsafe(gr_ctx_t ctx)
{
#if FLINT_USES_PTHREAD
    return gr_ctx_is_threadsafe(BASE(ctx));
#else
    return T_FALSE;
#endif
}

static gr_ptr _gr_poly_quotient_ctx_base(gr_ctx_t ctx) { return BASE(ctx); }

static int
_gr_poly_quotient_ctx_set_gen_name(gr_ctx_t ctx, const char * s)
{
    slong len = strlen(s);

    if (QCTX(ctx)->var == default_var)
        QCTX(ctx)->var = NULL;

    QCTX(ctx)->var = flint_realloc(QCTX(ctx)->var, len + 1);
    memcpy(QCTX(ctx)->var, s, len + 1);
    return GR_SUCCESS;
}

static int
_gr_poly_quotient_ctx_set_gen_names(gr_ctx_t ctx, const char ** s)
{
    return _gr_poly_quotient_ctx_set_gen_name(ctx, s[0]);
}

static int
_gr_poly_quotient_ctx_gen_name(char ** name, slong i, gr_ctx_t ctx)
{
    size_t len;

    if (i != 0)
        return GR_DOMAIN;

    len = strlen(QCTX(ctx)->var);
    *name = flint_malloc(len + 1);
    memcpy(*name, QCTX(ctx)->var, len + 1);
    return GR_SUCCESS;
}

/* -------------------------------------------------------------------- */
/* reduction                                                             */
/* -------------------------------------------------------------------- */

/*
    A polynomial may be in a non-normalized state if a refinement of the
    base ring (when the base ring is itself a quotient ring) has turned
    its leading coefficients into (non-canonical) zeros. Such leading
    coefficients are stripped before any operation which depends on the
    length, using the base ring's zero test.
*/
static int
_is_normalized(const gr_poly_struct * x, gr_ctx_t ctx)
{
    if (x->length == 0)
        return 1;

    /* only a quotient ring base can produce stale zeros */
    if (BASE(ctx)->which_ring != GR_CTX_GR_POLY_QUOTIENT)
        return 1;

    return gr_is_zero(GR_ENTRY(x->coeffs, x->length - 1, BASE(ctx)->sizeof_elem), BASE(ctx)) != T_TRUE;
}

/* Normalize and reduce in place modulo the current modulus. */
static int
_reduce(gr_poly_t x, gr_ctx_t ctx)
{
    const _gr_poly_quotient_mod_struct * M = QCTX(ctx)->mod;

    if (!_is_normalized(x, ctx))
        _gr_poly_normalise(x, BASE(ctx));

    if (x->length < M->modulus.length)
        return GR_SUCCESS;

    return gr_poly_preinv_rem(x, x, &M->preinv, BASE(ctx));
}

/* Return a pointer to a normalized, reduced version of x, using tmp as
   scratch if x is not already in that state. The caller must clear tmp. */
static const gr_poly_struct *
_reduced(const gr_poly_struct * x, gr_poly_t tmp, int * status, gr_ctx_t ctx)
{
    if (x->length < MODULUS(ctx)->length && _is_normalized(x, ctx))
        return x;

    *status |= gr_poly_set(tmp, x, BASE(ctx));
    *status |= _reduce(tmp, ctx);
    return tmp;
}

/* -------------------------------------------------------------------- */
/* elements                                                              */
/* -------------------------------------------------------------------- */

static void
_gr_poly_quotient_init(gr_poly_t x, gr_ctx_t ctx)
{
    gr_poly_init(x, BASE(ctx));
}

static void
_gr_poly_quotient_clear(gr_poly_t x, gr_ctx_t ctx)
{
    gr_poly_clear(x, BASE(ctx));
}

static void
_gr_poly_quotient_swap(gr_poly_t x, gr_poly_t y, gr_ctx_t ctx)
{
    gr_poly_swap(x, y, BASE(ctx));
}

static void
_gr_poly_quotient_set_shallow(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
{
    *res = *x;
}

static int
_gr_poly_quotient_randtest(gr_poly_t res, flint_rand_t state, gr_ctx_t ctx)
{
    int status;
    status = gr_poly_randtest(res, state, DEG(ctx), BASE(ctx));
    status |= _reduce(res, ctx);
    return status;
}

static int
_gr_poly_quotient_write(gr_stream_t out, const gr_poly_t x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_poly_t t;
    const gr_poly_struct * xr;

    gr_poly_init(t, BASE(ctx));
    xr = _reduced(x, t, &status, ctx);
    status |= gr_poly_write(out, xr, QCTX(ctx)->var, BASE(ctx));
    gr_poly_clear(t, BASE(ctx));
    return status;
}

static int
_gr_poly_quotient_zero(gr_poly_t res, gr_ctx_t ctx)
{
    return gr_poly_zero(res, BASE(ctx));
}

static int
_gr_poly_quotient_one(gr_poly_t res, gr_ctx_t ctx)
{
    int status = gr_poly_one(res, BASE(ctx));
    return status | _reduce(res, ctx);
}

static int
_gr_poly_quotient_gen(gr_poly_t res, gr_ctx_t ctx)
{
    int status = gr_poly_gen(res, BASE(ctx));
    return status | _reduce(res, ctx);
}

static int
_gr_poly_quotient_gens_recursive(gr_vec_t vec, gr_ctx_t ctx)
{
    int status;
    gr_vec_t vec1;
    slong i, n;

    gr_vec_init(vec1, 0, BASE(ctx));
    status = gr_gens_recursive(vec1, BASE(ctx));
    n = vec1->length;

    gr_vec_set_length(vec, n + 1, ctx);

    for (i = 0; i < n; i++)
        status |= gr_poly_set_scalar(gr_vec_entry_ptr(vec, i, ctx),
                gr_vec_entry_srcptr(vec1, i, BASE(ctx)), BASE(ctx));

    status |= _gr_poly_quotient_gen(gr_vec_entry_ptr(vec, n, ctx), ctx);

    gr_vec_clear(vec1, BASE(ctx));
    return status;
}

static int
_gr_poly_quotient_set(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
{
    int status = gr_poly_set(res, x, BASE(ctx));
    return status | _reduce(res, ctx);
}

static int
_gr_poly_quotient_set_si(gr_poly_t res, slong c, gr_ctx_t ctx)
{
    int status = gr_poly_set_si(res, c, BASE(ctx));
    return status | _reduce(res, ctx);
}

static int
_gr_poly_quotient_set_ui(gr_poly_t res, ulong c, gr_ctx_t ctx)
{
    int status = gr_poly_set_ui(res, c, BASE(ctx));
    return status | _reduce(res, ctx);
}

static int
_gr_poly_quotient_set_fmpz(gr_poly_t res, const fmpz_t c, gr_ctx_t ctx)
{
    int status = gr_poly_set_fmpz(res, c, BASE(ctx));
    return status | _reduce(res, ctx);
}

static int
_gr_poly_quotient_set_fmpq(gr_poly_t res, const fmpq_t c, gr_ctx_t ctx)
{
    int status = gr_poly_set_fmpq(res, c, BASE(ctx));
    return status | _reduce(res, ctx);
}

static int
_gr_poly_quotient_set_other(gr_poly_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    if (x_ctx == ctx)
    {
        return _gr_poly_quotient_set(res, x, ctx);
    }
    else if (x_ctx == BASE(ctx))
    {
        int status = gr_poly_set_scalar(res, x, BASE(ctx));
        return status | _reduce(res, ctx);
    }
    else if (x_ctx->which_ring == GR_CTX_GR_POLY && POLYNOMIAL_ELEM_CTX(x_ctx) == BASE(ctx))
    {
        int status = gr_poly_set(res, x, BASE(ctx));
        return status | _reduce(res, ctx);
    }
    else if (x_ctx->which_ring == GR_CTX_GR_POLY_QUOTIENT && BASE(x_ctx) == BASE(ctx))
    {
        /* Same base ring: allowed only if the moduli agree
           (e.g. two contexts representing the same ring), since
           a general map between quotient rings is not canonical. */
        if (gr_poly_equal(MODULUS(x_ctx), MODULUS(ctx), BASE(ctx)) == T_TRUE)
        {
            int status = gr_poly_set(res, x, BASE(ctx));
            return status | _reduce(res, ctx);
        }
        return GR_DOMAIN;
    }
    else
    {
        int status;

        gr_poly_fit_length(res, 1, BASE(ctx));
        status = gr_set_other(res->coeffs, x, x_ctx, BASE(ctx));

        if (status == GR_SUCCESS)
        {
            _gr_poly_set_length(res, 1, BASE(ctx));
            _gr_poly_normalise(res, BASE(ctx));
            status |= _reduce(res, ctx);
        }
        else
        {
            _gr_poly_set_length(res, 0, BASE(ctx));
        }

        return status;
    }
}

static truth_t
_gr_poly_quotient_is_zero(const gr_poly_t x, gr_ctx_t ctx)
{
    truth_t res;
    int status = GR_SUCCESS;
    gr_poly_t t;
    const gr_poly_struct * xr;

    if (x->length == 0)
        return T_TRUE;

    /* an unreduced element is zero iff all its coefficients are zero
       (stale leading zeros included), so normalization is not needed;
       avoiding it keeps the cost linear in the depth of a nested
       quotient (normalization and the coefficient test would otherwise
       both recurse) */
    if (x->length < MODULUS(ctx)->length)
        return _gr_vec_is_zero(x->coeffs, x->length, BASE(ctx));

    gr_poly_init(t, BASE(ctx));
    xr = _reduced(x, t, &status, ctx);
    res = (status == GR_SUCCESS) ? gr_poly_is_zero(xr, BASE(ctx)) : T_UNKNOWN;
    gr_poly_clear(t, BASE(ctx));
    return res;
}

static truth_t
_gr_poly_quotient_is_one(const gr_poly_t x, gr_ctx_t ctx)
{
    truth_t res;
    int status = GR_SUCCESS;
    gr_poly_t t, u;
    const gr_poly_struct * xr;

    gr_poly_init(t, BASE(ctx));
    gr_poly_init(u, BASE(ctx));
    xr = _reduced(x, t, &status, ctx);
    status |= _gr_poly_quotient_one(u, ctx);
    res = (status == GR_SUCCESS) ? gr_poly_equal(xr, u, BASE(ctx)) : T_UNKNOWN;
    gr_poly_clear(t, BASE(ctx));
    gr_poly_clear(u, BASE(ctx));
    return res;
}

static truth_t
_gr_poly_quotient_equal(const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx)
{
    truth_t res;
    int status = GR_SUCCESS;
    gr_poly_t t, u;
    const gr_poly_struct * xr, * yr;

    gr_poly_init(t, BASE(ctx));
    gr_poly_init(u, BASE(ctx));
    xr = _reduced(x, t, &status, ctx);
    yr = _reduced(y, u, &status, ctx);
    res = (status == GR_SUCCESS) ? gr_poly_equal(xr, yr, BASE(ctx)) : T_UNKNOWN;
    gr_poly_clear(t, BASE(ctx));
    gr_poly_clear(u, BASE(ctx));
    return res;
}

static int
_gr_poly_quotient_neg(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
{
    int status = gr_poly_neg(res, x, BASE(ctx));
    return status | _reduce(res, ctx);
}

static int
_gr_poly_quotient_add(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx)
{
    int status = gr_poly_add(res, x, y, BASE(ctx));
    return status | _reduce(res, ctx);
}

static int
_gr_poly_quotient_sub(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx)
{
    int status = gr_poly_sub(res, x, y, BASE(ctx));
    return status | _reduce(res, ctx);
}

/* Whether y is the generator (the residue class of x). */
static int
_is_gen(const gr_poly_t y, gr_ctx_t ctx)
{
    return y->length == 2 &&
        gr_is_one(gr_poly_coeff_srcptr(y, 1, BASE(ctx)), BASE(ctx)) == T_TRUE &&
        gr_is_zero(gr_poly_coeff_srcptr(y, 0, BASE(ctx)), BASE(ctx)) == T_TRUE;
}

static int
_gr_poly_quotient_mul(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx)
{
    int status;

    /* multiplication by the generator is a shift (a common operation in
       towers: evaluating polynomials at a generator, synthetic division),
       with at most one reduction step */
    if (_is_gen(y, ctx))
        status = gr_poly_shift_left(res, x, 1, BASE(ctx));
    else if (_is_gen(x, ctx))
        status = gr_poly_shift_left(res, y, 1, BASE(ctx));
    else
        return gr_poly_preinv_mulmod(res, x, y, PREINV(ctx), BASE(ctx));

    return status | _reduce(res, ctx);
}

static int
_gr_poly_quotient_sqr(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
{
    return gr_poly_preinv_mulmod(res, x, x, PREINV(ctx), BASE(ctx));
}

static int _gr_poly_quotient_inv(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx);

/* Powers by modular exponentiation with the precomputed modulus. */
static int
_gr_poly_quotient_pow_fmpz(gr_poly_t res, const gr_poly_t x, const fmpz_t e, gr_ctx_t ctx)
{
    gr_poly_t t;
    int status;

    gr_poly_init(t, BASE(ctx));

    if (fmpz_sgn(e) < 0)
    {
        fmpz_t f;
        fmpz_init(f);
        fmpz_neg(f, e);
        status = _gr_poly_quotient_inv(t, x, ctx);
        if (status == GR_SUCCESS)
            status = gr_poly_preinv_powmod_fmpz_sliding(res, t, f, 0, PREINV(ctx), BASE(ctx));
        fmpz_clear(f);
    }
    else
    {
        status = gr_poly_preinv_powmod_fmpz_sliding(t, x, e, 0, PREINV(ctx), BASE(ctx));
        gr_poly_swap(res, t, BASE(ctx));
    }

    gr_poly_clear(t, BASE(ctx));
    return status;
}

static int
_gr_poly_quotient_pow_ui(gr_poly_t res, const gr_poly_t x, ulong e, gr_ctx_t ctx)
{
    fmpz_t f;
    int status;
    fmpz_init_set_ui(f, e);
    status = _gr_poly_quotient_pow_fmpz(res, x, f, ctx);
    fmpz_clear(f);
    return status;
}

static int
_gr_poly_quotient_pow_si(gr_poly_t res, const gr_poly_t x, slong e, gr_ctx_t ctx)
{
    fmpz_t f;
    int status;
    fmpz_init(f);
    fmpz_set_si(f, e);
    status = _gr_poly_quotient_pow_fmpz(res, x, f, ctx);
    fmpz_clear(f);
    return status;
}

static int
_gr_poly_quotient_mul_other(gr_poly_t res, const gr_poly_t x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
{
    if (y_ctx == BASE(ctx))
    {
        int status = gr_poly_mul_scalar(res, x, y, BASE(ctx));
        return status | _reduce(res, ctx);
    }

    return gr_generic_mul_other(res, x, y, y_ctx, ctx);
}

/* Scalar operations act coefficientwise: through the base ring (which
   may itself be a quotient ring, recursively), never through an
   inversion in the quotient ring. */

#define SCALAR_MUL_OP(name, op, T) \
static int \
name(gr_poly_t res, const gr_poly_t x, T c, gr_ctx_t ctx) \
{ \
    int status = op(res, x, c, BASE(ctx)); \
    return status | _reduce(res, ctx); \
}

SCALAR_MUL_OP(_gr_poly_quotient_mul_ui, gr_poly_mul_ui, ulong)
SCALAR_MUL_OP(_gr_poly_quotient_mul_si, gr_poly_mul_si, slong)
SCALAR_MUL_OP(_gr_poly_quotient_mul_fmpz, gr_poly_mul_fmpz, const fmpz_t)
SCALAR_MUL_OP(_gr_poly_quotient_mul_fmpq, gr_poly_mul_fmpq, const fmpq_t)

#define SCALAR_DIV_OP(name, base_op, T) \
static int \
name(gr_poly_t res, const gr_poly_t x, T c, gr_ctx_t ctx) \
{ \
    slong i, len = x->length; \
    int status = GR_SUCCESS; \
    if (len == 0) \
    { \
        /* (the status for a zero dividend: division by zero) */ \
        gr_ptr t; \
        GR_TMP_INIT(t, BASE(ctx)); \
        status = base_op(t, t, c, BASE(ctx)); \
        GR_TMP_CLEAR(t, BASE(ctx)); \
        return status | gr_poly_zero(res, BASE(ctx)); \
    } \
    gr_poly_fit_length(res, len, BASE(ctx)); \
    for (i = 0; i < len; i++) \
        status |= base_op(gr_poly_coeff_ptr(res, i, BASE(ctx)), gr_poly_coeff_srcptr(x, i, BASE(ctx)), c, BASE(ctx)); \
    _gr_poly_set_length(res, len, BASE(ctx)); \
    _gr_poly_normalise(res, BASE(ctx)); \
    return status; \
}

SCALAR_DIV_OP(_gr_poly_quotient_div_ui, gr_div_ui, ulong)
SCALAR_DIV_OP(_gr_poly_quotient_div_si, gr_div_si, slong)
SCALAR_DIV_OP(_gr_poly_quotient_div_fmpz, gr_div_fmpz, const fmpz_t)
SCALAR_DIV_OP(_gr_poly_quotient_div_fmpq, gr_div_fmpq, const fmpq_t)

/* Record a nontrivial factor of the modulus found during an inversion. */
static void
_record_zero_divisor(const gr_poly_t g, gr_ctx_t ctx)
{
    gr_poly_quotient_ctx_struct * q = QCTX(ctx);
    slong i;

    _qctx_lock(ctx);

    /* avoid duplicates */
    for (i = 0; i < q->num_zero_divisors; i++)
    {
        if (gr_poly_equal(q->zero_divisors + i, g, q->base) == T_TRUE)
        {
            _qctx_unlock(ctx);
            return;
        }
    }

    if (q->num_zero_divisors == q->alloc_zero_divisors)
    {
        slong new_alloc = FLINT_MAX(4, 2 * q->alloc_zero_divisors);
        q->zero_divisors = flint_realloc(q->zero_divisors, new_alloc * sizeof(gr_poly_struct));
        q->alloc_zero_divisors = new_alloc;
    }

    gr_poly_init(q->zero_divisors + q->num_zero_divisors, q->base);
    GR_MUST_SUCCEED(gr_poly_set(q->zero_divisors + q->num_zero_divisors, g, q->base));
    q->num_zero_divisors++;

    _qctx_unlock(ctx);
}

/*
    Inversion using a single inversion in the base ring, via the
    characteristic polynomial of the multiplication-by-x matrix over the
    base ring (Cayley-Hamilton): if chi(t) = t^d + c_{d-1} t^{d-1} + ... + c_0,
    then x^{-1} = -(x^{d-1} + c_{d-1} x^{d-2} + ... + c_1) / c_0.

    This is preferable to the Euclidean algorithm when inversion in the
    base ring is expensive relative to multiplication (e.g. when the base
    ring is itself a quotient ring), since the Euclidean algorithm performs
    about d base ring inversions.

    Sets *invertible = 1 on success. If c_0 = 0 (x is zero or a zero
    divisor), sets *invertible = 0 without recording anything; the caller
    should fall back to the gcd computation to obtain a factor.
*/
static int
_gr_poly_quotient_inv_norm(gr_poly_t res, int * invertible, const gr_poly_t xr, gr_ctx_t ctx)
{
    gr_ctx_struct * F = BASE(ctx);
    slong d = DEG(ctx);
    slong sz = F->sizeof_elem;
    int status = GR_SUCCESS;

    *invertible = 0;

    if (d == 2)
    {
        /* m = t^2 + b t + c;  (u + v t)^-1 = (u - v b - v t) / (u^2 - b u v + c v^2) */
        const gr_poly_struct * m = MODULUS(ctx);
        gr_srcptr b = gr_poly_coeff_srcptr(m, 1, F);
        gr_srcptr c = gr_poly_coeff_srcptr(m, 0, F);
        gr_srcptr u = xr->coeffs;
        gr_srcptr v = (xr->length == 2) ? GR_ENTRY(xr->coeffs, 1, sz) : NULL;
        gr_ptr N, t;

        if (v == NULL)
        {
            /* constant */
            GR_TMP_INIT(t, F);
            status = gr_inv(t, u, F);
            if (status == GR_SUCCESS)
            {
                *invertible = 1;
                status |= gr_poly_set_scalar(res, t, F);
            }
            else if (status == GR_DOMAIN)
                status = GR_SUCCESS;
            /* (GR_UNABLE: the base ring met a zero divisor and has
               recorded it; this propagates) */
            GR_TMP_CLEAR(t, F);
            return status;
        }

        GR_TMP_INIT2(N, t, F);

        status |= gr_mul(N, u, u, F);
        status |= gr_mul(t, u, v, F);
        status |= gr_mul(t, t, b, F);
        status |= gr_sub(N, N, t, F);
        status |= gr_mul(t, v, v, F);
        status |= gr_mul(t, t, c, F);
        status |= gr_add(N, N, t, F);

        if (status == GR_SUCCESS)
        {
            status = gr_inv(N, N, F);

            if (status == GR_SUCCESS)
            {
                *invertible = 1;
                gr_poly_fit_length(res, 2, F);
                status |= gr_mul(t, v, b, F);
                status |= gr_sub(res->coeffs, u, t, F);
                status |= gr_mul(res->coeffs, res->coeffs, N, F);
                status |= gr_mul(GR_ENTRY(res->coeffs, 1, sz), v, N, F);
                status |= gr_neg(GR_ENTRY(res->coeffs, 1, sz), GR_ENTRY(res->coeffs, 1, sz), F);
                _gr_poly_set_length(res, 2, F);
                _gr_poly_normalise(res, F);
            }
            else if (status == GR_DOMAIN)
            {
                /* N is zero (or, in a base ring which does not pretend
                   to be a field, a non-unit). Report failure to the
                   caller, which falls back to the gcd. GR_UNABLE (a zero
                   divisor recorded by a base ring pretending to be a
                   field) propagates. */
                status = GR_SUCCESS;
            }
        }

        GR_TMP_CLEAR2(N, t, F);
        return status;
    }
    else
    {
        gr_mat_t M;
        gr_poly_t chi, xpow, t;
        slong i, j;
        gr_ptr c0;

        gr_mat_init(M, d, d, F);
        gr_poly_init(chi, F);
        gr_poly_init(xpow, F);
        gr_poly_init(t, F);
        GR_TMP_INIT(c0, F);

        /* column j of M = coefficients of x * t^j mod m */
        status |= gr_poly_set(xpow, xr, F);
        for (j = 0; j < d && status == GR_SUCCESS; j++)
        {
            for (i = 0; i < d; i++)
            {
                if (i < xpow->length)
                    status |= gr_set(gr_mat_entry_ptr(M, i, j, F), GR_ENTRY(xpow->coeffs, i, sz), F);
                else
                    status |= gr_zero(gr_mat_entry_ptr(M, i, j, F), F);
            }

            if (j + 1 < d)
            {
                status |= gr_poly_shift_left(xpow, xpow, 1, F);
                status |= _reduce(xpow, ctx);
            }
        }

        if (status == GR_SUCCESS)
            status = gr_mat_charpoly_berkowitz(chi, M, F);

        if (status == GR_SUCCESS)
        {
            /* c0 = chi(0) */
            status |= gr_poly_get_coeff_scalar(c0, chi, 0, F);
            status = gr_inv(c0, c0, F);

            if (status == GR_SUCCESS)
            {
                /* res = -(chi(x) - c0) / (x c0) = -(x^{d-1} + c_{d-1} x^{d-2} + ... + c_1) / c0 */
                status |= gr_poly_shift_right(t, chi, 1, F);       /* (chi - c0) / t */
                /* evaluate t at x via Horner in the quotient ring */
                status |= gr_poly_zero(res, F);
                for (i = t->length - 1; i >= 0 && status == GR_SUCCESS; i--)
                {
                    status |= gr_poly_preinv_mulmod(res, res, xr, PREINV(ctx), F);
                    status |= gr_poly_add_scalar(res, res, GR_ENTRY(t->coeffs, i, sz), F);
                }
                status |= gr_poly_mul_scalar(res, res, c0, F);
                status |= gr_poly_neg(res, res, F);
                *invertible = 1;
            }
            else if (status == GR_DOMAIN)
            {
                status = GR_SUCCESS;
            }
        }

        gr_mat_clear(M, F);
        gr_poly_clear(chi, F);
        gr_poly_clear(xpow, F);
        gr_poly_clear(t, F);
        GR_TMP_CLEAR(c0, F);
        return status;
    }
}

/*
    Extended gcd of x with the modulus. Sets *invertible to 1 and res
    to the inverse if x is a unit, and to 0 otherwise. If the ring
    pretends to be a field and x is a zero divisor, the nontrivial
    factor of the modulus that was found is recorded and GR_UNABLE is
    returned. Otherwise, a nonzero status indicates failure in the base
    ring (GR_UNABLE if a base ring pretending to be a field met a zero
    divisor, which that ring has recorded).
*/
static int _gr_poly_quotient_try_inv_noalias(gr_poly_t res, int * invertible, const gr_poly_t x, gr_ctx_t ctx);

static int
_gr_poly_quotient_try_inv(gr_poly_t res, int * invertible, const gr_poly_t x, gr_ctx_t ctx)
{
    /* the algorithms below read x while writing res */
    if (res == x)
    {
        gr_poly_t t;
        int status;
        gr_poly_init(t, BASE(ctx));
        status = _gr_poly_quotient_try_inv_noalias(t, invertible, x, ctx);
        if (*invertible)
            gr_poly_swap(res, t, BASE(ctx));
        gr_poly_clear(t, BASE(ctx));
        return status;
    }

    return _gr_poly_quotient_try_inv_noalias(res, invertible, x, ctx);
}

static int
_gr_poly_quotient_try_inv_noalias(gr_poly_t res, int * invertible, const gr_poly_t x, gr_ctx_t ctx)
{
    gr_poly_t g, s, t, xr_tmp;
    const gr_poly_struct * xr;
    const gr_poly_struct * m = MODULUS(ctx);
    int status = GR_SUCCESS;

    *invertible = 0;

    gr_poly_init(xr_tmp, BASE(ctx));
    xr = _reduced(x, xr_tmp, &status, ctx);

    if (status != GR_SUCCESS || xr->length == 0)
    {
        gr_poly_clear(xr_tmp, BASE(ctx));
        return status;
    }

    /* a constant is inverted in the base ring */
    if (xr->length == 1)
    {
        gr_ptr c;
        GR_TMP_INIT(c, BASE(ctx));
        status = gr_inv(c, gr_poly_coeff_srcptr(xr, 0, BASE(ctx)), BASE(ctx));
        if (status == GR_SUCCESS)
        {
            status = gr_poly_set_scalar(res, c, BASE(ctx));
            *invertible = 1;
        }
        else if (status == GR_DOMAIN)
            status = GR_SUCCESS;   /* a non-unit of the base ring */
        /* (GR_UNABLE from a base ring pretending to be a field propagates:
           the base ring has recorded a zero divisor) */
        GR_TMP_CLEAR(c, BASE(ctx));
        gr_poly_clear(xr_tmp, BASE(ctx));
        return status;
    }

    /* When the base ring is itself a quotient ring, base ring inversions
       are expensive: use a method with a single base ring inversion. */
    if (BASE(ctx)->which_ring == GR_CTX_GR_POLY_QUOTIENT && DEG(ctx) <= 8)
    {
        status = _gr_poly_quotient_inv_norm(res, invertible, xr, ctx);

        if (status != GR_SUCCESS || *invertible)
        {
            gr_poly_clear(xr_tmp, BASE(ctx));
            return status;
        }

        /* not invertible: fall through to the gcd to find a factor */
    }

    gr_poly_init(g, BASE(ctx));
    gr_poly_init(s, BASE(ctx));
    gr_poly_init(t, BASE(ctx));

    /* Over a field, the Euclidean algorithm is always applicable
       (and is the D5 path: a failure to invert a leading coefficient
       in a base ring that pretends to be a field propagates as
       GR_UNABLE, with that ring having recorded a zero divisor). */
    if (gr_ctx_is_pretend_field(BASE(ctx)) == T_TRUE)
        status = gr_poly_xgcd_euclidean(g, s, t, xr, m, BASE(ctx));
    else
        status = gr_poly_xgcd(g, s, t, xr, m, BASE(ctx));

    if (status == GR_SUCCESS)
    {
        if (g->length == 1)
        {
            /* s x + t m = g with g constant: x is invertible iff g is a unit */
            if (gr_is_one(g->coeffs, BASE(ctx)) == T_TRUE)
            {
                *invertible = 1;
                status |= gr_poly_set(res, s, BASE(ctx));
                status |= _reduce(res, ctx);
            }
            else
            {
                gr_ptr c;
                int st;

                GR_TMP_INIT(c, BASE(ctx));
                st = gr_inv(c, g->coeffs, BASE(ctx));

                if (st == GR_SUCCESS)
                {
                    *invertible = 1;
                    status |= gr_poly_mul_scalar(res, s, c, BASE(ctx));
                    status |= _reduce(res, ctx);
                }
                else if (st != GR_DOMAIN)
                {
                    status |= st;
                }

                GR_TMP_CLEAR(c, BASE(ctx));
            }
        }
        else if (g->length > 1 && g->length < m->length)
        {
            /* x is a zero divisor: a contradiction to the pretense,
               which is reported; in a genuine quotient ring, x is
               simply not invertible */
            if (_gr_poly_quotient_ctx_is_pretend_field(ctx) == T_TRUE)
            {
                _record_zero_divisor(g, ctx);
                status = GR_UNABLE;
            }
        }
        /* g->length == m->length would mean x = 0 mod m, excluded above
           unless the base ring is inexact/non-canonical */
    }

    gr_poly_clear(g, BASE(ctx));
    gr_poly_clear(s, BASE(ctx));
    gr_poly_clear(t, BASE(ctx));
    gr_poly_clear(xr_tmp, BASE(ctx));

    return status;
}

static int
_gr_poly_quotient_inv(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
{
    int status, invertible;

    if (DEG(ctx) == 0)   /* zero ring */
        return GR_DOMAIN;

    status = _gr_poly_quotient_try_inv(res, &invertible, x, ctx);

    if (status != GR_SUCCESS)
        return status;

    return invertible ? GR_SUCCESS : GR_DOMAIN;
}

static truth_t
_gr_poly_quotient_is_invertible(const gr_poly_t x, gr_ctx_t ctx)
{
    int status, invertible;
    gr_poly_t t;
    truth_t res;

    if (DEG(ctx) == 0)
        return T_FALSE;

    gr_poly_init(t, BASE(ctx));
    status = _gr_poly_quotient_try_inv(t, &invertible, x, ctx);
    res = (status != GR_SUCCESS) ? T_UNKNOWN : (invertible ? T_TRUE : T_FALSE);
    gr_poly_clear(t, BASE(ctx));
    return res;
}

static int
_gr_poly_quotient_div(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx)
{
    gr_poly_t t;
    int status;

    gr_poly_init(t, BASE(ctx));
    status = _gr_poly_quotient_inv(t, y, ctx);
    if (status == GR_SUCCESS)
        status = _gr_poly_quotient_mul(res, x, t, ctx);
    gr_poly_clear(t, BASE(ctx));
    return status;
}

/* -------------------------------------------------------------------- */
/* public helpers                                                        */
/* -------------------------------------------------------------------- */

const gr_poly_struct *
gr_poly_quotient_ctx_modulus(gr_ctx_t ctx)
{
    return MODULUS(ctx);
}

gr_ctx_struct *
gr_poly_quotient_ctx_base(gr_ctx_t ctx)
{
    return BASE(ctx);
}

slong
gr_poly_quotient_ctx_degree(gr_ctx_t ctx)
{
    return DEG(ctx);
}

slong
gr_poly_quotient_ctx_num_zero_divisors(gr_ctx_t ctx)
{
    return QCTX(ctx)->num_zero_divisors;
}

const gr_poly_struct *
gr_poly_quotient_ctx_zero_divisor(gr_ctx_t ctx, slong i)
{
    return QCTX(ctx)->zero_divisors + i;
}

static void
_gr_poly_quotient_clear_zero_divisors(gr_ctx_t ctx)
{
    gr_poly_quotient_ctx_clear_zero_divisors(ctx);
}

static int
_gr_poly_quotient_ctx_recover_zero_divisor(gr_poly_t res, gr_ctx_t ctx)
{
    gr_poly_quotient_ctx_struct * q = QCTX(ctx);
    int status = GR_UNABLE;

    _qctx_lock(ctx);
    if (q->num_zero_divisors > 0)
        status = gr_poly_set(res, q->zero_divisors, q->base);
    _qctx_unlock(ctx);

    return status;
}

void
gr_poly_quotient_ctx_clear_zero_divisors(gr_ctx_t ctx)
{
    gr_poly_quotient_ctx_struct * q = QCTX(ctx);
    slong i;

    _qctx_lock(ctx);
    for (i = 0; i < q->num_zero_divisors; i++)
        gr_poly_clear(q->zero_divisors + i, q->base);
    q->num_zero_divisors = 0;
    _qctx_unlock(ctx);
}

int
gr_poly_quotient_ctx_refine(gr_ctx_t ctx, const gr_poly_t new_modulus)
{
    gr_poly_quotient_ctx_struct * q = QCTX(ctx);
    _gr_poly_quotient_mod_struct * m;

    if (new_modulus->length < 2)
        return GR_DOMAIN;

    if (gr_is_one(GR_ENTRY(new_modulus->coeffs, new_modulus->length - 1, q->base->sizeof_elem), q->base) != T_TRUE)
        return GR_DOMAIN;

    m = _mod_new(new_modulus, q->base);

    _qctx_lock(ctx);

    /* keep the old modulus alive: concurrent readers may hold a pointer */
    if (q->num_old_moduli == q->alloc_old_moduli)
    {
        slong new_alloc = FLINT_MAX(4, 2 * q->alloc_old_moduli);
        q->old_moduli = flint_realloc(q->old_moduli, new_alloc * sizeof(_gr_poly_quotient_mod_struct *));
        q->alloc_old_moduli = new_alloc;
    }
    q->old_moduli[q->num_old_moduli++] = q->mod;
    q->mod = m;
    q->version++;

    _qctx_unlock(ctx);

    return GR_SUCCESS;
}

ulong
gr_poly_quotient_ctx_version(gr_ctx_t ctx)
{
    return QCTX(ctx)->version;
}

int
gr_poly_quotient_get_poly(gr_poly_t res, gr_srcptr x, gr_ctx_t ctx)
{
    int status = gr_poly_set(res, x, BASE(ctx));
    return status | _reduce(res, ctx);
}

int
gr_poly_quotient_set_poly(gr_ptr res, const gr_poly_t x, gr_ctx_t ctx)
{
    int status = gr_poly_set(res, x, BASE(ctx));
    return status | _reduce(res, ctx);
}

/* -------------------------------------------------------------------- */
/* method table and init                                                 */
/* -------------------------------------------------------------------- */

int _gr_poly_quotient_methods_initialized = 0;

gr_static_method_table _gr_poly_quotient_methods;

gr_method_tab_input _gr_poly_quotient_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_poly_quotient_ctx_write},
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_poly_quotient_ctx_clear},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) _gr_poly_quotient_ctx_is_ring},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) _gr_poly_quotient_ctx_is_commutative_ring},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) _gr_poly_quotient_ctx_is_integral_domain},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) _gr_poly_quotient_ctx_is_field},
    {GR_METHOD_CTX_SET_IS_FIELD,        (gr_funcptr) _gr_poly_quotient_ctx_set_is_field},
    {GR_METHOD_CTX_IS_PRETEND_FIELD,    (gr_funcptr) _gr_poly_quotient_ctx_is_pretend_field},
    {GR_METHOD_CTX_SET_IS_PRETEND_FIELD, (gr_funcptr) _gr_poly_quotient_ctx_set_is_pretend_field},
    {GR_METHOD_CTX_RECOVER_ZERO_DIVISOR, (gr_funcptr) _gr_poly_quotient_ctx_recover_zero_divisor},
    {GR_METHOD_CTX_IS_RATIONAL_VECTOR_SPACE, (gr_funcptr) _gr_poly_quotient_ctx_is_rational_vector_space},
    {GR_METHOD_CTX_IS_REAL_VECTOR_SPACE, (gr_funcptr) _gr_poly_quotient_ctx_is_real_vector_space},
    {GR_METHOD_CTX_IS_COMPLEX_VECTOR_SPACE, (gr_funcptr) _gr_poly_quotient_ctx_is_complex_vector_space},
    {GR_METHOD_CTX_IS_THREADSAFE,       (gr_funcptr) _gr_poly_quotient_ctx_is_threadsafe},
    {GR_METHOD_CTX_IS_FINITE,           (gr_funcptr) _gr_poly_quotient_ctx_is_finite},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC, (gr_funcptr) _gr_poly_quotient_ctx_is_finite_characteristic},
    {GR_METHOD_CTX_IS_EXACT,            (gr_funcptr) _gr_poly_quotient_ctx_is_exact},
    {GR_METHOD_CTX_IS_CANONICAL,        (gr_funcptr) _gr_poly_quotient_ctx_is_canonical},
    {GR_METHOD_CTX_SET_GEN_NAME,        (gr_funcptr) _gr_poly_quotient_ctx_set_gen_name},
    {GR_METHOD_CTX_SET_GEN_NAMES,       (gr_funcptr) _gr_poly_quotient_ctx_set_gen_names},
    {GR_METHOD_CTX_NGENS,               (gr_funcptr) gr_generic_ctx_ngens_1},
    {GR_METHOD_CTX_GEN_NAME,            (gr_funcptr) _gr_poly_quotient_ctx_gen_name},
    {GR_METHOD_CTX_BASE,                (gr_funcptr) _gr_poly_quotient_ctx_base},

    {GR_METHOD_INIT,            (gr_funcptr) _gr_poly_quotient_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_poly_quotient_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_poly_quotient_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _gr_poly_quotient_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_poly_quotient_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_poly_quotient_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_poly_quotient_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_poly_quotient_one},
    {GR_METHOD_GEN,             (gr_funcptr) _gr_poly_quotient_gen},
    {GR_METHOD_GENS,            (gr_funcptr) gr_generic_gens_single},
    {GR_METHOD_GENS_RECURSIVE,  (gr_funcptr) _gr_poly_quotient_gens_recursive},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_poly_quotient_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_poly_quotient_is_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_poly_quotient_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_poly_quotient_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_poly_quotient_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_poly_quotient_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_poly_quotient_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) _gr_poly_quotient_set_fmpq},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_poly_quotient_set_other},
    {GR_METHOD_SET_STR,         (gr_funcptr) gr_generic_set_str_balance_additions},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_poly_quotient_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_poly_quotient_add},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_poly_quotient_sub},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_poly_quotient_mul},
    {GR_METHOD_MUL_OTHER,       (gr_funcptr) _gr_poly_quotient_mul_other},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_poly_quotient_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_poly_quotient_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_poly_quotient_mul_fmpz},
    {GR_METHOD_MUL_FMPQ,        (gr_funcptr) _gr_poly_quotient_mul_fmpq},
    {GR_METHOD_DIV_UI,          (gr_funcptr) _gr_poly_quotient_div_ui},
    {GR_METHOD_DIV_SI,          (gr_funcptr) _gr_poly_quotient_div_si},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) _gr_poly_quotient_div_fmpz},
    {GR_METHOD_DIV_FMPQ,        (gr_funcptr) _gr_poly_quotient_div_fmpq},
    {GR_METHOD_SQR,             (gr_funcptr) _gr_poly_quotient_sqr},
    {GR_METHOD_POW_UI,          (gr_funcptr) _gr_poly_quotient_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) _gr_poly_quotient_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) _gr_poly_quotient_pow_fmpz},
    {GR_METHOD_INV,             (gr_funcptr) _gr_poly_quotient_inv},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_poly_quotient_is_invertible},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_poly_quotient_div},
    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_gr_poly_quotient(gr_ctx_t ctx, gr_ctx_t base, const gr_poly_t modulus)
{
    gr_poly_quotient_ctx_struct * q;
    slong sz = base->sizeof_elem;

    if (modulus->length < 1)
        flint_throw(FLINT_ERROR, "(%s): the modulus must be nonzero\n", __func__);

    if (gr_is_one(GR_ENTRY(modulus->coeffs, modulus->length - 1, sz), base) != T_TRUE)
        flint_throw(FLINT_ERROR, "(%s): the modulus must be monic\n", __func__);

    q = flint_calloc(1, sizeof(gr_poly_quotient_ctx_struct));
    q->base = base;
    q->mod = _mod_new(modulus, base);
    q->is_field = T_UNKNOWN;
    q->is_pretend_field = T_FALSE;
    q->var = (char *) default_var;
    q->version = 0;

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&q->mutex, NULL);
#endif

    ctx->which_ring = GR_CTX_GR_POLY_QUOTIENT;
    ctx->sizeof_elem = sizeof(gr_poly_struct);
    ctx->size_limit = WORD_MAX;
    GR_CTX_DATA_AS_PTR(ctx) = q;

    ctx->methods = _gr_poly_quotient_methods;

    if (!_gr_poly_quotient_methods_initialized)
    {
        gr_method_tab_init(_gr_poly_quotient_methods, _gr_poly_quotient_methods_input);
        _gr_poly_quotient_methods_initialized = 1;
    }
}
