/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_factor.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"
#include "fmpz_mod_poly.h"
#include "fmpz_mod_poly_factor.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"

typedef struct
{
    fmpz_mod_ctx_struct * ctx;
    truth_t is_prime;
}
fmpz_mod_ctx_extended_struct;

#define FMPZ_MOD_CTX(ring_ctx) ((((fmpz_mod_ctx_extended_struct *)(ring_ctx))->ctx))
#define FMPZ_MOD_IS_PRIME(ring_ctx) (((fmpz_mod_ctx_extended_struct *)(ring_ctx))->is_prime)

int
_gr_fmpz_mod_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Integers mod ");
    gr_stream_write_fmpz(out, FMPZ_MOD_CTX(ctx)->n);
    gr_stream_write(out, " (fmpz)");
    return GR_SUCCESS;
}

void
_gr_fmpz_mod_ctx_clear(gr_ctx_t ctx)
{
    fmpz_mod_ctx_clear(FMPZ_MOD_CTX(ctx));
    flint_free(FMPZ_MOD_CTX(ctx));
}

truth_t
_gr_fmpz_mod_ctx_is_field(gr_ctx_t ctx)
{
/*
    if (FMPZ_MOD_IS_PRIME(ctx) == T_UNKNOWN)
        FMPZ_MOD_IS_PRIME(ctx) = fmpz_is_prime(FMPZ_MOD_CTX(ctx)->n) ? T_TRUE : T_FALSE;
*/

    return FMPZ_MOD_IS_PRIME(ctx);
}

void
_gr_fmpz_mod_init(fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_init(x);
}

void
_gr_fmpz_mod_clear(fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_clear(x);
}

void
_gr_fmpz_mod_swap(fmpz_t x, fmpz_t y, const gr_ctx_t ctx)
{
    fmpz_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

void
_gr_fmpz_mod_set_shallow(fmpz_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    *res = *x;
}

int
_gr_fmpz_mod_randtest(fmpz_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    fmpz_mod_rand(res, state, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_write(gr_stream_t out, const fmpz_t x, const gr_ctx_t ctx)
{
    gr_stream_write_fmpz(out, x);
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_zero(fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_zero(x);
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_one(fmpz_t x, const gr_ctx_t ctx)
{
    if (fmpz_is_one(FMPZ_MOD_CTX(ctx)->n))
        fmpz_zero(x);
    else
        fmpz_one(x);
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_set_si(fmpz_t res, slong v, const gr_ctx_t ctx)
{
    fmpz_mod_set_si(res, v, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_set_ui(fmpz_t res, ulong v, const gr_ctx_t ctx)
{
    fmpz_mod_set_ui(res, v, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

/* todo: public interface */
#define NMOD_CTX_REF(ring_ctx) (((nmod_t *)((ring_ctx))))
#define NMOD_CTX(ring_ctx) (*NMOD_CTX_REF(ring_ctx))

int
_gr_fmpz_mod_set_other(fmpz_t res, gr_ptr v, gr_ctx_t v_ctx, const gr_ctx_t ctx)
{
    if (v_ctx->which_ring == GR_CTX_FMPZ_MOD)
    {
        if (!fmpz_equal(FMPZ_MOD_CTX(ctx)->n, FMPZ_MOD_CTX(v_ctx)->n))
            return GR_DOMAIN;

        fmpz_set(res, v);
        return GR_SUCCESS;
    }

    if (v_ctx->which_ring == GR_CTX_NMOD)
    {
        if (!fmpz_equal_ui(FMPZ_MOD_CTX(ctx)->n, NMOD_CTX(v_ctx).n))
            return GR_DOMAIN;

        fmpz_set_ui(res, * (ulong *) v);
        return GR_SUCCESS;
    }

    return GR_UNABLE;
}

int
_gr_fmpz_mod_set_fmpz(fmpz_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    fmpz_mod_set_fmpz(res, v, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

truth_t
_gr_fmpz_mod_is_zero(const fmpz_t x, const gr_ctx_t ctx)
{
    return fmpz_is_zero(x) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fmpz_mod_is_one(const fmpz_t x, const gr_ctx_t ctx)
{
    return fmpz_mod_is_one(x, FMPZ_MOD_CTX(ctx)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fmpz_mod_is_neg_one(const fmpz_t x, const gr_ctx_t ctx)
{
    truth_t res;
    fmpz_t t;
    fmpz_init(t);
    fmpz_mod_set_si(t, -1, FMPZ_MOD_CTX(ctx));
    res = fmpz_equal(t, x) ? T_TRUE : T_FALSE;
    fmpz_clear(t);
    return res;
}

truth_t
_gr_fmpz_mod_equal(const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    return fmpz_equal(x, y) ? T_TRUE : T_FALSE;
}

int
_gr_fmpz_mod_set(fmpz_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_set(res, x);
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_neg(fmpz_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_mod_neg(res, x, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_add(fmpz_t res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fmpz_mod_add(res, x, y, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_add_si(fmpz_t res, const fmpz_t x, slong y, const gr_ctx_t ctx)
{
    fmpz_mod_add_si(res, x, y, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_add_ui(fmpz_t res, const fmpz_t x, ulong y, const gr_ctx_t ctx)
{
    fmpz_mod_add_ui(res, x, y, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_sub(fmpz_t res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fmpz_mod_sub(res, x, y, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_mul(fmpz_t res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
#if 1
    fmpz_mod_mul(res, x, y, FMPZ_MOD_CTX(ctx));
#else
    fmpz_mul(res, x, y);
    fmpz_mod(res, res, FMPZ_MOD_CTX(ctx)->n);
#endif
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_mul_si(fmpz_t res, const fmpz_t x, slong y, const gr_ctx_t ctx)
{
    fmpz_mod_mul_si(res, x, y, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_addmul(fmpz_t res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_mul(t, x, y);
    fmpz_add(t, t, res);
    fmpz_mod_set_fmpz(res, t, FMPZ_MOD_CTX(ctx));
    fmpz_clear(t);
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_submul(fmpz_t res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_mul(t, x, y);
    fmpz_sub(t, res, t);
    fmpz_mod_set_fmpz(res, t, FMPZ_MOD_CTX(ctx));
    fmpz_clear(t);
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_mul_two(fmpz_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_mod_add(res, x, x, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_sqr(fmpz_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    fmpz_mod_mul(res, x, x, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_inv(fmpz_t res, const fmpz_t x, const gr_ctx_t ctx)
{
    /* todo: also check for -1 when fast? */
    if (fmpz_is_one(x))
    {
        fmpz_one(res);
        return GR_SUCCESS;
    }
    else if (fmpz_is_zero(x))
    {
        fmpz_zero(res);
        return fmpz_is_one(FMPZ_MOD_CTX(ctx)->n) ? GR_SUCCESS : GR_DOMAIN;
    }
    else
    {
        int status;

        fmpz_t d;
        fmpz_init(d);
        fmpz_gcdinv(d, res, x, FMPZ_MOD_CTX(ctx)->n);

        if (fmpz_is_one(d))
            status = GR_SUCCESS;
        else
            status = GR_DOMAIN;

        fmpz_clear(d);
        return status;
    }
}

int
_gr_fmpz_mod_div(fmpz_t res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    int status;
    fmpz_t t;

    fmpz_init(t);
    status = _gr_fmpz_mod_inv(t, y, ctx);
    if (status == GR_SUCCESS)
        fmpz_mod_mul(res, x, t, FMPZ_MOD_CTX(ctx));
    else
        fmpz_zero(res);
    fmpz_clear(t);

    return status;
}

int
_gr_fmpz_mod_div_nonunique(fmpz_t res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    int status;

#if 1
    status = fmpz_mod_divides(res, x, y, FMPZ_MOD_CTX(ctx)) ? GR_SUCCESS : GR_DOMAIN;
#else
    if (FMPZ_MOD_IS_PRIME(ctx) != T_TRUE)
    {
        status = fmpz_mod_divides(res, x, y, FMPZ_MOD_CTX(ctx)) ? GR_SUCCESS : GR_DOMAIN;
    }
    else
    {
        fmpz_t t;
        fmpz_init(t);
        status = _gr_fmpz_mod_inv(t, y, ctx);
        if (status == GR_SUCCESS)
            fmpz_mod_mul(res, x, t, FMPZ_MOD_CTX(ctx));
        else
            fmpz_zero(res);
        fmpz_clear(t);
    }
#endif

    return status;
}

truth_t
_gr_fmpz_mod_divides(const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    truth_t res;
    fmpz_t t;
    fmpz_init(t);
    res = fmpz_mod_divides(t, y, x, FMPZ_MOD_CTX(ctx)) ? T_TRUE : T_FALSE;
    fmpz_clear(t);
    return res;
}

truth_t
_gr_fmpz_mod_is_invertible(const fmpz_t x, const gr_ctx_t ctx)
{
    return fmpz_mod_is_invertible(x, FMPZ_MOD_CTX(ctx)) ? T_TRUE : T_FALSE;
}

int
_gr_fmpz_mod_pow_ui(fmpz_t res, const fmpz_t x, ulong exp, const gr_ctx_t ctx)
{
    fmpz_mod_pow_ui(res, x, exp, FMPZ_MOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mod_pow_fmpz(fmpz_t res, const fmpz_t x, const fmpz_t exp, const gr_ctx_t ctx)
{
    return fmpz_mod_pow_fmpz(res, x, exp, FMPZ_MOD_CTX(ctx)) ? GR_SUCCESS : GR_DOMAIN;
}

/* todo: len 1 */
int
_gr_fmpz_mod_vec_dot(fmpz_t res, const fmpz_t initial, int subtract, const fmpz * vec1, const fmpz * vec2, slong len, gr_ctx_t ctx)
{
    slong i;

    if (len <= 0)
    {
        if (initial == NULL)
            fmpz_zero(res);
        else
            fmpz_set(res, initial);
        return GR_SUCCESS;
    }

    if (initial == NULL)
    {
        fmpz_mul(res, vec1, vec2);
    }
    else
    {
        if (subtract)
            fmpz_neg(res, initial);
        else
            fmpz_set(res, initial);

        fmpz_addmul(res, vec1, vec2);
    }

    for (i = 1; i < len; i++)
        fmpz_addmul(res, vec1 + i, vec2 + i);

    if (subtract)
        fmpz_neg(res, res);

    fmpz_mod_set_fmpz(res, res, FMPZ_MOD_CTX(ctx));

    return GR_SUCCESS;
}

/* todo: len 1 */
int
_gr_fmpz_mod_vec_dot_rev(fmpz_t res, const fmpz_t initial, int subtract, const fmpz * vec1, const fmpz * vec2, slong len, gr_ctx_t ctx)
{
    slong i;

    if (len <= 0)
    {
        if (initial == NULL)
            fmpz_zero(res);
        else
            fmpz_set(res, initial);
        return GR_SUCCESS;
    }

    if (initial == NULL)
    {
        fmpz_mul(res, vec1, vec2 + len - 1);
    }
    else
    {
        if (subtract)
            fmpz_neg(res, initial);
        else
            fmpz_set(res, initial);

        fmpz_addmul(res, vec1, vec2 + len - 1);
    }

    for (i = 1; i < len; i++)
        fmpz_addmul(res, vec1 + i, vec2 + len - 1 - i);

    if (subtract)
        fmpz_neg(res, res);

    fmpz_mod_set_fmpz(res, res, FMPZ_MOD_CTX(ctx));

    return GR_SUCCESS;
}

int
_gr_fmpz_mod_poly_mullow(fmpz * res,
    const fmpz * poly1, slong len1,
    const fmpz * poly2, slong len2, slong n, gr_ctx_t ctx)
{
    if (len1 >= len2)
        _fmpz_mod_poly_mullow(res, poly1, len1, poly2, len2, n, FMPZ_MOD_CTX(ctx));
    else
        _fmpz_mod_poly_mullow(res, poly2, len2, poly1, len1, n, FMPZ_MOD_CTX(ctx));

    return GR_SUCCESS;
}

/* fixme: duplicate _fmpz_mod_poly methods for error handling */
int
_gr_fmpz_mod_poly_divrem(fmpz * Q, fmpz * R, const fmpz * A, slong lenA,
                                  const fmpz * B, slong lenB, gr_ctx_t ctx)
{
    if (lenB <= 30 || lenA - lenB <= 5)
    {
        fmpz_t invB;
        int status;

        fmpz_init(invB);
        status = _gr_fmpz_mod_inv(invB, B + lenB - 1, ctx);
        if (status == GR_SUCCESS)
            _fmpz_mod_poly_divrem_basecase(Q, R, A, lenA, B, lenB, invB, FMPZ_MOD_CTX(ctx));

        fmpz_clear(invB);
        return status;
    }
    else
    {
        return _gr_poly_divrem_newton(Q, R, A, lenA, B, lenB, ctx);
    }
}

#define TUNE_TAB_SIZE 23

static const int tuning_bit_steps[TUNE_TAB_SIZE] = { 32, 45, 64, 91, 128, 181, 256, 362, 512, 724, 1024, 1448, 2048, 2896, 4096, 5793, 8192, 11585, 16384, 23170, 32768, 46341, 65536};
static const short inv_series_cutoff_tab[TUNE_TAB_SIZE] = {21, 14, 40, 39, 48, 60, 89, 72, 72, 54, 48, 39, 32, 24, 24, 20, 17, 18, 16, 15, 13, 12, 14, };
static const short div_series_cutoff_tab[TUNE_TAB_SIZE] = {23, 21, 52, 50, 66, 101, 106, 97, 106, 72, 60, 50, 44, 35, 38, 30, 26, 22, 20, 18, 16, 14, 22, };

static const slong find_cutoff(const short * tab, slong b)
{
    slong i;

    i = 0;
    while (i + 1 < TUNE_TAB_SIZE && tuning_bit_steps[i + 1] <= b)
        i++;

    return tab[i];
}

int
_gr_fmpz_mod_poly_inv_series(fmpz * Q, const fmpz * B, slong lenB, slong len, gr_ctx_t ctx)
{
    slong cutoff, bits;
    lenB = FLINT_MIN(len, lenB);

    if (lenB <= 20)
        return _gr_poly_inv_series_basecase(Q, B, lenB, len, ctx);

    bits = fmpz_bits(fmpz_mod_ctx_modulus(FMPZ_MOD_CTX(ctx)));
    cutoff  = find_cutoff(inv_series_cutoff_tab, bits);

    if (lenB <= cutoff)
        return _gr_poly_inv_series_basecase(Q, B, lenB, len, ctx);
    else
        return _gr_poly_inv_series_newton(Q, B, lenB, len, cutoff, ctx);
}

/* todo: the fmpz_mod_poly module has better basecase code */
int
_gr_fmpz_mod_poly_div_series(fmpz * Q, const fmpz * A, slong lenA, const fmpz * B, slong lenB, slong len, gr_ctx_t ctx)
{
    slong cutoff, bits;

    lenA = FLINT_MIN(len, lenA);
    lenB = FLINT_MIN(len, lenB);

    if (lenB <= 20)
        return _gr_poly_div_series_basecase(Q, A, lenA, B, lenB, len, ctx);

    bits = fmpz_bits(fmpz_mod_ctx_modulus(FMPZ_MOD_CTX(ctx)));
    cutoff  = find_cutoff(div_series_cutoff_tab, bits);

    if (lenB <= cutoff)
        return _gr_poly_div_series_basecase(Q, A, lenA, B, lenB, len, ctx);
    else
        return _gr_poly_div_series_newton(Q, A, lenA, B, lenB, len, cutoff, ctx);
}


/* todo: also need the _other version ... ? */
/* todo: implement generically */

int
_gr_fmpz_mod_roots_gr_poly(gr_vec_t roots, gr_vec_t mult, const fmpz_mod_poly_t poly, int flags, gr_ctx_t ctx)
{
    if (poly->length == 0)
        return GR_DOMAIN;

    {
        gr_ctx_t ZZ;
        fmpz_mod_poly_factor_t fac;
        slong i, num;
        int status = GR_SUCCESS;

        gr_ctx_init_fmpz(ZZ);

        fmpz_mod_poly_factor_init(fac, FMPZ_MOD_CTX(ctx));

        if (gr_ctx_is_field(ctx) == T_TRUE)
        {
            fmpz_mod_poly_roots(fac, poly, 1, FMPZ_MOD_CTX(ctx));
        }
        else
        {
            fmpz_factor_t nfac;
            fmpz_factor_init(nfac);
            fmpz_factor(nfac, FMPZ_MOD_CTX(ctx)->n);

            num = 0;
            for (i = 0; i < nfac->num; i++)
                num += nfac->exp[i];

            if (num > 20)
            {
                status = GR_UNABLE;
            }
            else
            {
                if (!fmpz_mod_poly_roots_factored_with_length_limit(fac, poly, 1, 1000000, nfac, FMPZ_MOD_CTX(ctx)))
                    status = GR_UNABLE;
            }

            fmpz_factor_clear(nfac);
        }

        if (status == GR_SUCCESS)
        {
            num = fac->num;

            gr_vec_set_length(roots, num, ctx);
            gr_vec_set_length(mult, num, ZZ);

            for (i = 0; i < num; i++)
            {
                fmpz_mod_neg(gr_vec_entry_ptr(roots, i, ctx), fac->poly[i].coeffs, FMPZ_MOD_CTX(ctx));

                /* work around flint bug: factors can be non-monic */
                if (!fmpz_mod_is_one(fac->poly[i].coeffs + 1, FMPZ_MOD_CTX(ctx)))
                    status |= _gr_fmpz_mod_div(gr_vec_entry_ptr(roots, i, ctx), gr_vec_entry_ptr(roots, i, ctx), fac->poly[i].coeffs + 1, ctx);

                fmpz_set_ui(((fmpz *) mult->entries) + i, fac->exp[i]);
            }
        }

        fmpz_mod_poly_factor_clear(fac, FMPZ_MOD_CTX(ctx));
        gr_ctx_clear(ZZ);

        return status;
    }
}

int
_gr_fmpz_mod_mat_mul(fmpz_mat_t res, const fmpz_mat_t x, const fmpz_mat_t y, gr_ctx_t ctx)
{
    fmpz_mod_mat_t MM;

    fmpz_mat_mul(res, x, y);

    *MM->mat = *res;
    *MM->mod = *(FMPZ_MOD_CTX(ctx)->n);

    _fmpz_mod_mat_reduce(MM);

    return GR_SUCCESS;
}

int _fmpz_mod_methods_initialized = 0;

gr_static_method_table _fmpz_mod_methods;

gr_method_tab_input _fmpz_mod_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_fmpz_mod_ctx_write},
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_fmpz_mod_ctx_clear},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) _gr_fmpz_mod_ctx_is_field},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) _gr_fmpz_mod_ctx_is_field},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,
                                (gr_funcptr) _gr_fmpz_mod_ctx_is_field},
    {GR_METHOD_CTX_IS_FINITE,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_INIT,            (gr_funcptr) _gr_fmpz_mod_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_fmpz_mod_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_fmpz_mod_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _gr_fmpz_mod_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_fmpz_mod_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_fmpz_mod_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_fmpz_mod_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_fmpz_mod_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_fmpz_mod_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_fmpz_mod_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_fmpz_mod_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_fmpz_mod_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_fmpz_mod_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_fmpz_mod_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_fmpz_mod_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_fmpz_mod_set_fmpz},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_fmpz_mod_set_other},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_fmpz_mod_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_fmpz_mod_add},
    {GR_METHOD_ADD_UI,          (gr_funcptr) _gr_fmpz_mod_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) _gr_fmpz_mod_add_si},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_fmpz_mod_sub},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_fmpz_mod_mul},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_fmpz_mod_mul_si},
    {GR_METHOD_ADDMUL,          (gr_funcptr) _gr_fmpz_mod_addmul},
    {GR_METHOD_SUBMUL,          (gr_funcptr) _gr_fmpz_mod_submul},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) _gr_fmpz_mod_mul_two},
    {GR_METHOD_SQR,             (gr_funcptr) _gr_fmpz_mod_sqr},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_fmpz_mod_div},
    {GR_METHOD_DIV_NONUNIQUE,   (gr_funcptr) _gr_fmpz_mod_div_nonunique},
    {GR_METHOD_DIVIDES,         (gr_funcptr) _gr_fmpz_mod_divides},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_fmpz_mod_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) _gr_fmpz_mod_inv},
    {GR_METHOD_POW_UI,          (gr_funcptr) _gr_fmpz_mod_pow_ui},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) _gr_fmpz_mod_pow_fmpz},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) _gr_fmpz_mod_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) _gr_fmpz_mod_vec_dot_rev},
    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _gr_fmpz_mod_poly_mullow},
    {GR_METHOD_POLY_INV_SERIES, (gr_funcptr) _gr_fmpz_mod_poly_inv_series},
    {GR_METHOD_POLY_DIV_SERIES, (gr_funcptr) _gr_fmpz_mod_poly_div_series},
    {GR_METHOD_POLY_DIVREM,     (gr_funcptr) _gr_fmpz_mod_poly_divrem},
    {GR_METHOD_POLY_ROOTS,      (gr_funcptr) _gr_fmpz_mod_roots_gr_poly},
    {GR_METHOD_MAT_MUL,         (gr_funcptr) _gr_fmpz_mod_mat_mul},
    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_fmpz_mod(gr_ctx_t ctx, const fmpz_t n)
{
    ctx->which_ring = GR_CTX_FMPZ_MOD;
    ctx->sizeof_elem = sizeof(fmpz);

    FMPZ_MOD_CTX(ctx) = flint_malloc(sizeof(fmpz_mod_ctx_struct));
    fmpz_mod_ctx_init(FMPZ_MOD_CTX(ctx), n);
    FMPZ_MOD_IS_PRIME(ctx) = T_UNKNOWN;

    ctx->size_limit = WORD_MAX;

    ctx->methods = _fmpz_mod_methods;

    if (!_fmpz_mod_methods_initialized)
    {
        gr_method_tab_init(_fmpz_mod_methods, _fmpz_mod_methods_input);
        _fmpz_mod_methods_initialized = 1;
    }
}

void
_gr_ctx_init_fmpz_mod_from_ref(gr_ctx_t ctx, const void * fctx)
{
    ctx->which_ring = GR_CTX_FMPZ_MOD;
    ctx->sizeof_elem = sizeof(fmpz);

    FMPZ_MOD_CTX(ctx) = (fmpz_mod_ctx_struct *) fctx;
    FMPZ_MOD_IS_PRIME(ctx) = T_UNKNOWN;

    ctx->size_limit = WORD_MAX;

    ctx->methods = _fmpz_mod_methods;

    if (!_fmpz_mod_methods_initialized)
    {
        gr_method_tab_init(_fmpz_mod_methods, _fmpz_mod_methods_input);
        _fmpz_mod_methods_initialized = 1;
    }
}

void
gr_ctx_fmpz_mod_set_primality(gr_ctx_t ctx, truth_t is_prime)
{
    FMPZ_MOD_IS_PRIME(ctx) = is_prime;
}
