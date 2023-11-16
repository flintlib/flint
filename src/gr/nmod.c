/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpq.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_mat.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_mat.h"
#include "ulong_extras.h"

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

void
_gr_nmod_set_shallow(ulong * res, const ulong * x, const gr_ctx_t ctx)
{
    *res = *x;
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
    gr_stream_write_ui(out, x[0]);
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
    res[0] = nmod_set_si(v, NMOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nmod_set_ui(ulong * res, ulong v, const gr_ctx_t ctx)
{
    res[0] = nmod_set_ui(v, NMOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nmod_set_fmpz(ulong * res, const fmpz_t v, const gr_ctx_t ctx)
{
    nmod_t mod = NMOD_CTX(ctx);
    res[0] = fmpz_get_nmod(v, mod);
    return GR_SUCCESS;
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

#include "fmpz_mod.h"

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

    if (v_ctx->which_ring == GR_CTX_FMPZ)
    {
        res[0] = fmpz_get_nmod(v, NMOD_CTX(ctx));
        return GR_SUCCESS;
    }

    if (v_ctx->which_ring == GR_CTX_FMPQ)
    {
        ulong a, b;
        int status;
        const fmpq * q = v;

        if (fmpz_is_one(fmpq_denref(q)))
        {
            res[0] = fmpz_get_nmod(fmpq_numref(q), NMOD_CTX(ctx));
        }
        else
        {
            b = fmpz_get_nmod(fmpq_denref(q), NMOD_CTX(ctx));
            status = _gr_nmod_inv(&b, &b, ctx);
            if (status != GR_SUCCESS)
                return status;

            a = fmpz_get_nmod(fmpq_numref(q), NMOD_CTX(ctx));
            res[0] = nmod_mul(a, b, NMOD_CTX(ctx));
        }

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
_gr_nmod_add_ui(ulong * res, const ulong * x, ulong y, const gr_ctx_t ctx)
{
    res[0] = nmod_add(x[0], nmod_set_ui(y, NMOD_CTX(ctx)), NMOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nmod_add_si(ulong * res, const ulong * x, slong y, const gr_ctx_t ctx)
{
    res[0] = nmod_add(x[0], nmod_set_si(y, NMOD_CTX(ctx)), NMOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nmod_sub(ulong * res, const ulong * x, const ulong * y, const gr_ctx_t ctx)
{
    res[0] = nmod_sub(x[0], y[0], NMOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nmod_sub_ui(ulong * res, const ulong * x, ulong y, const gr_ctx_t ctx)
{
    res[0] = nmod_sub(x[0], nmod_set_ui(y, NMOD_CTX(ctx)), NMOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nmod_sub_si(ulong * res, const ulong * x, slong y, const gr_ctx_t ctx)
{
    res[0] = nmod_sub(x[0], nmod_set_si(y, NMOD_CTX(ctx)), NMOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nmod_mul(ulong * res, const ulong * x, const ulong * y, const gr_ctx_t ctx)
{
    res[0] = nmod_mul(x[0], y[0], NMOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nmod_mul_ui(ulong * res, const ulong * x, ulong y, const gr_ctx_t ctx)
{
/*
    res[0] = nmod_mul(x[0], nmod_set_ui(y, NMOD_CTX(ctx)), NMOD_CTX(ctx));
*/
    res[0] = n_mulmod2_preinv(x[0], y, NMOD_CTX(ctx).n, NMOD_CTX(ctx).ninv);
    return GR_SUCCESS;
}

int
_gr_nmod_mul_si(ulong * res, const ulong * x, slong y, const gr_ctx_t ctx)
{
/*
    res[0] = nmod_mul(x[0], nmod_set_si(y, NMOD_CTX(ctx)), NMOD_CTX(ctx));
*/
    if (y >= 0)
        res[0] = n_mulmod2_preinv(x[0], y, NMOD_CTX(ctx).n, NMOD_CTX(ctx).ninv);
    else
        res[0] = nmod_neg(n_mulmod2_preinv(x[0], -y, NMOD_CTX(ctx).n, NMOD_CTX(ctx).ninv), NMOD_CTX(ctx));

    return GR_SUCCESS;
}

int
_gr_nmod_mul_fmpz(ulong * res, const ulong * x, const fmpz_t y, const gr_ctx_t ctx)
{
    if (!COEFF_IS_MPZ(*y))
    {
        if (*y >= 0)
            res[0] = n_mulmod2_preinv(x[0], *y, NMOD_CTX(ctx).n, NMOD_CTX(ctx).ninv);
        else
            res[0] = nmod_neg(n_mulmod2_preinv(x[0], -*y, NMOD_CTX(ctx).n, NMOD_CTX(ctx).ninv), NMOD_CTX(ctx));
    }
    else
    {
        res[0] = nmod_mul(x[0], fmpz_get_nmod(y, NMOD_CTX(ctx)), NMOD_CTX(ctx));
    }

    return GR_SUCCESS;
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
_gr_nmod_div(ulong * res, const ulong * x, const ulong * y, const gr_ctx_t ctx)
{
    ulong t;
    int status;

    status = _gr_nmod_inv(&t, y, ctx);

    if (status == GR_SUCCESS)
        _gr_nmod_mul(res, x, &t, ctx);

    return status;
}

int
_gr_nmod_div_si(ulong * res, const ulong * x, slong y, const gr_ctx_t ctx)
{
    ulong c = nmod_set_si(y, NMOD_CTX(ctx));

    return _gr_nmod_div(res, x, &c, ctx);
}

int
_gr_nmod_div_ui(ulong * res, const ulong * x, ulong y, const gr_ctx_t ctx)
{
    ulong c = nmod_set_ui(y, NMOD_CTX(ctx));

    return _gr_nmod_div(res, x, &c, ctx);
}

int
_gr_nmod_div_fmpz(ulong * res, const ulong * x, const fmpz_t y, const gr_ctx_t ctx)
{
    ulong c = fmpz_get_nmod(y, NMOD_CTX(ctx));

    return _gr_nmod_div(res, x, &c, ctx);
}

truth_t
_gr_nmod_is_invertible(const ulong * x, const gr_ctx_t ctx)
{
    ulong r, g;
    g = n_gcdinv(&r, x[0], NMOD_CTX(ctx).n);
    return (g == 1) ? T_TRUE : T_FALSE;
}

truth_t
_gr_nmod_divides(const ulong * x, const ulong * y, const gr_ctx_t ctx)
{
    ulong t;
    return nmod_divides(&t, y[0], x[0], NMOD_CTX(ctx)) ? T_TRUE : T_FALSE;
}

int
_gr_nmod_div_nonunique(ulong * res, const ulong * x, const ulong * y, const gr_ctx_t ctx)
{
    ulong t;
    int status;

    status = _gr_nmod_inv(&t, y, ctx);

    if (status == GR_SUCCESS)
    {
        _gr_nmod_mul(res, x, &t, ctx);
    }
    else
    {
        status = nmod_divides(res, *x, *y, NMOD_CTX(ctx)) ? GR_SUCCESS : GR_DOMAIN;
    }

    return status;
}

int
_gr_nmod_mul_2exp_si(ulong * res, ulong * x, slong y, const gr_ctx_t ctx)
{
    ulong c, m = NMOD_CTX(ctx).n;

    if (y >= 0)
    {
        if (y < FLINT_BITS)
        {
            c = UWORD(1) << y;
            if (c >= m)
                NMOD_RED(c, c, NMOD_CTX(ctx));
        }
        else
        {
            /* accidentally also works when mod <= 2 */
            c = nmod_pow_ui(2, y, NMOD_CTX(ctx));
        }
    }
    else
    {
        if (m % 2 == 0)
        {
            if (m == 1)
            {
                res[0] = 0;
                return GR_SUCCESS;
            }

            return GR_DOMAIN;
        }

        /* quickly construct 1/2 */
        c = (m - 1) / 2 + 1;

        c = nmod_pow_ui(c, -y, NMOD_CTX(ctx));
    }

    res[0] = nmod_mul(x[0], c, NMOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nmod_sqrt(ulong * res, const ulong * x, gr_ctx_t ctx)
{
    if (x[0] <= 1)
    {
        res[0] = x[0];
        return GR_SUCCESS;
    }

    /* todo: caching prime status */
    /* todo: handle non-primes */
    if (!n_is_prime(NMOD_CTX(ctx).n))
        return GR_UNABLE;

    res[0] = n_sqrtmod(x[0], NMOD_CTX(ctx).n);

    if (res[0] == 0)
        return GR_DOMAIN;
    else
        return GR_SUCCESS;
}

/* todo: pow_ui, ... */

void
_gr_nmod_vec_init(ulong * res, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        res[i] = 0;
}

void
_gr_nmod_vec_clear(ulong * res, slong len, gr_ctx_t ctx)
{
}

int
_gr_nmod_vec_set(ulong * res, const ulong * vec, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        res[i] = vec[i];

    return GR_SUCCESS;
}

int
_gr_nmod_vec_normalise(slong * res, const ulong * vec, slong len, gr_ctx_t ctx)
{
    while (len > 0 && vec[len - 1] == 0)
        len--;

    res[0] = len;
    return GR_SUCCESS;
}

slong
_gr_nmod_vec_normalise_weak(const ulong * vec, slong len, gr_ctx_t ctx)
{
    while (len > 0 && vec[len - 1] == 0)
        len--;

    return len;
}

int
_gr_nmod_vec_neg(ulong * res, const ulong * vec, slong len, gr_ctx_t ctx)
{
    slong i;
    nmod_t mod = NMOD_CTX(ctx);

    for (i = 0; i < len; i++)
        res[i] = nmod_neg(vec[i], mod);

    return GR_SUCCESS;
}

int
_gr_nmod_vec_add(ulong * res, const ulong * vec1, const ulong * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    nmod_t mod = NMOD_CTX(ctx);

    for (i = 0; i < len; i++)
        res[i] = nmod_add(vec1[i], vec2[i], mod);

    return GR_SUCCESS;
}

int
_gr_nmod_vec_sub(ulong * res, const ulong * vec1, const ulong * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    nmod_t mod = NMOD_CTX(ctx);

    for (i = 0; i < len; i++)
        res[i] = nmod_sub(vec1[i], vec2[i], mod);

    return GR_SUCCESS;
}


static inline void _nmod_vec_scalar_mul_nmod_fullword_inline(mp_ptr res, mp_srcptr vec,
                               slong len, mp_limb_t c, nmod_t mod)
{
    slong i;

    for (i = 0; i < len; i++)
        NMOD_MUL_FULLWORD(res[i], vec[i], c, mod);
}

static inline void _nmod_vec_scalar_mul_nmod_generic_inline(mp_ptr res, mp_srcptr vec,
                               slong len, mp_limb_t c, nmod_t mod)
{
    slong i;

    for (i = 0; i < len; i++)
        NMOD_MUL_PRENORM(res[i], vec[i], c << mod.norm, mod);
}

static inline void _nmod_vec_scalar_mul_nmod_inline(mp_ptr res, mp_srcptr vec,
                               slong len, mp_limb_t c, nmod_t mod)
{
    if (NMOD_BITS(mod) == FLINT_BITS)
        _nmod_vec_scalar_mul_nmod_fullword_inline(res, vec, len, c, mod);
    else if (len > 10)
        _nmod_vec_scalar_mul_nmod_shoup(res, vec, len, c, mod);
    else
        _nmod_vec_scalar_mul_nmod_generic_inline(res, vec, len, c, mod);
}

int
_gr_nmod_vec_mul_scalar(ulong * res, const ulong * vec1, slong len, const ulong * c, gr_ctx_t ctx)
{
    nmod_t mod = NMOD_CTX(ctx);
    _nmod_vec_scalar_mul_nmod_inline(res, vec1, len, c[0], mod);
    return GR_SUCCESS;
}

int
_gr_nmod_vec_mul_scalar_si(ulong * res, const ulong * vec1, slong len, slong c, gr_ctx_t ctx)
{
    nmod_t mod = NMOD_CTX(ctx);
    _nmod_vec_scalar_mul_nmod_inline(res, vec1, len, nmod_set_si(c, mod), NMOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nmod_vec_mul_scalar_ui(ulong * res, const ulong * vec1, slong len, ulong c, gr_ctx_t ctx)
{
    nmod_t mod = NMOD_CTX(ctx);
    _nmod_vec_scalar_mul_nmod_inline(res, vec1, len, nmod_set_ui(c, mod), NMOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nmod_vec_mul_scalar_fmpz(ulong * res, const ulong * vec1, slong len, const fmpz_t c, gr_ctx_t ctx)
{
    nmod_t mod = NMOD_CTX(ctx);
    _nmod_vec_scalar_mul_nmod(res, vec1, len, fmpz_get_nmod(c, mod), NMOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nmod_vec_mul_scalar_2exp_si(ulong * res, const ulong * vec1, slong len, slong c, gr_ctx_t ctx)
{
    ulong t[1];
    int status = GR_SUCCESS;

    if (c == 1)
        return _gr_nmod_vec_add(res, vec1, vec1, len, ctx);

    status |= _gr_nmod_one(t, ctx);
    status |= _gr_nmod_mul_2exp_si(t, t, c, ctx);
    status |= _gr_nmod_vec_mul_scalar(res, vec1, len, t, ctx);
    return status;
}

int
_gr_nmod_vec_addmul_scalar(ulong * res, const ulong * vec1, slong len, const ulong * c, gr_ctx_t ctx)
{
    nmod_t mod = NMOD_CTX(ctx);
    _nmod_vec_scalar_addmul_nmod(res, vec1, len, c[0], mod);
    return GR_SUCCESS;
}

int
_gr_nmod_vec_submul_scalar(ulong * res, const ulong * vec1, slong len, const ulong * c, gr_ctx_t ctx)
{
    nmod_t mod = NMOD_CTX(ctx);
    _nmod_vec_scalar_addmul_nmod(res, vec1, len, nmod_neg(c[0], mod), mod);
    return GR_SUCCESS;
}

int
_gr_nmod_vec_addmul_scalar_si(ulong * res, const ulong * vec1, slong len, slong c, gr_ctx_t ctx)
{
    nmod_t mod = NMOD_CTX(ctx);
    _nmod_vec_scalar_addmul_nmod(res, vec1, len, nmod_set_si(c, mod), mod);
    return GR_SUCCESS;
}

int
_gr_nmod_vec_submul_scalar_si(ulong * res, const ulong * vec1, slong len, slong c, gr_ctx_t ctx)
{
    nmod_t mod = NMOD_CTX(ctx);
    _nmod_vec_scalar_addmul_nmod(res, vec1, len, nmod_neg(nmod_set_si(c, mod), mod), mod);
    return GR_SUCCESS;
}

int
_gr_nmod_vec_sum(ulong * res, const ulong * vec, slong len, gr_ctx_t ctx)
{
    ulong hi, lo;
    slong i;
    nmod_t mod = NMOD_CTX(ctx);

    if (len < 10)
    {
        lo = 0;
        for (i = 0; i < len; i++)
            lo = nmod_add(lo, vec[i], mod);
    }
    else
    {
        umul_ppmm(hi, lo, mod.n, len);

        if (hi == 0)
        {
            lo = vec[0];
            for (i = 1; i < len; i++)
                lo += vec[i];

            NMOD_RED(lo, lo, mod);
        }
        else
        {
            lo = vec[0];
            hi = 0;

            for (i = 1; i < len; i++)
                add_ssaaaa(hi, lo, hi, lo, 0, vec[i]);

            NMOD2_RED2(lo, hi, lo, mod);
        }
    }

    res[0] = lo;
    return GR_SUCCESS;
}

int
_gr_nmod_vec_product(ulong * res, const ulong * vec, slong len, gr_ctx_t ctx)
{
    nmod_t mod = NMOD_CTX(ctx);

    if (len <= 2)
    {
        if (len == 2)
            res[0] = nmod_mul(vec[0], vec[1], mod);
        else if (len == 1)
            res[0] = vec[0];
        else
            res[0] = (mod.n != 1);
    }
    else
    {
        ulong p;
        slong i;

        if (mod.norm == 0)
        {
            p = _nmod_mul_fullword(vec[0], vec[1], mod);

            for (i = 2; i < len; i++)
                p = _nmod_mul_fullword(p, vec[i], mod);
        }
        else
        {
            p = nmod_mul(vec[0], vec[1], mod);

            for (i = 2; i < len; i++)
                p = nmod_mul(p, vec[i], mod);
        }

        res[0] = p;
    }

    return GR_SUCCESS;
}


int
__gr_nmod_vec_dot(ulong * res, const ulong * initial, int subtract, const ulong * vec1, const ulong * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    ulong s;
    int nlimbs;
    nmod_t mod;

    if (len <= 1)
    {
        if (len == 2)   /* todo: fmma */
        {
            mod = NMOD_CTX(ctx);
            s = nmod_mul(vec1[0], vec2[0], mod);
            s = nmod_addmul(s, vec1[1], vec2[1], mod);
        }
        else if (len == 1)
        {
            mod = NMOD_CTX(ctx);
            s = nmod_mul(vec1[0], vec2[0], mod);
        }
        else
        {
            if (initial == NULL)
                _gr_nmod_zero(res, ctx);
            else
                _gr_nmod_set(res, initial, ctx);
            return GR_SUCCESS;
        }
    }
    else
    {
        mod = NMOD_CTX(ctx);

        if (len <= 16)
        {
            if (mod.n <= UWORD(1) << (FLINT_BITS / 2 - 2))
                nlimbs = 1;
            if (mod.n <= UWORD(1) << (FLINT_BITS - 2))
                nlimbs = 2;
            else
                nlimbs = 3;
        }
        else
        {
            nlimbs = _nmod_vec_dot_bound_limbs(len, mod);
        }

        NMOD_VEC_DOT(s, i, len, vec1[i], vec2[i], mod, nlimbs);
    }

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

    if (len <= 1)
    {
        if (len == 2)   /* todo: fmma */
        {
            mod = NMOD_CTX(ctx);
            s = nmod_mul(vec1[0], vec2[1], mod);
            s = nmod_addmul(s, vec1[1], vec2[0], mod);
        }
        else if (len == 1)
        {
            mod = NMOD_CTX(ctx);
            s = nmod_mul(vec1[0], vec2[0], mod);
        }
        else
        {
            if (initial == NULL)
                _gr_nmod_zero(res, ctx);
            else
                _gr_nmod_set(res, initial, ctx);
            return GR_SUCCESS;
        }
    }
    else
    {
        mod = NMOD_CTX(ctx);

        if (len <= 16)
        {
            if (mod.n <= UWORD(1) << (FLINT_BITS / 2 - 2))
                nlimbs = 1;
            if (mod.n <= UWORD(1) << (FLINT_BITS - 2))
                nlimbs = 2;
            else
                nlimbs = 3;
        }
        else
        {
            nlimbs = _nmod_vec_dot_bound_limbs(len, mod);
        }

        NMOD_VEC_DOT(s, i, len, vec1[i], vec2[len - 1 - i], mod, nlimbs);
    }

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

/* todo: better algorithms for large len */
int
_gr_nmod_vec_reciprocals(ulong * res, slong len, gr_ctx_t ctx)
{
    nmod_t mod = NMOD_CTX(ctx);
    slong k;
    ulong c2;

    if (len <= 1)
    {
        res[0] = (mod.n != 1);
        return GR_SUCCESS;
    }

    if (mod.n <= len || mod.n % 2 == 0)
        return GR_DOMAIN;

    res[0] = 1;
    res[1] = c2 = (mod.n - 1) / 2 + 1;;

    for (k = 3; k <= len; k += 2)
    {
        if (n_gcdinv(res + k - 1, k, mod.n) != 1)
            return GR_DOMAIN;
    }

    for (k = 4; k <= len; k += 2)
        res[k - 1] = nmod_mul(res[k / 2 - 1], c2, mod);

    return GR_SUCCESS;
}

int
_gr_nmod_poly_mullow(ulong * res,
    const ulong * poly1, slong len1,
    const ulong * poly2, slong len2, slong n, gr_ctx_t ctx)
{
    if (len1 + len2 - 1 == n)
    {
        if (len1 >= len2)
            _nmod_poly_mul(res, poly1, len1, poly2, len2, NMOD_CTX(ctx));
        else
            _nmod_poly_mul(res, poly2, len2, poly1, len1, NMOD_CTX(ctx));
    }
    else
    {
        if (len1 >= len2)
            _nmod_poly_mullow(res, poly1, len1, poly2, len2, n, NMOD_CTX(ctx));
        else
            _nmod_poly_mullow(res, poly2, len2, poly1, len1, n, NMOD_CTX(ctx));
    }

    return GR_SUCCESS;
}

/* fixme: duplicates _nmod_poly_divrem for error handling */
int
_gr_nmod_poly_divrem(mp_ptr Q, mp_ptr R, mp_srcptr A, slong lenA,
                                  mp_srcptr B, slong lenB, gr_ctx_t ctx)
{
    if (lenA <= 20 || lenB <= 8 || lenA - lenB <= 6 ||
            (NMOD_BITS(NMOD_CTX(ctx)) <= 61 && lenA <= 40) ||
            (NMOD_BITS(NMOD_CTX(ctx)) <= 29 && lenA <= 70))
    {
        mp_limb_t invB;
        int status;

        status = _gr_nmod_inv(&invB, &B[lenB - 1], ctx);
        if (status != GR_SUCCESS)
            return status;

        _nmod_poly_divrem_basecase_preinv1(Q, R, A, lenA, B, lenB, invB, NMOD_CTX(ctx));

        return status;
    }
    else
    {
#ifdef FLINT_HAVE_FFT_SMALL
        return _gr_poly_divrem_newton(Q, R, A, lenA, B, lenB, ctx);
#else
        if (NMOD_BITS(NMOD_CTX(ctx)) >= 16 && lenB >= 1024 && lenA <= 16384)
            return _gr_poly_divrem_divconquer(Q, R, A, lenA, B, lenB, 16, ctx);
        else
            return _gr_poly_divrem_newton(Q, R, A, lenA, B, lenB, ctx);
#endif
    }
}

/* basecase -> Newton cutoffs */
/* todo: better unbalanced tuning */

#ifdef FLINT_HAVE_FFT_SMALL

static const short inv_series_cutoff_tab[64] = {36, 37, 40, 44, 42, 46, 50, 81, 97,
  106, 133, 152, 292, 279, 279, 191, 191, 191, 191, 200, 191, 407, 388, 407, 407,
  407, 407, 407, 266, 292, 292, 292, 231, 242, 279, 279, 292, 292, 279, 292,
  292, 279, 292, 292, 292, 370, 370, 370, 370, 370, 370, 388, 407, 448,
  427, 407, 427, 427, 427, 427, 388, 370, 370, 353, };

#else

static const short inv_series_cutoff_tab[64] = {38, 36, 38, 36, 41, 48, 49, 54, 60,
  102, 112, 150, 165, 172, 210, 272, 339, 378, 385, 442, 468, 557, 596,
  621, 710, 746, 756, 978, 768, 679, 700, 696, 620, 619, 642, 766, 901,
  883, 924, 997, 979, 1028, 1101, 1094, 1152, 1169, 1279, 1311, 1284,
  1381, 1418, 1513, 1540, 1598, 1692, 1846, 1883, 1942, 1963, 1803,
  1788, 1861, 1881, 1920, };

#endif

void _nmod_poly_inv_series_basecase_preinv1(mp_ptr Qinv, mp_srcptr Q, slong Qlen, slong n, mp_limb_t q, nmod_t mod);

int
_gr_nmod_poly_inv_series_basecase(ulong * res,
    const ulong * f, slong flen, slong n, gr_ctx_t ctx)
{
    mp_limb_t q;

    q = f[0];
    if (q != 1)
    {
        if (n_gcdinv(&q, q, NMOD_CTX(ctx).n) != 1)
            return GR_DOMAIN;
    }

    _nmod_poly_inv_series_basecase_preinv1(res, f, flen, n, q, NMOD_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nmod_poly_inv_series(ulong * res,
    const ulong * f, slong flen, slong n, gr_ctx_t ctx)
{
    slong cutoff;

    flen = FLINT_MIN(flen, n);

    cutoff = inv_series_cutoff_tab[NMOD_BITS(NMOD_CTX(ctx)) - 1];

    if (flen < cutoff)
        return _gr_poly_inv_series_basecase(res, f, flen, n, ctx);
    else
        return _gr_poly_inv_series_newton(res, f, flen, n, cutoff, ctx);
}


void _nmod_poly_div_series_basecase_preinv1(mp_ptr Qinv, mp_srcptr P, slong Plen, mp_srcptr Q, slong Qlen, slong n, mp_limb_t q, nmod_t mod);

int
_gr_nmod_poly_div_series_basecase(ulong * res,
    const ulong * f, slong flen, const ulong * g, slong glen, slong n, gr_ctx_t ctx)
{
    mp_limb_t q;

    q = g[0];
    if (q != 1)
    {
        if (n_gcdinv(&q, q, NMOD_CTX(ctx).n) != 1)
            return GR_DOMAIN;
    }

    _nmod_poly_div_series_basecase_preinv1(res, f, flen, g, glen, n, q, NMOD_CTX(ctx));
    return GR_SUCCESS;
}

/* todo: unbalanced cutoffs */

#ifdef FLINT_HAVE_FFT_SMALL

static const short div_series_cutoff_tab[64] = {
  85, 81, 85, 89, 106, 145, 166, 182, 210, 337, 321, 321, 306, 321, 321, 306, 321, 321,
  321, 306, 626, 597, 597, 597, 597, 597, 597, 597, 370, 388, 388, 407, 321, 321, 321,
  321, 321, 337, 321, 337, 321, 321, 337, 321, 337, 427, 427, 448, 427, 427, 448, 448,
  470, 470, 470, 470, 470, 470, 493, 470, 427, 470, 448, 448, };

#else

static const short div_series_cutoff_tab[64] = {
  64, 64, 63, 65, 77, 84, 140, 158, 197, 208, 276, 243, 348, 358, 450, 560, 585, 692,
  749, 787, 864, 922, 914, 1239, 1454, 1424, 1449, 1327, 1115, 1174, 1250, 1282, 961,
  981, 1000, 1091, 1138, 1222, 1205, 1374, 1397, 1420, 1446, 1486, 1513, 1763, 1782,
  1890, 2047, 2071, 2075, 2271, 2192, 2307, 2407, 2324, 2308, 2958, 3001, 2816,
  2747, 2985, 2970, 2871, };

#endif

int
_gr_nmod_poly_div_series(ulong * res,
    const ulong * f, slong flen, const ulong * g, slong glen, slong n, gr_ctx_t ctx)
{
    slong cutoff;

    flen = FLINT_MIN(flen, n);
    glen = FLINT_MIN(glen, n);

    cutoff = div_series_cutoff_tab[NMOD_BITS(NMOD_CTX(ctx)) - 1];

    if (glen < cutoff)
        return _gr_poly_div_series_basecase(res, f, flen, g, glen, n, ctx);
    else
        return _gr_poly_div_series_newton(res, f, flen, g, glen, n, cutoff, ctx);
}


/* todo: unbalanced cutoffs */

#ifdef FLINT_HAVE_FFT_SMALL

static const short rsqrt_series_cutoff_tab[64] = {3, 22, 24, 42, 36, 40, 44,
    52, 89, 106, 139, 166, 174, 200, 191, 200, 191, 191, 191, 200, 306, 370,
    353, 337, 337, 370, 370, 370, 292, 292, 279, 292, 231, 231, 220, 279,
    292, 292, 279, 292, 292, 279, 292, 292, 292, 370, 370, 370, 370, 370,
    370, 388, 407, 448, 427, 407, 427, 427, 427, 427, 388, 370, 370, 353, };

#else

static const short rsqrt_series_cutoff_tab[64] = {6, 22, 22, 24, 27, 28, 28, 58,
  77, 96, 116, 160, 232, 270, 315, 387, 402, 472, 502, 580, 627, 760, 824, 940,
  988, 1018, 1155, 1182, 938, 932, 925, 1016, 836, 891, 915, 960, 1038, 1101,
  1203, 1236, 1255, 1311, 1386, 1422, 1489, 1592, 1624, 1879, 1828, 2055, 2227,
  2369, 2156, 2361, 2415, 2472, 2581, 2719, 2679, 2302, 2199, 2455, 2440, 2356, };

#endif

int
_gr_nmod_poly_rsqrt_series(ulong * res,
    const ulong * f, slong flen, slong n, gr_ctx_t ctx)
{
    slong cutoff;

    flen = FLINT_MIN(flen, n);

    cutoff = rsqrt_series_cutoff_tab[NMOD_BITS(NMOD_CTX(ctx)) - 1];

    if (flen < cutoff)
        return _gr_poly_rsqrt_series_basecase(res, f, flen, n, ctx);
    else
        return _gr_poly_rsqrt_series_newton(res, f, flen, n, cutoff, ctx);
}

#ifdef FLINT_HAVE_FFT_SMALL

static const short sqrt_series_cutoff_tab[] = { 32767, 1353, 1353, 919,
    1289, 1228, 1491, 1228, 1491, 1289, 1420, 1420, 1228, 1289, 1420,
    1353, 1289, 1289, 1353, 4345, 4345, 2423, 2544, 4345, 4562, 2423,
    4345, 2671, 1725, 1811, 1725, 3577, 1643, 1643, 1725, 1725, 1725,
    3407, 3407, 1725, 1643, 1811, 3407, 3407, 4562, 3755, 3755, 3577,
    3577, 4562, 3755, 3942, 3755, 3755, 4562, 4562, 4562, 4562, 3755,
    3245, 3245, 2944, 3091, 2423, };

#else

static const short sqrt_series_cutoff_tab[] = { 32767, 632, 732, 928, 1443,
  1731, 2364, 2490, 2893, 3173, 5316, 5412, 5727, 6123, 6613, 7290, 7572,
  8023, 9114, 9105, 8656, 10645, 11290, 13223, 11507, 15489, 12328, 9338,
  9517, 9795, 9596, 13162, 10168, 8013, 9949, 10654, 9932, 12222, 11287,
  11623, 11971, 12577, 12207, 13886, 14160, 12200, 13207, 15943, 15320,
  14290, 15933, 15463, 14281, 15457, 15302, 17929, 18106, 17058, 14844,
  17740, 17916, 18640, 18093, 18638, };

#endif

/* todo: unbalanced cutoffs */
int
_gr_nmod_poly_sqrt_series(ulong * res,
    const ulong * f, slong flen, slong n, gr_ctx_t ctx)
{
    slong cutoff;

    flen = FLINT_MIN(flen, n);

    cutoff = sqrt_series_cutoff_tab[NMOD_BITS(NMOD_CTX(ctx)) - 1];

    if (flen < cutoff)
        return _gr_poly_sqrt_series_basecase(res, f, flen, n, ctx);
    else
        return _gr_poly_sqrt_series_newton(res, f, flen, n, cutoff, ctx);
}

#ifdef FLINT_HAVE_FFT_SMALL

static const short exp_series_mul_cutoff_tab[] = { 470, 470, 470, 470,
    470, 470, 470, 470, 470, 266, 266, 231, 279, 292, 292, 266, 266,
    254, 266, 266, 448, 470, 470, 448, 448, 470, 493, 470, 337, 321,
    427, 321, 292, 292, 292, 292, 292, 292, 279, 292, 292, 292, 279,
    279, 292, 448, 448, 427, 448, 448, 448, 448, 448, 448, 427, 448,
    427, 448, 448, 448, 353, 407, 407, 353, };

#else

static const short exp_series_mul_cutoff_tab[] = {
    4, 4, 6, 12, 18, 38, 56, 57, 97, 121, 144, 145, 187, 228, 216, 286,
    315, 426, 505, 547, 620, 631, 649, 781, 931, 810, 950, 1066, 838, 714,
    880, 931, 698, 598, 590, 620, 731, 782, 780, 848, 976, 1005, 978, 1006,
    1046, 1040, 1203, 1296, 1342, 1327, 1389, 1555, 1568, 1644, 1676, 1725,
    1692, 1793, 2014, 1904, 2021, 2086, 2030, 2198,
};

#endif

#ifdef FLINT_HAVE_FFT_SMALL

static const short exp_series_newton_cutoff_tab[64] = { 470, 470, 470,
    470, 470, 470, 517, 597, 689, 759, 689, 723, 626, 723, 723, 723,
    723, 723, 723, 723, 1289, 1491, 1353, 1353, 1353, 1353, 1353, 1353,
    964, 964, 1115, 1115, 835, 876, 876, 919, 876, 919, 919, 876, 919,
    876, 835, 919, 919, 1228, 1170, 1289, 1228, 1228, 1353, 1353,
    1353, 1228, 1228, 1228, 1289, 1353, 1228, 1170, 1115, 1115, 1170,
    1353, };

#else

static const short exp_series_newton_cutoff_tab[64] = {
    4, 4, 6, 12, 18, 38, 68, 132, 258, 522, 1033, 1288, 1322, 1494, 1775,
    2074, 2404, 2534, 2887, 3017, 3343, 3234, 3550, 4263, 3813, 4534,
    5398, 4048, 3381, 3607, 3811, 3724, 2962, 3095, 3268, 3253, 3469, 3664,
    4108, 4219, 4043, 4402, 4434, 4377, 4583, 4802, 5030, 5430, 5366, 5407,
    5341, 6363, 6415, 6343, 6395, 6444, 6433, 6358, 6939, 6794, 6782, 6837,
    6758, 6680,
};

#endif

int
_gr_nmod_poly_exp_series(ulong * res,
    const ulong * f, slong flen, slong n, gr_ctx_t ctx)
{
    slong cutoff1, cutoff2;

    flen = FLINT_MIN(flen, n);

    cutoff1 = exp_series_mul_cutoff_tab[NMOD_BITS(NMOD_CTX(ctx)) - 1];

    if (flen < cutoff1)
        return _gr_poly_exp_series_basecase(res, f, flen, n, ctx);

    cutoff2 = exp_series_newton_cutoff_tab[NMOD_BITS(NMOD_CTX(ctx)) - 1];

    if (flen < cutoff2)
        return _gr_poly_exp_series_basecase_mul(res, f, flen, n, ctx);

    return _gr_poly_exp_series_newton(res, NULL, f, flen, n, cutoff2, ctx);
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

int
_gr_nmod_mat_mul(gr_mat_t res, const gr_mat_t x, const gr_mat_t y, gr_ctx_t ctx)
{
    nmod_mat_t R, X, Y;
    nmod_mat_struct *XX, *YY;

    R->entries = res->entries;
    R->rows = (mp_ptr *) res->rows;
    R->r = res->r;
    R->c = res->c;
    R->mod = NMOD_CTX(ctx);

    if (res == x)
    {
        XX = R;
    }
    else
    {
        X->entries = x->entries;
        X->rows = (mp_ptr *) x->rows;
        X->r = x->r;
        X->c = x->c;
        X->mod = NMOD_CTX(ctx);
        XX = X;
    }

    if (res == y)
    {
        YY = R;
    }
    else if (x == y)
    {
        YY = XX;
    }
    else
    {
        Y->entries = y->entries;
        Y->rows = (mp_ptr *) y->rows;
        Y->r = y->r;
        Y->c = y->c;
        Y->mod = NMOD_CTX(ctx);
        YY = Y;
    }

    nmod_mat_mul(R, XX, YY);

    return GR_SUCCESS;
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
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _gr_nmod_set_shallow},
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
    {GR_METHOD_ADD_UI,          (gr_funcptr) _gr_nmod_add_ui},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_nmod_sub},
    {GR_METHOD_SUB_SI,          (gr_funcptr) _gr_nmod_sub_si},
    {GR_METHOD_SUB_UI,          (gr_funcptr) _gr_nmod_sub_ui},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_nmod_mul},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_nmod_mul_si},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_nmod_mul_ui},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_nmod_mul_fmpz},
    {GR_METHOD_ADDMUL,          (gr_funcptr) _gr_nmod_addmul},
    {GR_METHOD_SUBMUL,          (gr_funcptr) _gr_nmod_submul},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) _gr_nmod_mul_two},
    {GR_METHOD_MUL_2EXP_SI,     (gr_funcptr) _gr_nmod_mul_2exp_si},
    {GR_METHOD_SQR,             (gr_funcptr) _gr_nmod_sqr},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_nmod_div},
    {GR_METHOD_DIV_SI,          (gr_funcptr) _gr_nmod_div_si},
    {GR_METHOD_DIV_UI,          (gr_funcptr) _gr_nmod_div_ui},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) _gr_nmod_div_fmpz},
    {GR_METHOD_DIV_NONUNIQUE,   (gr_funcptr) _gr_nmod_div_nonunique},
    {GR_METHOD_DIVIDES,         (gr_funcptr) _gr_nmod_divides},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_nmod_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) _gr_nmod_inv},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_nmod_sqrt},
    {GR_METHOD_VEC_INIT,        (gr_funcptr) _gr_nmod_vec_init},
    {GR_METHOD_VEC_CLEAR,       (gr_funcptr) _gr_nmod_vec_clear},
    {GR_METHOD_VEC_SET,         (gr_funcptr) _gr_nmod_vec_set},
    {GR_METHOD_VEC_NORMALISE,   (gr_funcptr) _gr_nmod_vec_normalise},
    {GR_METHOD_VEC_NORMALISE_WEAK,   (gr_funcptr) _gr_nmod_vec_normalise_weak},
    {GR_METHOD_VEC_NEG,         (gr_funcptr) _gr_nmod_vec_neg},
    {GR_METHOD_VEC_ADD,         (gr_funcptr) _gr_nmod_vec_add},
    {GR_METHOD_VEC_SUB,         (gr_funcptr) _gr_nmod_vec_sub},
    {GR_METHOD_VEC_MUL_SCALAR,      (gr_funcptr) _gr_nmod_vec_mul_scalar},
    {GR_METHOD_VEC_MUL_SCALAR_SI,   (gr_funcptr) _gr_nmod_vec_mul_scalar_si},
    {GR_METHOD_VEC_MUL_SCALAR_UI,   (gr_funcptr) _gr_nmod_vec_mul_scalar_ui},
    {GR_METHOD_VEC_MUL_SCALAR_FMPZ, (gr_funcptr) _gr_nmod_vec_mul_scalar_fmpz},
    {GR_METHOD_VEC_MUL_SCALAR_2EXP_SI,   (gr_funcptr) _gr_nmod_vec_mul_scalar_2exp_si},
    {GR_METHOD_VEC_ADDMUL_SCALAR,        (gr_funcptr) _gr_nmod_vec_addmul_scalar},
    {GR_METHOD_VEC_ADDMUL_SCALAR_SI,     (gr_funcptr) _gr_nmod_vec_addmul_scalar_si},
    {GR_METHOD_VEC_SUBMUL_SCALAR,        (gr_funcptr) _gr_nmod_vec_submul_scalar},
    {GR_METHOD_VEC_SUBMUL_SCALAR_SI,     (gr_funcptr) _gr_nmod_vec_submul_scalar_si},
    {GR_METHOD_VEC_SUB,         (gr_funcptr) _gr_nmod_vec_sub},
    {GR_METHOD_VEC_SUM,         (gr_funcptr) _gr_nmod_vec_sum},
    {GR_METHOD_VEC_PRODUCT,     (gr_funcptr) _gr_nmod_vec_product},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) __gr_nmod_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) __gr_nmod_vec_dot_rev},
    {GR_METHOD_VEC_RECIPROCALS, (gr_funcptr) _gr_nmod_vec_reciprocals},
    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _gr_nmod_poly_mullow},
    {GR_METHOD_POLY_DIVREM,     (gr_funcptr) _gr_nmod_poly_divrem},
    {GR_METHOD_POLY_INV_SERIES, (gr_funcptr) _gr_nmod_poly_inv_series},
    {GR_METHOD_POLY_INV_SERIES_BASECASE, (gr_funcptr) _gr_nmod_poly_inv_series_basecase},
    {GR_METHOD_POLY_DIV_SERIES, (gr_funcptr) _gr_nmod_poly_div_series},
    {GR_METHOD_POLY_DIV_SERIES_BASECASE, (gr_funcptr) _gr_nmod_poly_div_series_basecase},
    {GR_METHOD_POLY_RSQRT_SERIES, (gr_funcptr) _gr_nmod_poly_rsqrt_series},
    {GR_METHOD_POLY_SQRT_SERIES,  (gr_funcptr) _gr_nmod_poly_sqrt_series},
    {GR_METHOD_POLY_EXP_SERIES,  (gr_funcptr) _gr_nmod_poly_exp_series},
    {GR_METHOD_POLY_ROOTS,      (gr_funcptr) _gr_nmod_roots_gr_poly},
    {GR_METHOD_MAT_MUL,         (gr_funcptr) _gr_nmod_mat_mul},

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

void
_gr_ctx_init_nmod(gr_ctx_t ctx, void * nmod_t_ref)
{
    ctx->which_ring = GR_CTX_NMOD;
    ctx->sizeof_elem = sizeof(ulong);
    ctx->size_limit = WORD_MAX;

    *NMOD_CTX_REF(ctx) = ((nmod_t *) nmod_t_ref)[0];
    ctx->methods = __gr_nmod_methods;

    if (!__gr_nmod_methods_initialized)
    {
        gr_method_tab_init(__gr_nmod_methods, __gr_nmod_methods_input);
        __gr_nmod_methods_initialized = 1;
    }
}
