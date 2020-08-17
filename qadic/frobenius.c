/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "ulong_extras.h"
#include "qadic.h"

/*
    Assumes that \code{len1} and \code{len2} are positive but at 
    most~$d$, and also that \code{len1} is at least $6$.

    The latter assumption guarantees that $\ceil{n/B} \geq 2$, 
    i.e.\ $n \geq 2B$ so $n \geq 2 \ceil{\sqrt{n}}$.
 */

static void 
_fmpz_mod_poly_compose_smod_rectangular(fmpz *rop, 
                           const fmpz *op1, slong len1, 
                           const fmpz *op2, slong len2, 
                           const fmpz *a, const slong *j, slong lena, 
                           const fmpz_t p)
{
    const slong d = j[lena - 1];

    if (len2 == 1)
    {
        _fmpz_mod_poly_evaluate_fmpz(rop, op1, len1, op2, p);
        _fmpz_vec_zero(rop + 1, d - 1);
    }
    else
    {
        const slong B = n_sqrt(len1);
        slong i, k;
        fmpz *pows, *t;

        pows = _fmpz_vec_init((B + 2) * d);
        t    = _fmpz_vec_init(2 * d - 1);

        fmpz_one(pows + 0 * d + 0);
        _fmpz_vec_set(pows + 1 * d, op2, len2);
        for (i = 2; i <= B; i++)
        {
            _fmpz_poly_mul(pows + i * d, pows + (i - 1) * d, d, op2, len2);
            _fmpz_poly_reduce(pows + i * d, d + len2 - 1, a, j, lena);
            _fmpz_vec_scalar_mod_fmpz(pows + i * d, pows + i * d, d, p);
        }

        _fmpz_vec_zero(rop, d);

        for (i = (len1 + B - 1) / B - 1; i >= 0; i--)
        {
            _fmpz_poly_mul(t, rop, d, pows + B * d, d);
            _fmpz_poly_reduce(t, 2 * d - 1, a, j, lena);

            _fmpz_vec_set(rop, t, d);
            fmpz_add(rop + 0, rop + 0, op1 + i*B);
            for (k = FLINT_MIN(B, len1 - i*B) - 1; k > 0; k--)
            {
                _fmpz_vec_scalar_addmul_fmpz(rop, pows + k * d, d, op1 + (i*B + k));
            }

            _fmpz_vec_scalar_mod_fmpz(rop, rop, d, p);
        }

        _fmpz_vec_clear(pows, (B + 2) * d);
        _fmpz_vec_clear(t, 2 * d - 1);
    }
}

static void 
_fmpz_mod_poly_compose_smod_horner(fmpz *rop, 
                           const fmpz *op1, slong len1, 
                           const fmpz *op2, slong len2, 
                           const fmpz *a, const slong *j, slong lena, 
                           const fmpz_t p)
{
    const slong d = j[lena - 1];

    if (len1 == 1)
    {
        fmpz_set(rop, op1);
        _fmpz_vec_zero(rop + 1, d - 1);
    }
    else if (len2 == 1)
    {
        _fmpz_mod_poly_evaluate_fmpz(rop, op1, len1, op2, p);
        _fmpz_vec_zero(rop + 1, d - 1);
    }
    else
    {
        slong i;
        fmpz *t;

        t = _fmpz_vec_init(2*d - 1);

        _fmpz_vec_zero(rop, d);

        for (i = len1 - 1; i >= 0; i--)
        {
            _fmpz_poly_mul(t, rop, d, op2, len2);
            _fmpz_poly_reduce(t, d + len2 - 1, a, j, lena);
            _fmpz_poly_add(rop, t, d, op1 + i, 1);
            _fmpz_vec_scalar_mod_fmpz(rop, rop, d, p);
        }

        _fmpz_vec_clear(t, 2*d - 1);
    }
}

/* 
    Computes the composition $f(g(X))$ modulo the sparse polynomial 
    given by the data \code{(a, j, lena)}, which is assumed to be 
    of degree~$d \geq 2$.

    Sets the vector \code{(rop, d)}.

    Assumes that \code{len1} and \code{len2} are positive but at 
    most~$d$.

    Does not support aliasing.
 */

void 
_fmpz_mod_poly_compose_smod(fmpz *rop, 
                           const fmpz *op1, slong len1, 
                           const fmpz *op2, slong len2, 
                           const fmpz *a, const slong *j, slong lena, 
                           const fmpz_t p)
{
    if (len1 < 6)
    {
        _fmpz_mod_poly_compose_smod_horner(rop, op1, len1, op2, len2, a, j, lena, p);
    }
    else
    {
        _fmpz_mod_poly_compose_smod_rectangular(rop, op1, len1, op2, len2, a, j, lena, p);
    }
}

void _qadic_frobenius_a(fmpz *rop, slong exp, 
                    const fmpz *a, const slong *j, slong lena, 
                    const fmpz_t p, slong N)
{
    const slong d = j[lena - 1];

    slong *e, i, n;
    fmpz *pow, *f1, *f2, *inv, *s, *t;

    n = FLINT_CLOG2(N) + 1;

    e = flint_malloc(n * sizeof(slong));
    for (e[i = 0] = N; e[i] > 1; i++)
        e[i + 1] = (e[i] + 1) / 2;

    pow = _fmpz_vec_init(n);
    f1  = _fmpz_vec_init(d + 1);
    f2  = _fmpz_vec_init(d);
    inv = _fmpz_vec_init(2*d - 1);
    s   = _fmpz_vec_init(2*d - 1);
    t   = _fmpz_vec_init(2*d - 1);

    /* Compute powers of p */
    {
        fmpz_one(t);
        fmpz_set(pow + i, p);
    }
    for (i--; i >= 1; i--)
    {
        if (e[i] & WORD(1))
        {
            fmpz_mul(pow + i, t, pow + (i + 1));
            fmpz_mul(t, t, t);
        }
        else
        {
            fmpz_mul(t, t, pow + (i + 1));
            fmpz_mul(pow + i, pow + (i + 1), pow + (i + 1));
        }
    }
    {
        if (e[i] & WORD(1))
            fmpz_mul(pow + i, t, pow + (i + 1));
        else
            fmpz_mul(pow + i, pow + (i + 1), pow + (i + 1));
    }

    /* Dense representation of f and f' */
    {
        slong k;

        for (k = 0; k < lena; k++)
            fmpz_set(f1 + j[k], a + k);
        for (k = 1; k < lena; k++)
            fmpz_mul_ui(f2 + (j[k] - 1), a + k, j[k]);
    }

    /* Run Newton iteration */
    i = n - 1;
    {
        fmpz op[2] = {WORD(0), WORD(1)};

        fmpz_pow_ui(t, p, exp);
        _qadic_pow(rop, op, 2, t, a, j, lena, pow + i);
        _fmpz_mod_poly_compose_smod(t, f2, d, rop, d, a, j, lena, pow + i);
        _qadic_inv(inv, t, d, a, j, lena, p, 1);
    }
    for (i--; i >= 0; i--)
    {
        _fmpz_mod_poly_compose_smod(s, f1, d + 1, rop, d, a, j, lena, pow + i);
        _fmpz_mod_poly_mul(t, s, d, inv, d, pow + i);
        _fmpz_mod_poly_reduce(t, 2*d - 1, a, j, lena, pow + i);
        _fmpz_mod_poly_sub(rop, rop, d, t, d, pow + i);

        if (i > 0)
        {
            _fmpz_mod_poly_compose_smod(s, f2, d, rop, d, a, j, lena, pow + i);
            _fmpz_mod_poly_mul(t, inv, d, s, d, pow + i);
            _fmpz_mod_poly_reduce(t, 2*d - 1, a, j, lena, pow + i);
            fmpz_sub_ui(t, t, 2);
            if (fmpz_sgn(t) < 0)
                fmpz_add(t, t, pow + i);
            _fmpz_mod_poly_neg(t, t, d, pow + i);
            _fmpz_mod_poly_mul(s, inv, d, t, d, pow + i);
            _fmpz_mod_poly_reduce(s, 2*d - 1, a, j, lena, pow + i);

            /* SWAP(inv, s).  Requires the arrays to be of the same size. */
            {
                fmpz *__t;

                __t = inv;
                inv = s;
                s   = __t;
            }
        }
    }

    _fmpz_vec_clear(pow, n);
    _fmpz_vec_clear(f1, d + 1);
    _fmpz_vec_clear(f2, d);
    _fmpz_vec_clear(inv, 2*d - 1);
    _fmpz_vec_clear(s, 2*d - 1);
    _fmpz_vec_clear(t, 2*d - 1);
    flint_free(e);
}

/*
    Sets (rop, 2d-1) to the image of (op, len) under the Frobenius operator 
    raised to the e-th power.
 */

void _qadic_frobenius(fmpz *rop, const fmpz *op, slong len, slong e, 
                  const fmpz *a, const slong *j, slong lena, 
                  const fmpz_t p, slong N)
{
    const slong d = j[lena - 1];

    if (len == 1)  /* op is in Zp, not just Zq */
    {
        _fmpz_vec_set(rop, op, len);
        _fmpz_vec_zero(rop + len, (2*d - 1)  - len);
    }
    else if (N == 1)
    {
        fmpz_t t;

        fmpz_init(t);
        fmpz_pow_ui(t, p, e);
        _qadic_pow(rop, op, len, t, a, j, lena, p);
        fmpz_clear(t);
    }
    else
    {
        fmpz *t;
        fmpz_t pow;

        t = _fmpz_vec_init(2*d - 1);
        fmpz_init(pow);
        fmpz_pow_ui(pow, p, N);

        _qadic_frobenius_a(t, e, a, j, lena, p, N);

        _fmpz_mod_poly_compose_smod(rop, op, len, t, d, a, j, lena, pow);
        _fmpz_vec_zero(rop + d, d - 1);

        _fmpz_vec_clear(t, 2*d - 1);
        fmpz_clear(pow);
    }
}

void qadic_frobenius(qadic_t rop, const qadic_t op, slong e, const qadic_ctx_t ctx)
{
    const slong N = qadic_prec(rop);
    const slong d = qadic_ctx_degree(ctx);

    e = e % d;
    if (e < 0)
        e += d;

    if (qadic_is_zero(op) || op->val >= N)
    {
        qadic_zero(rop);
    }
    else if (e == 0)
    {
        padic_poly_set(rop, op, &ctx->pctx);
    }
    else
    {
        fmpz *t;

        if (rop == op)
        {
            t = _fmpz_vec_init(2 * d - 1);
        }
        else
        {
            padic_poly_fit_length(rop, 2 * d - 1);
            t = rop->coeffs;
        }

        _qadic_frobenius(t, op->coeffs, op->length, e, 
                     ctx->a, ctx->j, ctx->len, (&ctx->pctx)->p, N - op->val);

        if (rop == op)
        {
            _fmpz_vec_clear(rop->coeffs, rop->alloc);
            rop->coeffs = t;
            rop->alloc  = 2 * d - 1;
            rop->length = d;
        }
        else
        {
            rop->val = op->val;
            _padic_poly_set_length(rop, d);
        }
        _padic_poly_normalise(rop);
    }
}

