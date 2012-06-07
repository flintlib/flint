/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Sebastian Pancratz
 
******************************************************************************/

#include "fmpz_mod_poly.h"
#include "qadic.h"

/*
    Returns whether the non-zero element \code{(op, len)} 
    has a square root in $\mathbf{F}_q$, and if so sets 
    sets \code{(rop, 2 * d - 1)} to its square root.

    Assumes that $p$ is an odd prime.

    Assumes that $d \geq 1$.

    Does not support aliasing.
 */

static int 
_fmpz_mod_poly_sqrtmod_p(fmpz *rop, const fmpz *op, long len, 
                         const fmpz *a, const long *j, long lena, 
                         const fmpz_t p)
{
    const long d = j[lena - 1];
    int ans;

    /*
        In the case when $q \equiv 3 \pmod{4}$, if $x$ is a square 
        it is easy to show that one of its square roots is equal to 
        $x^{(q+1)/4}$.

        To find out whether $x$ is a square in the first place, we 
        would like to compute $x^{(q-1)/2}$ and check that it is 
        equal to $1$.

        To avoid unnecessary recomputation of powers of $x$, we 
        compute $x^{(q-3)/4}$, multiply this by $x$ to obtain 
        the potential square root $x^{(q+1)/4}$, and then combine 
        these two powers to find $x^{(q-1)/2}$.
     */
    if (fmpz_fdiv_ui(p, 4) == 3 && (d & 1L))
    {
        fmpz_t z;
        fmpz *v, *w;

        fmpz_init(z);
        v = _fmpz_vec_init(4 * d - 2);
        w = v + (2 * d - 1);

        fmpz_pow_ui(z, p, d);
        fmpz_sub_ui(z, z, 3);
        fmpz_fdiv_q_2exp(z, z, 2);

        _qadic_pow(v, op, len, z, a, j, lena, p);

        _fmpz_mod_poly_mul(rop, v, d, op, len, p);
        _fmpz_vec_zero(rop + d + len - 1, d - len);
        _fmpz_mod_poly_reduce(rop, d + len - 1, a, j, lena, p);

        _fmpz_mod_poly_mul(w, rop, d, v, d, p);
        _fmpz_mod_poly_reduce(w, 2 * d - 1, a, j, lena, p);
        ans = fmpz_is_one(w + 0);

        fmpz_clear(z);
        _fmpz_vec_clear(v, 4 * d - 2);
    }

    /*
        In this case, $q \equiv 1 \pmod{4}$, the procedure for computing 
        the square root is more involved and we test up front whether 
        the element is indeed a square by computing $x^{(q-1)/2}$.

        TODO:  Find a better way to integrate the check for square-ness 
        into the computation of a potential square root.
     */
    else
    {
        long i, s;

        fmpz_t t, pm1, qm1, z;

        fmpz *b, *g, *bpow, *gpow, *w;

        fmpz_init(t);
        fmpz_init(pm1);
        fmpz_init(qm1);
        fmpz_init(z);

        fmpz_sub_ui(pm1, p, 1);
        fmpz_pow_ui(qm1, p, d);
        fmpz_sub_ui(qm1, qm1, 1);

        b    = _fmpz_vec_init(2 * d - 1);
        g    = _fmpz_vec_init(2 * d - 1);
        bpow = _fmpz_vec_init(2 * d - 1);
        gpow = _fmpz_vec_init(2 * d - 1);
        w    = _fmpz_vec_init(2 * d - 1);

        /* Check whether op is a square, i.e. op^{(q-1}/2} == 1 */
        fmpz_fdiv_q_2exp(z, qm1, 1);
        _qadic_pow(w, op, len, z, a, j, lena, p);
        ans = fmpz_is_one(w);
        if (!ans)
            goto exit;

        /*
            Find a non-residue g; 
            uses stack-based depth first traversal starting from [0,...,0,1]
         */
        fmpz_fdiv_q_2exp(z, qm1, 1);
        {
            fmpz_one(g + (d - 1));
            i = d;

            while (i > 0)
            {
                if (i == d)
                {
                    /* Consider this element, g^{(q-1)/2} == -1 ? */
                    _qadic_pow(w, g, d, z, a, j, lena, p);
                    if (fmpz_equal(w + 0, pm1))
                        break;

                    /* Backtrace, find the next element */
                    for (i--; i >= 0 && fmpz_equal(g + i, pm1); i--) ;
                    if (i >= 0)
                    {
                        fmpz_add_ui(g + i, g + i, 1);
                        i++;
                    }
                }
                else
                {
                    _fmpz_vec_zero(g + i, d - i);
                    i = d;
                }
            }
        }

        /* Write q - 1 = 2^s t */
        for (s = 0, fmpz_set(t, qm1); fmpz_is_even(t); s++)
            fmpz_fdiv_q_2exp(t, t, 1);

        /* Set g = g^t */
        _qadic_pow(w, g, d, t, a, j, lena, p);
        _fmpz_vec_set(g, w, d);

        /* Set rop = op^{(t+1)/2} */
        fmpz_add_ui(z, t, 1);
        fmpz_fdiv_q_2exp(z, z, 1);
        _qadic_pow(rop, op, len, z, a, j, lena, p);

        /* Set b = op^t */
        _qadic_pow(b, op, len, t, a, j, lena, p);

        while (!_fmpz_poly_is_one(b, d))
        {
            long k;

            _fmpz_vec_set(bpow, b, d);
            for (k = 1; (k < s) && !_fmpz_poly_is_one(bpow, d); k++)
            {
                _fmpz_poly_sqr(w, bpow, d);
                _fmpz_mod_poly_reduce(w, 2 * d - 1, a, j, lena, p);
                _fmpz_vec_scalar_mod_fmpz(bpow, w, d, p);
            }

            _fmpz_vec_set(gpow, g, d);
            for (i = 1; i < s - k; i++)
            {
                _fmpz_poly_sqr(w, gpow, d);
                _fmpz_mod_poly_reduce(w, 2 * d - 1, a, j, lena, p);
                _fmpz_vec_scalar_mod_fmpz(gpow, w, d, p);
            }

            _fmpz_poly_mul(w, rop, d, gpow, d);
            _fmpz_mod_poly_reduce(w, 2 * d - 1, a, j, lena, p);
            _fmpz_vec_scalar_mod_fmpz(rop, w, d, p);

            _fmpz_poly_sqr(w, gpow, d);
            _fmpz_mod_poly_reduce(w, 2 * d - 1, a, j, lena, p);
            _fmpz_vec_scalar_mod_fmpz(gpow, w, d, p);

            _fmpz_poly_mul(w, b, d, gpow, d);
            _fmpz_mod_poly_reduce(w, 2 * d - 1, a, j, lena, p);
            _fmpz_vec_scalar_mod_fmpz(b, w, d, p);

            s = k;
        }

      exit:

        fmpz_clear(t);
        fmpz_clear(pm1);
        fmpz_clear(qm1);
        fmpz_clear(z);
        _fmpz_vec_clear(b,    2 * d - 1);
        _fmpz_vec_clear(g,    2 * d - 1);
        _fmpz_vec_clear(bpow, 2 * d - 1);
        _fmpz_vec_clear(gpow, 2 * d - 1);
        _fmpz_vec_clear(w,    2 * d - 1);
    }
    return ans;
}

/*
    Sets \code{(rop, 2 * d - 1)} to the square root of 
    \code{(op, len)}.

    Note that the group of units of $\mathbf{F}_q$ is cyclic 
    of order $q - 1$, which is odd if $p = 2$.  In this case, 
    computing square roots amounts to exponentiating by 
    $q / 2$.

    Note that $x^{q-1} = 1$ for every non-zero $x$, so $x^q = x$ 
    for every $x$ and $u = x^{q/2}$ satisfies $u^2 = x$.

    In particular, every element is a square, and every element 
    has a unique square root.

    Assumes that $d \geq 1$.
 */

static void 
_fmpz_mod_poly_sqrtmod_2(fmpz *rop, const fmpz *op, long len, 
                         const fmpz *a, const long *j, long lena)
{
    const fmpz_t p = {2L};
    const long d   = j[lena - 1];

    fmpz_t z;

    fmpz_init(z);
    fmpz_setbit(z, d - 1);
    _qadic_pow(rop, op, len, z, a, j, lena, p);
    fmpz_clear(z);
}

/*
    Returns whether \code{(op, len)} is a square, and if so 
    sets \code{(rop, 2 * d - 1)} to a square root, reduced 
    modulo $p^N$.

    Assumes that \code{(op, len)} is a unit.
 */

static int 
_qadic_sqrt_p(fmpz *rop, const fmpz *op, long len, 
              const fmpz *a, const long *j, long lena, 
              const fmpz_t p, long N)
{
    const long d = j[lena - 1];
    int ans;

    if (N == 1)
    {
        ans = _fmpz_mod_poly_sqrtmod_p(rop, op, len, a, j, lena, p);
        return ans;
    }
    else
    {
        long *e, i, k, n;
        fmpz *pow, *u;
        fmpz *r, *s, *t;

        n = FLINT_CLOG2(N) + 1;

        /* Compute sequence of exponents */
        e = flint_malloc(n * sizeof(long));
        for (e[i = 0] = N; e[i] > 1; i++)
            e[i + 1] = (e[i] + 1) / 2;

        pow = _fmpz_vec_init(n);
        u   = _fmpz_vec_init(len * n);
        r   = _fmpz_vec_init(2 * d - 1);
        s   = _fmpz_vec_init(2 * d - 1);
        t   = _fmpz_vec_init(2 * d - 1);

        /* Compute powers of p */
        {
            fmpz_one(t);
            fmpz_set(pow + i, p);
        }
        for (i--; i >= 1; i--)
        {
            if (e[i] & 1L)
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
            if (e[i] & 1L)
                fmpz_mul(pow + i, t, pow + (i + 1));
            else
                fmpz_mul(pow + i, pow + (i + 1), pow + (i + 1));
        }

        /* Compute reduced units */
        {
            _fmpz_vec_scalar_mod_fmpz(u + 0 * len, op, len, pow + 0);
        }
        for (i = 1; i < n; i++)
        {
            _fmpz_vec_scalar_mod_fmpz(u + i * len, u + (i - 1) * len, len, pow + i);
        }

        /* Run Newton iteration */
        i = n - 1;
        {
            ans = _fmpz_mod_poly_sqrtmod_p(t, u + i * len, len, a, j, lena, p);
            if (!ans)
                goto exit;

            /* Dense copy of f, used for inversion */
            for (k = 0; k < lena; k++)
                fmpz_set(s + j[k], a + k);
            _fmpz_mod_poly_invmod(rop, t, d, s, d + 1, p);
        }
        for (i--; i >= 1; i--)  /* z := z - z (a z^2 - 1) / 2 */
        {
            _fmpz_poly_sqr(s, rop, d);
            _fmpz_poly_reduce(s, 2 * d - 1, a, j, lena);
            _fmpz_poly_mul(t, s, d, u + i * len, len);
            _fmpz_poly_reduce(t, d + len - 1, a, j, lena);
            fmpz_sub_ui(t, t, 1);

            for (k = 0; k < d; k++)
            {
                if (fmpz_is_odd(t + k))
                    fmpz_add(t + k, t + k, pow + i);
                fmpz_fdiv_q_2exp(t + k, t + k, 1);
            }

            _fmpz_poly_mul(s, t, d, rop, d);
            _fmpz_poly_reduce(s, 2 * d - 1, a, j, lena);
            _fmpz_poly_sub(rop, rop, d, s, d);
            _fmpz_vec_scalar_mod_fmpz(rop, rop, d, pow + i);
        }
        {
            _fmpz_poly_mul(s, rop, d, u + 1 * len, len);
            _fmpz_poly_reduce(s, d + len - 1, a, j, lena);
            _fmpz_poly_sqr(t, s, d);
            _fmpz_poly_reduce(t, 2 * d - 1, a, j, lena);
            _fmpz_poly_sub(t, u + 0 * len, len, t, d);

            for (k = 0; k < d; k++)
            {
                if (fmpz_is_odd(t + k))
                    fmpz_add(t + k, t + k, pow + 0);
                fmpz_fdiv_q_2exp(t + k, t + k, 1);
            }

            _fmpz_poly_mul(r, rop, d, t, d);
            _fmpz_poly_reduce(r, 2 * d - 1, a, j, lena);
            _fmpz_poly_add(rop, r, d, s, d);
            _fmpz_vec_scalar_mod_fmpz(rop, rop, d, pow + 0);
        }

      exit:

        _fmpz_vec_clear(pow, n);
        _fmpz_vec_clear(u, len * n);
        _fmpz_vec_clear(r, 2 * d - 1);
        _fmpz_vec_clear(s, 2 * d - 1);
        _fmpz_vec_clear(t, 2 * d - 1);
        flint_free(e);

        return ans;
    }
}

/*
    Returns whether \code{(op, len)} is a square, and if so 
    sets \code{(rop, 2 * d - 1)} to a square root, reduced 
    modulo $2^N$.

    Assumes that \code{(op, len)} is a unit.
 */

int _qadic_sqrt(fmpz *rop, const fmpz *op, long len, 
                const fmpz *a, const long *j, long lena, 
                const fmpz_t p, long N)
{
    if (*p == 2L)
    {
        printf("Exception (_qadic_sqrt).  Case p = 2.\n");
        abort();
    }
    else
    {
        return _qadic_sqrt_p(rop, op, len, a, j, lena, p, N);
    }
}

int qadic_sqrt(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
{
    const fmpz *p = (&ctx->pctx)->p;
    const long d  = qadic_ctx_degree(ctx);
    const long N  = (&ctx->pctx)->N;

    fmpz *t;
    int ans;

    if (qadic_is_zero(op))
    {
        qadic_zero(rop);
        return 1;
    }
    if (op->val & 1L)
    {
        return 0;
    }

    rop->val = op->val / 2;

    /*
        FIXME:  In this case, we don't actually 
        check whether the element is a square!
     */
    if (rop->val >= N)
    {
        qadic_zero(rop);
        return 1;
    }

    if (rop != op)
    {
        padic_poly_fit_length(rop, 2 * d - 1);
        t = rop->coeffs;
    }
    else
    {
        t = _fmpz_vec_init(2 * d - 1);
    }

    ans = _qadic_sqrt(t, op->coeffs, op->length, ctx->a, ctx->j, ctx->len, p, N - rop->val);

    if (rop == op)
    {
        _fmpz_vec_clear(rop->coeffs, rop->alloc);
        rop->coeffs = t;
        rop->alloc  = 2 * d - 1;
        rop->length = d;
    }
    _padic_poly_set_length(rop, d);
    _padic_poly_normalise(rop);

    return ans;
}

