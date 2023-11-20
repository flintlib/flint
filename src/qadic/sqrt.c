/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mat.h"
#include "fmpz_poly.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "qadic.h"

static int __fmpz_mod_poly_invmod(fmpz *A,
                          const fmpz *B, slong lenB,
                          const fmpz *P, slong lenP, const fmpz_t p)
{
    fmpz * t, * u;
    fmpz_mod_ctx_t mod;
    int res;

    t = _fmpz_vec_init(lenB);
    u = _fmpz_vec_init(lenP);

    fmpz_mod_ctx_init(mod, p);

    _fmpz_vec_scalar_mod_fmpz(t, B, lenB, p);
    _fmpz_vec_scalar_mod_fmpz(u, P, lenP, p);

    res = _fmpz_mod_poly_invmod(A, t, lenB, u, lenP, mod);

    fmpz_mod_ctx_clear(mod);

    _fmpz_vec_clear(t, lenB);
    _fmpz_vec_clear(u, lenP);

    return res;
}

/*
    FILE DOCUMENTATION.

    This file contains routines for computing square roots in
    finite fields and Hensel lifting.

        - _find_nonresidue()
        - _artin_schreier_preimage()
        - _fmpz_mod_poly_sqrtmod_p()
        - _fmpz_mod_poly_sqrtmod_2()
        - _qadic_sqrt_p()
        - _qadic_sqrt()
        - qadic_sqrt()
 */

/*
    Sets \code{(rop,d)} to a non-residue in $\mathbf{F_q}$ for odd $p$
    and $d \geq 2$.
 */
static void _find_nonresidue(fmpz *rop,
                             const fmpz *a, const slong *j, slong lena,
                             const fmpz_t p)
{
    const slong d = j[lena - 1];

    fmpz *w;
    fmpz_t pm1, z;
    slong i;

    w = _fmpz_vec_init(2 * d - 1);
    fmpz_init(pm1);
    fmpz_init(z);

    fmpz_sub_ui(pm1, p, 1);
    fmpz_pow_ui(z, p, d);
    fmpz_sub_ui(z, z, 1);
    fmpz_fdiv_q_2exp(z, z, 1);

    /*
        Find a non-residue g;
        uses stack-based depth first traversal starting from [0,...,0,1]
     */
    _fmpz_vec_zero(rop, d);
    fmpz_one(rop + (d - 1));

    for (i = d; i > 0; )
    {
        if (i == d)
        {
            /* Consider this element, rop^{(q-1)/2} == -1 ? */
            _qadic_pow(w, rop, d, z, a, j, lena, p);
            if (fmpz_equal(w + 0, pm1) && _fmpz_vec_is_zero(w + 1, d - 1))
                break;

            /* Backtrace, find the next element */
            for (i--; i >= 0 && fmpz_equal(rop + i, pm1); i--) ;
            if (i >= 0)
            {
                fmpz_add_ui(rop + i, rop + i, 1);
                i++;
            }
        }
        else
        {
            _fmpz_vec_zero(rop + i, d - i);
            i = d;
        }
    }

    _fmpz_vec_clear(w, 2 * d - 1);
    fmpz_clear(pm1);
    fmpz_clear(z);
}

/*
    Given an element $d$ as \code{(op,len)}, returns whether it has a
    preimage $u$ under the Artin Schreier map $Y \mapsto Y^2 + Y$, and
    if so sets \code{(rop,d)} to $u$.

    In this case, the other preimage is given by $u + 1$.  This
    completes the set of preimages as the kernel of the Artin Schreier
    map is $\mathbf{F}_2$.

    The value of \code{(rop,d)}$ is undefined when the return value
    is zero.
 */
int
_artin_schreier_preimage(fmpz *rop, const fmpz *op, slong len,
                         const fmpz *a, const slong *j, slong lena)
{
    const slong d   = j[lena - 1];
    const fmpz_t p = {WORD(2)};

    int ans;

    fmpz *e, *f;
    nmod_mat_t A;
    slong i, k, *P;

    e = _fmpz_vec_init(d);
    f = _fmpz_vec_init(2 * d - 1);

    nmod_mat_init(A, d, d, 2);
    P  = flint_malloc(d * sizeof(slong));

    /* Construct Artin Schreier matrix ------------------------------------- */

    for (i = 0; i < d; i++)
    {
        fmpz_one(e + i);

        _fmpz_poly_sqr(f, e, i+1);
        _fmpz_poly_reduce(f, 2*(i+1)-1, a, j, lena);
        fmpz_add_ui(f + i, f + i, 1);
        _fmpz_vec_scalar_mod_fmpz(f, f, d, p);

        for (k = 0; k < d; k++)
        {
            nmod_mat_entry(A, k, i) = (mp_limb_t) f[k];
        }
        fmpz_zero(e + i);
    }

    /* Perform LUP decomposition ------------------------------------------- */

    /*
        Write LU = PA and r = rank(A) so that U is in row echelon form
        with only the top r rows non-zero, and L is lower unit triangular
        truncated to r columns.

        We know that Y^2 + Y = 0 if and only if Y is in the base field,
        i.e., dim ker(A) = 1 and rank(A) = d-1.

        Consider the linear system A x = b, which we can then write as
        LU x = Pb, hence L y = Pb and U x = y.
     */

#ifdef FLINT_WANT_ASSERT
    {
        slong res = nmod_mat_lu(P, A, 0);
        FLINT_ASSERT(res == d - 1);
    }
#else
    nmod_mat_lu(P, A, 0);
#endif

    /* Solve for a preimage of (op,len) ------------------------------------ */

    _fmpz_vec_zero(rop, d);

    for (i = 0; i < d; i++)
    {
        rop[i] = P[i] < len ? op[P[i]] : 0;
        for (k = 0; k < i; k++)
            rop[i] ^= nmod_mat_entry(A, i, k) & rop[k];
    }

    ans = (rop[d-1] == 0);

    if (ans)
    {
        slong c;

        for (c = 0; c < d; c++)
            if (nmod_mat_entry(A, c, c) == 0)
                break;

        for (k = d - 1; k > c; k--)
        {
            rop[k] = rop[k-1];

            if ((rop[k]))
                for (i = k - 2; i >= 0; i--)
                    rop[i] ^= nmod_mat_entry(A, i, k);
        }
        rop[k] = 0;
        for (k--; k >= 0; k--)
        {
            if ((rop[k]))
                for (i = k - 1; i >= 0; i--)
                    rop[i] ^= nmod_mat_entry(A, i, k);
        }
    }

    /* Clean-up ------------------------------------------------------------ */

    _fmpz_vec_clear(e, d);
    _fmpz_vec_clear(f, 2 * d - 1);
    nmod_mat_clear(A);
    flint_free(P);

    return ans;
}

/*
    Returns whether the non-zero element \code{(op, len)} is a square,
    and if so sets \code{(rop, 2 * d - 1)} to its square root.

    Assumes that $p$ is an odd prime.

    Assumes that $d \geq 1$.

    Does not support aliasing.
 */
static int
_fmpz_mod_poly_sqrtmod_p(fmpz *rop, const fmpz *op, slong len,
                         const fmpz *a, const slong *j, slong lena,
                         const fmpz_t p)
{
    const slong d = j[lena - 1];
    int ans;

    /*
        When $q \equiv 3 \pmod{4}$...

        A non-zero element $x$ is a square if and only if $x^{(q-1)/2} = 1$,
        and in this case one of its square roots is given by $x^{(q+1)/4}$.

        To avoid recomputation of powers of $x$, we compute $x^{(q-3)/4}$,
        multiply this by $x$ to obtain the potential square root $x^{(q+1)/4}$,
        and then combine these two powers to find $x^{(q-1)/2}$.
     */
    if (fmpz_fdiv_ui(p, 4) == 3 && (d & WORD(1)))
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

        _fmpz_poly_mul(rop, v, d, op, len);
        _fmpz_vec_zero(rop + d + len - 1, d - len);
        _fmpz_mod_poly_reduce(rop, d + len - 1, a, j, lena, p);

        _fmpz_poly_mul(w, rop, d, v, d);
        _fmpz_mod_poly_reduce(w, 2 * d - 1, a, j, lena, p);
        ans = fmpz_is_one(w + 0);

        fmpz_clear(z);
        _fmpz_vec_clear(v, 4 * d - 2);
    }

    /*
        When $q \equiv 5 \pmod{8}$...

        A non-zero element $x$ is a square if and only if $y = x^{(q-1)/4}$
        is $\pm 1$.  If $y = +1$, a square root is given by $x^{(q+3)/8}$,
        otherwise it is $(2x) (4x)^{(q-5)/8}$.

        We begin by computing $v = x^{(q-5)/8}$, $w = x v = x^{(q+3)/8}$
        and $y = v w = x^{(q-1)/4}$.  If $y = +1$ we return $w$, and if
        $y = -1$ we return $2^{(q-1)/4} x^{(q+3)/8} = 2^{(q-1)/4} w$.
     */
    else if (fmpz_fdiv_ui(p, 8) == 5 && (d & WORD(1)))
    {
        fmpz_t f, g, pm1;

        fmpz *v;  /* v = x^{(q-5)/8} */
        fmpz *w;  /* w = x^{(q+3)/8} */
        fmpz *y;  /* y = x^{(q-1)/4} */

        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(pm1);
        fmpz_sub_ui(pm1, p, 1);

        v = _fmpz_vec_init(2 * d - 1);
        w = _fmpz_vec_init(2 * d - 1);
        y = _fmpz_vec_init(2 * d - 1);

        fmpz_pow_ui(g, p, d);
        fmpz_sub_ui(g, g, 5);
        fmpz_fdiv_q_2exp(g, g, 3);
        _qadic_pow(v, op, len, g, a, j, lena, p);

        _fmpz_poly_mul(w, v, d, op, len);
        _fmpz_mod_poly_reduce(w, d + len - 1, a, j, lena, p);

        _fmpz_poly_mul(y, v, d, w, d);
        _fmpz_mod_poly_reduce(y, 2 * d - 1, a, j, lena, p);

        if ((fmpz_is_one(y + 0) || fmpz_equal(y + 0, pm1))  /* q.r. */
            && _fmpz_vec_is_zero(y + 1, d - 1))
        {
            if (fmpz_is_one(y + 0))
            {
                _fmpz_vec_set(rop, w, d);
            }
            else
            {
                fmpz_t two = {WORD(2)};

                fmpz_zero(g);
                fmpz_pow_ui(g, p, d);
                fmpz_sub_ui(g, g, 1);
                fmpz_fdiv_q_2exp(g, g, 2);
                _fmpz_vec_set(rop, v, d);

                fmpz_powm(f, two, g, p);

                _fmpz_vec_scalar_mul_fmpz(rop, w, d, f);
                _fmpz_vec_scalar_mod_fmpz(rop, rop, d, p);
            }
            _fmpz_vec_zero(rop + d, d - 1);
            ans = 1;
        }
        else  /* q.n.r. */
        {
            ans = 0;
        }

        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(pm1);

        _fmpz_vec_clear(v, 2 * d - 1);
        _fmpz_vec_clear(w, 2 * d - 1);
        _fmpz_vec_clear(y, 2 * d - 1);
    }

    /*
        General case for odd $q$...

        TODO:  Find a better way to integrate the check for square-ness
        into the computation of a potential square root.
     */
    else
    {
        slong i, s;

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

        if (ans)
        {
            /* Find a non-residue g */
            _find_nonresidue(g, a, j, lena, p);

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
                slong k;

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
        }

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

    Note that the group of units of $\mathbf{F}_q$ is cyclic of order $q - 1$,
    which is odd if $p = 2$.  We have $x^{q-1} = 1$ for every non-zero $x$, so
    $x^q = x$ for every $x$ and $u = x^{q/2}$ satisfies $u^2 = x$.

    Assumes that $d \geq 1$.
 */
static void
_fmpz_mod_poly_sqrtmod_2(fmpz *rop, const fmpz *op, slong len,
                         const fmpz *a, const slong *j, slong lena)
{
    const fmpz_t p = {WORD(2)};
    const slong d   = j[lena - 1];

    fmpz_t z;

    fmpz_init(z);
    fmpz_setbit(z, d - 1);
    _qadic_pow(rop, op, len, z, a, j, lena, p);
    fmpz_clear(z);
}

/*
    Returns whether \code{(op, len)} is a square, and if so
    sets \code{(rop, 2 * d - 1)} to a square root mod $p^N$.

    Assumes that \code{(op, len)} has valuation $0$.
 */
static int
_qadic_sqrt_p(fmpz *rop, const fmpz *op, slong len,
              const fmpz *a, const slong *j, slong lena,
              const fmpz_t p, slong N)
{
    const slong d = j[lena - 1];
    int ans;

    if (N == 1)
    {
        ans = _fmpz_mod_poly_sqrtmod_p(rop, op, len, a, j, lena, p);
        return ans;
    }
    else
    {
        slong *e, i, k, n;
        fmpz *pow, *u;
        fmpz *r, *s, *t, *f;

        n = FLINT_CLOG2(N) + 1;

        /* Compute sequence of exponents */
        e = flint_malloc(n * sizeof(slong));
        for (e[i = 0] = N; e[i] > 1; i++)
            e[i + 1] = (e[i] + 1) / 2;

        pow = _fmpz_vec_init(n);
        u   = _fmpz_vec_init(len * n);
        r   = _fmpz_vec_init(2 * d - 1);
        s   = _fmpz_vec_init(2 * d - 1);
        t   = _fmpz_vec_init(2 * d - 1);
        f   = _fmpz_vec_init(d + 1);

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
                fmpz_set(f + j[k], a + k);
            __fmpz_mod_poly_invmod(rop, t, d, f, d + 1, p);
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
        _fmpz_vec_clear(f, d + 1);
        flint_free(e);

        return ans;
    }
}

/*
    Returns whether \code{(op, len)} is a square, and if so
    sets \code{(rop, 2 * d - 1)} to a square root mod $2^N$.

    Assumes that \code{(op, len)} has valuation $0$.
 */
static int
_qadic_sqrt_2(fmpz *rop, const fmpz *op, slong len,
              const fmpz *a, const slong *j, slong lena, slong N)
{
    int ans;

    const slong d    = j[lena - 1];
    const fmpz_t p  = {WORD(2)};

    fmpz *g, *r, *s, *t;
    slong i;

    g = _fmpz_vec_init(2 * d - 1);
    r = _fmpz_vec_init(2 * d - 1);
    s = _fmpz_vec_init(2 * d - 1);
    t = _fmpz_vec_init(2 * d - 1);

    /* Compute r = 1/op (mod 8) */
    _qadic_inv(r, op, len, a, j, lena, p, 3);

    /* Compute s = invsqrt(op) (mod 2) */
    _fmpz_vec_scalar_fdiv_r_2exp(t, r, d, 1);
    _fmpz_mod_poly_sqrtmod_2(s, r, d, a, j, lena);

    /* Compute t = (1/op - s^2) / 4 (mod 2) if that makes sense */
    _fmpz_poly_sqr(t, s, d);
    _fmpz_poly_reduce(t, 2*d - 1, a, j, lena);
    _fmpz_vec_sub(t, r, t, d);

    _fmpz_vec_zero(rop, 2 * d - 1);
    ans = 1;
    for (i = 0; i < d; i++)
    {
        if (fmpz_val2(t + i) == 1)
            ans = 0;
    }

    if (ans) {
        _fmpz_vec_scalar_fdiv_q_2exp(t, t, d, 2);
        _fmpz_vec_scalar_fdiv_r_2exp(t, t, d, 1);

        /* Check whether X^2 + s X = t (mod 2) has a solution */
        /*
            u = (op,len) is a square in Zq iff it is a square modulo 8,
            which is the case iff X^2 + s X + t == 0 (mod 2) is soluble.
            Let g be t/s^2.  Then the above quadratic is soluble iff
            Y^2 + Y + d = 0 is soluble.

            We can observe that 1/s^2 = op (mod 2).
         */
        _fmpz_vec_scalar_fdiv_r_2exp(r, op, len, 1);
        _fmpz_poly_mul(g, t, d, r, len);
        _fmpz_mod_poly_reduce(g, 2 * d - 1, a, j, lena, p);

        ans = _artin_schreier_preimage(r, g, d, a, j, lena);

        if (ans)
        {
            if (N == 1)
            {
                _fmpz_mod_poly_sqrtmod_2(rop, op, len, a, j, lena);
            }
            else
            {

                /* Set t = r s, a root of X^2 + s X = t (mod 2) */
                _fmpz_poly_mul(t, r, d, s, d);
                _fmpz_mod_poly_reduce(t, 2 * d - 1, a, j, lena, p);

                /* Now s + 2t is an invsqrt of (op, len) to precision 4. */
                _fmpz_vec_scalar_addmul_si(s, t, d, 2);

                if (N == 2)
                {
                    _qadic_inv(rop, s, d, a, j, lena, p, 2);
                }
                else  /* N >= 3 */
                {
                    slong *e, n;
                    fmpz *u;

                    n = FLINT_CLOG2(N - 1) + 1;

                    /* Compute sequence of exponents */
                    e = flint_malloc(n * sizeof(slong));
                    for (e[i = 0] = N; e[i] > 2; i++)
                        e[i + 1] = e[i] / 2 + 1;

                    u = _fmpz_vec_init(len * n);

                    /* Compute reduced units */
                    {
                        _fmpz_vec_scalar_fdiv_r_2exp(u + 0 * len, op, len, e[0]);
                    }
                    for (i = 1; i < n; i++)
                    {
                        _fmpz_vec_scalar_fdiv_r_2exp(u + i * len, u + (i - 1) * len, len, e[i] + 1);
                    }

                    /* Run Newton iteration */
                    {
                        _fmpz_vec_set(rop, s, d);
                    }
                    for (i = n - 2; i >= 1; i--)  /* z := z - z (a z^2 - 1) / 2 */
                    {
                        _fmpz_poly_sqr(s, rop, d);
                        _fmpz_poly_reduce(s, 2 * d - 1, a, j, lena);
                        _fmpz_poly_mul(t, s, d, u + i * len, len);
                        _fmpz_poly_reduce(t, d + len - 1, a, j, lena);
                        fmpz_sub_ui(t, t, 1);
                        _fmpz_vec_scalar_fdiv_q_2exp(t, t, d, 1);
                        _fmpz_poly_mul(s, t, d, rop, d);
                        _fmpz_poly_reduce(s, 2 * d - 1, a, j, lena);
                        _fmpz_poly_sub(rop, rop, d, s, d);
                        _fmpz_vec_scalar_fdiv_r_2exp(rop, rop, d, e[i]);
                    }
                    /* z := a z + z (a - (a z)^2) / 2 */
                    {
                        _fmpz_poly_mul(s, rop, d, u + len, len);
                        _fmpz_poly_reduce(s, d + len - 1, a, j, lena);
                        _fmpz_poly_sqr(t, s, d);
                        _fmpz_poly_reduce(t, 2 * d - 1, a, j, lena);
                        _fmpz_poly_sub(t, u + 0 * len, len, t, d);
                        _fmpz_vec_scalar_fdiv_q_2exp(t, t, d, 1);
                        _fmpz_poly_mul(r, rop, d, t, d);
                        _fmpz_poly_reduce(r, 2 * d - 1, a, j, lena);
                        _fmpz_poly_add(rop, r, d, s, d);
                        _fmpz_vec_scalar_fdiv_r_2exp(rop, rop, d, e[0]);
                    }

                    _fmpz_vec_clear(u, len * n);
                    flint_free(e);
                }
            }
        }
    }

    _fmpz_vec_clear(g, 2 * d - 1);
    _fmpz_vec_clear(r, 2 * d - 1);
    _fmpz_vec_clear(s, 2 * d - 1);
    _fmpz_vec_clear(t, 2 * d - 1);

    return ans;
}

/*
    Returns whether \code{(op, len)} is a square, and if so
    sets \code{(rop, 2 * d - 1)} to a square root, reduced
    modulo $2^N$.

    Assumes that \code{(op, len)} has valuation $0$.
 */
int _qadic_sqrt(fmpz *rop, const fmpz *op, slong len,
                const fmpz *a, const slong *j, slong lena,
                const fmpz_t p, slong N)
{
    if (*p == WORD(2))
    {
        return _qadic_sqrt_2(rop, op, len, a, j, lena, N);
    }
    else
    {
        return _qadic_sqrt_p(rop, op, len, a, j, lena, p, N);
    }
}

int qadic_sqrt(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
{
    const fmpz *p   = (&ctx->pctx)->p;
    const slong d   = qadic_ctx_degree(ctx);
    const slong N   = qadic_prec(rop);
    const slong vop = qadic_val(op);

    fmpz *t;
    int ans;

    if (qadic_is_zero(op))
    {
        qadic_zero(rop);
        return 1;
    }
    if (op->val & WORD(1))
    {
        return 0;
    }

    rop->val = op->val / 2;

    if (rop == op)
    {
        t = _fmpz_vec_init(2*d - 1);
    }
    else
    {
        padic_poly_fit_length(rop, 2*d - 1);
        t = rop->coeffs;
    }

    if (rop->val >= N) /* must check square even if rop prec will be <= valuation */
    {
        slong prec = FLINT_MIN(qadic_prec(op) - vop, 1 + (fmpz_cmp_ui(p, 2) == 0));
        ans = prec <= 0 ? 1 : _qadic_sqrt(t, op->coeffs, op->length, ctx->a, ctx->j, ctx->len, p, prec);
    }
    else
        ans = _qadic_sqrt(t, op->coeffs, op->length, ctx->a, ctx->j, ctx->len, p, N - rop->val);

    if (rop == op)
    {
        _fmpz_vec_clear(rop->coeffs, rop->alloc);

        rop->coeffs = t;
        rop->alloc  = 2 * d - 1;
        rop->length = d;
    }

    if (rop->val >= N)
        qadic_zero(rop);
    else
    {
        _padic_poly_set_length(rop, d);
        _padic_poly_normalise(rop);
        if (padic_poly_length(rop) == 0)
            padic_poly_val(rop) = 0;
    }

    return ans;
}
