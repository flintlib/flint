#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "padic.h"
#include "ulong_extras.h"

extern long _padic_log_bound(long v, long N, long p);

/*
    Carries out the finite series evaluation for the logarithm 
    \begin{equation*}
    \sum_{i=1}^{n} a_i x^i
    = \sum_{j=0}^{\ceil{n/b} - 1} \Bigl( \sum_{i=1}^b a_{i+jb} x^i \Bigr) x^{jb}
    \end{equation*}
    where $a_i = 1/i$ with the choice $b = \floor{\sqrt{n}}$, 
    all modulo $p^N$, where also $P = p^N$.

    Does not support aliasing.
 */
static void 
_padic_log_rectangular_series(fmpz_t z, const fmpz_t y, long n, 
                              const fmpz_t p, long N, const fmpz_t P)
{
    if (n <= 2)
    {
        if (n == 1)
        {
            fmpz_mod(z, y, P);
        }
        else  /* n == 2;  z = y(1 + y/2) */
        {
            if (fmpz_is_even(y))
            {
                fmpz_fdiv_q_2exp(z, y, 1);
            }
            else  /* => p and y are odd */
            {
                fmpz_add(z, y, P);
                fmpz_fdiv_q_2exp(z, z, 1);
            }
            fmpz_add_ui(z, z, 1);
            fmpz_mul(z, z, y);
            fmpz_mod(z, z, P);
        }
    }
    else
    {
        const long b = n_sqrt(n);

        long i, j;
        fmpz_t c, f, m, t;
        fmpz *ypow;

        ypow = _fmpz_vec_init(b + 1);
        fmpz_init(c);
        fmpz_init(f);
        fmpz_init(m);
        fmpz_init(t);

        fmpz_one(ypow + 0);
        for (i = 1; i <= b; i++)
            fmpz_mul(ypow + i, ypow + (i - 1), y);

        fmpz_zero(z);
        fmpz_one(m);

        for (j = 0; j < (n + (b-1)) / b; j++)
        {
            const long hi = FLINT_MIN(b, n - j*b);

            fmpz_rfac_uiui(f, 1 + j*b, hi);

            fmpz_zero(c);
            for (i = 1; i <= hi; i++)
            {
                fmpz_divexact_ui(t, f, i+j*b);
                fmpz_addmul(c, t, ypow + i);
            }

            i = fmpz_remove(f, f, p);
            _padic_inv(f, f, p, N);
            fmpz_pow_ui(t, p, i);
            fmpz_divexact(c, c, t);
            fmpz_mul(c, c, f);

            fmpz_mul(c, c, m);
            fmpz_add(z, z, c);
            fmpz_mod(z, z, P);

            fmpz_mul(m, m, ypow + b);
        }

        fmpz_clear(c);
        fmpz_clear(f);
        fmpz_clear(m);
        fmpz_clear(t);
        _fmpz_vec_clear(ypow, b + 1);
    }
}

/*
    Computes 
    \begin{equation*}
    z = \sum_{i = 1}^{\infty} \frac{y^i}{i} \pmod{p^N}.
    \end{equation*}

    Note that this can be used to compute the $p$-adic logarithm 
    via the equation 
    \begin{align*}
    \log(x) & = \sum_{i=1}^{\infty} (-1)^{i-1} \frac{(x-1)^i}{i} \\
            & = - \sum_{i=1}^{\infty} \frac{(1-x)^i}{i}.
    \end{align*}

    Assumes that $y = 1 - x$ is non-zero and that $v = \ord_p(y)$ 
    is at least $1$ when $p$ is odd and at least $2$ when $p = 2$ 
    so that the series converges.

    Assumes that $v < N$.

    Does not support aliasing between $y$ and $z$.
 */
static void _padic_log_rectangular(fmpz_t z, const fmpz_t y, long v, const padic_ctx_t ctx)
{
    if (fmpz_fits_si(ctx->p))
    {
        const long p = fmpz_get_si(ctx->p);
        const long i = _padic_log_bound(v, ctx->N, p) - 1;
        const long k = n_flog(i, p);

        fmpz_t P;

        fmpz_init(P);
        fmpz_pow_ui(P, ctx->p, ctx->N + k);
        _padic_log_rectangular_series(z, y, i, ctx->p, ctx->N + k, P);
        fmpz_pow_ui(P, ctx->p, ctx->N);
        fmpz_mod(z, z, P);
        fmpz_clear(P);
    }
    else
    {
        /*
            When p does not fit into a signed long, 
            p does not divide the index i.

            Assumes that (N - 1) / v is a small 
            fmpz integer.

            TODO:  Fix this part.
         */
        long i;
        fmpz_t m, t;

        i = (ctx->N - 1) / v;

        fmpz_init(m);
        fmpz_init(t);

        fmpz_pow_ui(m, ctx->p, ctx->N);

        fmpz_zero(z);

        for ( ; i > 0; i--)
        {
            fmpz_mul(t, z, y);

            _padic_inv(z, (fmpz *) &i, ctx->p, ctx->N);

            fmpz_add(z, z, t);
            fmpz_mod(z, z, m);
        }

        fmpz_mul(z, z, y);

        fmpz_clear(m);
        fmpz_clear(t);
    }
}

int padic_log_rectangular(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (padic_val(op) < 0)
    {
        return 0;
    }
    else
    {
        fmpz_t x;
        int ans;

        fmpz_init(x);

        padic_get_fmpz(x, op, ctx);
        fmpz_sub_ui(x, x, 1);
        fmpz_neg(x, x);

        if (fmpz_is_zero(x))
        {
            padic_zero(rop);
            ans = 1;
        }
        else
        {
            fmpz_t t;
            long v;

            fmpz_init(t);
            v = fmpz_remove(t, x, ctx->p);
            fmpz_clear(t);

            if ((*(ctx->p) == 2L && v >= 2) || v >= 1)
            {
                if (v >= ctx->N)
                {
                    padic_zero(rop);
                }
                else
                {
                    _padic_log_rectangular(padic_unit(rop), x, v, ctx);
                    fmpz_neg(padic_unit(rop), padic_unit(rop));
                    padic_val(rop) = 0;
                    padic_reduce(rop, ctx);
                }
                ans = 1;
            }
            else
            {
                ans = 0;
            }
        }

        fmpz_clear(x);
        return ans;
    }
}

