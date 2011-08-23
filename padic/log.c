#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "padic.h"

/*
    Returns $b$ such that for all $i \geq b$ we have 
    \begin{equation*}
    i v - \ord_p(i) \geq N
    \end{equation*}
    where $v \geq 1$.

    Write $p^{k-1} < N/v <= p^{k}$ and then let $i = p^{k} + x$ 
    with $x < (p - 1) p^k$ so that $i < p^{k+1}$.  Then 
    \begin{equation*}
    i v - \ord_p(i) - N \geq p^k v + x v - k - N \geq x v - k
    \end{equation*}
    so we need $x >= k / v$.

    Thus, have the condition that $\ceil{k/v} < (p-1)p^k$, 
    which is always satisfied.

    Assumes that $v \geq 1$.

    Assumes that the result fits into a \code{long}.  Note that 
    the result is usually a little greater than $\ceil{N/v}$.

    Does not guarantee that the return value is the minimial $b$ 
    with the above property.

    N.B.  If $x$ is a $p$-adic integer s.t. $x - 1$ is divisible 
    by $p$ so that $\log(x)$ converges, then in order to compute 
    $\log(x)$ modulo $p^N$ it suffices to evaluate 
    \begin{equation*}
    \sum_{i=1}^{b-1} (-1)^{i-1} \frac{(x-1)^i}{i}.
    \end{equation*}
 */

static long bound(long v, long N, long p)
{
/*    printf("v N p = %ld %ld %ld\n", v, N, p); fflush(stdout); */
    if (v >= N)
    {
        return 1;
    }
    else
    {
        long k = n_clog((N + (v - 1)) / v, p);

        return n_pow(p, k) + (k + (v - 1)) / v;
    }
}

/*
    Computes 
    \begin{align*}
    \log(x) & = \sum_{i=1}^{\infty} (-1)^{i-1} \frac{(x-1)^i}{i} \\
            & = - \sum_{i=1}^{\infty} \frac{(1-x)^i}{i}
    \end{align*}
 */
void padic_log(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (_padic_is_one(op))
    {
        _padic_zero(rop);
        return;
    }

    if (padic_val(op) < 0)
    {
        printf("ERROR (padic_log).  op has negative valuation.\n");
        abort();
    }

    if (fmpz_cmp_ui(ctx->p, 2) == 0)
    {
        printf("ERROR (padic_log).  p = 2 case not implemented yet.\n");
        abort();
    }
    else
    {
        fmpz_t y;
        long v = fmpz_remove(padic_unit(op), y, ctx->p);

        fmpz_init(y);

        /* y = 1 - op */
        padic_get_fmpz(y, op, ctx);
        fmpz_sub_ui(y, y, 1);
        fmpz_neg(y, y);

        if (v <= 0)
        {
            printf("ERROR (padic_log).  Series does not converge.\n");
            abort();
        }

        if (ctx->N <= v)
        {
            _padic_zero(rop);
        }
        else if (fmpz_fits_si(ctx->p))
        {
            const long p = fmpz_get_si(ctx->p);
            long i = bound(padic_val(y), ctx->N, p) - 1;
            long k = n_flog(i, p);
            long j;
            long e;

            fmpz *q;
            fmpz_t m, s, t;

            q = _fmpz_vec_init(k + 1);
            fmpz_set_ui(q + 0, 1);
            for (j = 1; j <= k; j++)
                fmpz_mul_ui(q + j, q + (j - 1), p);
            fmpz_init(m);
            fmpz_pow_ui(m, ctx->p, ctx->N + k);
            fmpz_init(s);
            fmpz_init(t);

            /* rop := i^{-1} */
            j = i;
            e = n_remove((mp_limb_t *) &j, p);
            _padic_inv(s, (fmpz *) &j, ctx->p, ctx->N + k);
            fmpz_mul(padic_unit(rop), s, q + (k - e));

            for (i--; i > 0; i--)
            {
                /* t := rop * y */
                fmpz_mul(t, padic_unit(rop), y);

                /* rop := i^{-1} + t */
                j = i;
                e = n_remove((mp_limb_t *) &j, p);
                _padic_inv(s, (fmpz *) &j, ctx->p, ctx->N + k);
                fmpz_mul(padic_unit(rop), s, q + (k - e));

                fmpz_add(padic_unit(rop), padic_unit(rop), t);
                fmpz_mod(padic_unit(rop), padic_unit(rop), m);
            }

            fmpz_mul(padic_unit(rop), padic_unit(rop), y);
            fmpz_divexact(padic_unit(rop), padic_unit(rop), q + k);
            padic_val(rop) = 0;
            fmpz_neg(rop, rop);
            padic_reduce(rop, ctx);

            _fmpz_vec_clear(q, k + 1);
            fmpz_clear(m);
            fmpz_clear(s);
            fmpz_clear(t);
        }
        else
        {
            printf("ERROR (padic_log).  Only implemented for small p.\n");
            abort();
        }

        fmpz_clear(y);
    }
}

