#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "padic.h"
#include "ulong_extras.h"

/*
    Returns $b$ such that for all $i \geq b$ we have 
    \begin{equation*}
    i v - \ord_p(i) \geq N
    \end{equation*}
    where $v \geq 1$.

    Assumes that $N + (v - 1)$, $b v$, and $N + e$ 
    do not overflow, where $e = \floor{\log_{p}{b}}$.
 */
static long bound(long v, long N, const fmpz_t p)
{
    if (v >= N)
    {
        return 1;
    }
    else
    {
        long i = (N + (v - 1)) / v;

        if (fmpz_fits_si(p))
        {
            long e, j, q = fmpz_get_si(p);

            --i;
            do 
            {
                j = ++i;
                e = n_remove(&j, q);
            }
            while (i * v < N + e);
        }
        return i;
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
        long v;

        fmpz_init(y);

        /* y = 1 - op */
        padic_get_fmpz(y, op, ctx);
        fmpz_sub_ui(y, y, 1);
        fmpz_neg(y, y);

        if (fmpz_is_zero(y))
        {
            padic_zero(rop, ctx);
            fmpz_clear(y);
            return;
        }

        v = fmpz_remove(padic_unit(rop), y, ctx->p);

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
            long e, i, j, k, p;
            fmpz_t m, s, t;
            fmpz *q;

            p = fmpz_get_si(ctx->p);
            i = bound(v, ctx->N, ctx->p) - 1;
            k = n_flog(i, p);

            fmpz_init(m);
            fmpz_init(s);
            fmpz_init(t);
            q = _fmpz_vec_init(k + 1);

            fmpz_pow_ui(m, ctx->p, ctx->N + k);
            fmpz_set_ui(q + 0, 1);
            for (j = 1; j <= k; j++)
                fmpz_mul_ui(q + j, q + (j - 1), p);

            fmpz_zero(padic_unit(rop));

            for ( ; i > 0; i--)
            {
                fmpz_mul(t, padic_unit(rop), y);

                j = i;
                e = n_remove((mp_limb_t *) &j, p);
                _padic_inv(s, (fmpz *) &j, ctx->p, ctx->N + k);
                fmpz_mul(padic_unit(rop), s, q + (k - e));

                fmpz_add(padic_unit(rop), padic_unit(rop), t);
                fmpz_mod(padic_unit(rop), padic_unit(rop), m);
            }

            fmpz_divexact(padic_unit(rop), padic_unit(rop), q + k);
            fmpz_mul(padic_unit(rop), padic_unit(rop), y);
            fmpz_neg(padic_unit(rop), padic_unit(rop));
            padic_val(rop) = 0;
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

