#include "padic.h"

/*
    Assumes that x is non-zero, that p is greater than 1 
    and that pinv is the precomputed double inverse of p 
    whenever p is a small integer.
 */

static 
long _padic_remove(fmpz_t x, const fmpz_t p, double pinv)
{
    long e;
    fmpz y = *x;
    fmpz q = *p;

    if (!COEFF_IS_MPZ(y))  /* x is small */
    {
        if (!COEFF_IS_MPZ(q))  /* p is small */
        {
            mp_limb_t quo, rem;
            int sgn = 1;

            if (y < 0)
            {
                sgn = -1;
                y = -y;
            }

            e = 0;
            rem = n_divrem2_precomp(&quo, y, q, pinv);
            while (rem == 0)
            {
                y = quo;
                e++;
                rem = n_divrem2_precomp(&quo, y, q, pinv);
            }

            *x = sgn * y;

            return e;
        }
        else  /* p is large */
        {
            return 0;
        }
    }
    else  /* x is big */
    {
        if (!COEFF_IS_MPZ(q))  /* p is small */
        {
            if (!mpz_divisible_ui_p(COEFF_TO_PTR(y), q))
            {
                return 0;
            }
            else
            {
                mpz_t prime;

                mpz_init_set_ui(prime, q);
                e = mpz_remove(COEFF_TO_PTR(y), COEFF_TO_PTR(y), prime);
                _fmpz_demote_val(x);
                mpz_clear(prime);

                return e;
            }
        }
        else  /* p is large */
        {
            __mpz_struct *a = COEFF_TO_PTR(y);
            __mpz_struct *b = COEFF_TO_PTR(q);

            if (!mpz_divisible_p(a, b))
            {
                return 0;
            }
            else
            {
                e = mpz_remove(a, a, b);
                _fmpz_demote_val(x);
                return e;
            }
        }
    }
}

void padic_normalise(padic_t rop, const padic_ctx_t ctx)
{
    if (!fmpz_is_zero(rop) && (rop[1] < ctx->N))
    {
        rop[1] += _padic_remove(rop, ctx->p, ctx->pinv);

        if (rop[1] < ctx->N)
            _padic_reduce_unit(rop, ctx);
        else
            padic_zero(rop, ctx);
    }
    else
        padic_zero(rop, ctx);
}

