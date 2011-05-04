#include "padic.h"

void padic_sub(padic_t rop, const padic_t op1, const padic_t op2, 
               const padic_ctx_t ctx)
{
    if (padic_is_zero(op2, ctx))
    {
        padic_set(rop, op1, ctx);
        return;
    }

    if (padic_is_zero(op1, ctx))
    {
        padic_neg(rop, op2, ctx);
        return;
    }

    /* Aliasing */
    if (rop == op1 || rop == op2)
    {
        padic_t t;

        padic_init(t, ctx);
        padic_sub(t, op1, op2, ctx);
        padic_swap(rop, t, ctx);
        padic_clear(t, ctx);
        return;
    }

    if (op1[1] > op2[1])
    {
        /* u1 p^v1 - u2 p^v2 = p^v2 (u1 p^{v1-v2} - u2) */
        fmpz_t pow;

        fmpz_init(pow);
        fmpz_pow_ui(pow, ctx->p, op1[1] - op2[1]);
        fmpz_neg(rop, op2);
        fmpz_addmul(rop, op1, pow);
        fmpz_clear(pow);

        rop[1] = op2[1];

        _padic_reduce_unit(rop, ctx);
    }
    else if (op1[1] < op2[1])
    {
        /* u1 p^v1 - u2 p^v2 = p^v1 (u1 - u2 p^{v2-v1}) */
        fmpz_t pow;

        fmpz_init(pow);
        fmpz_pow_ui(pow, ctx->p, op2[1] - op1[1]);
        fmpz_set(rop, op1);
        fmpz_submul(rop, op2, pow);
        fmpz_clear(pow);

        rop[1] = op1[1];

        _padic_reduce_unit(rop, ctx);
    }
    else
    {
        fmpz_sub(rop, op1, op2);
        rop[1] = op1[1];

        padic_normalise(rop, ctx);
    }
}

