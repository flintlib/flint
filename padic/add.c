#include "padic.h"

/*
    Computes $u_1 p^{v_1} + u_2 p^{v_2} = p^{v_2} (u_1 p^{v_1 - v_2} + u_2)$.

    Assumes that \code{op1} , \code{op2} are non-zero and $v_1 > v_2$. 
 */

void _padic_add(padic_t rop, const padic_t op1, const padic_t op2, 
                const padic_ctx_t ctx)
{
    /* TODO:  Faster powering using the context. */

    fmpz_t pow;

    fmpz_init(pow);
    fmpz_pow_ui(pow, ctx->p, op1[1] - op2[1]);

    if (rop == op1)
    {
        fmpz_t t;

        fmpz_init(t);
        fmpz_set(t, op2);
        fmpz_addmul(t, op1, pow);
        fmpz_swap(rop, t);
        fmpz_clear(t);
    }
    else
    {
        fmpz_set(rop, op2);
        fmpz_addmul(rop, op1, pow);
    }

    fmpz_clear(pow);

    rop[1] = op2[1];
}

/*
    Adding u_1 p^{v_1} + u_2 p^{v_2}.
 */

void padic_add(padic_t rop, const padic_t op1, const padic_t op2, 
               const padic_ctx_t ctx)
{
    if (padic_is_zero(op1, ctx))
    {
        padic_set(rop, op2, ctx);
        return;
    }
    if (padic_is_zero(op2, ctx))
    {
        padic_set(rop, op1, ctx);
        return;
    }

    if (op1[1] > op2[1])
    {
        _padic_add(rop, op1, op2, ctx);
        _padic_reduce_unit(rop, ctx);
    }
    else if (op1[1] < op2[1])
    {
        _padic_add(rop, op2, op1, ctx);
        _padic_reduce_unit(rop, ctx);
    }
    else
    {
        fmpz_add(rop, op1, op2);
        rop[1] = op1[1];

        padic_normalise(rop, ctx);
    }
}

