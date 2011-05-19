#include "padic.h"

/*
    Computes $u_1 p^{v_1} + u_2 p^{v_2} = p^{v_2} (u_1 p^{v_1 - v_2} + u_2)$.

    Assumes that \code{op1} , \code{op2} are non-zero and $v_1 > v_2$. 

    Does not reduce the unit part.

    Supports aliasing.
 */

void _padic_add(padic_t rop, const padic_t op1, const padic_t op2, 
                const padic_ctx_t ctx)
{
    fmpz_t pow;
    int alloc;

    _padic_ctx_pow_ui(pow, &alloc, padic_val(op1) - padic_val(op2), ctx);

    if (rop == op1)
    {
        fmpz_t t;

        fmpz_init(t);
        fmpz_set(t, padic_unit(op2));
        fmpz_addmul(t, padic_unit(op1), pow);
        fmpz_swap(padic_unit(rop), t);
        fmpz_clear(t);
    }
    else
    {
        fmpz_set(padic_unit(rop), padic_unit(op2));
        fmpz_addmul(padic_unit(rop), padic_unit(op1), pow);
    }

    padic_val(rop) = padic_val(op2);

    if (alloc)
        fmpz_clear(pow);
}

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

    if (padic_val(op1) > padic_val(op2))
    {
        _padic_add(rop, op1, op2, ctx);
        _padic_reduce_unit(rop, ctx);
    }
    else if (padic_val(op1) < padic_val(op2))
    {
        _padic_add(rop, op2, op1, ctx);
        _padic_reduce_unit(rop, ctx);
    }
    else
    {
        fmpz_add(padic_unit(rop), padic_unit(op1), padic_unit(op2));
        padic_val(rop) = padic_val(op1);

        padic_normalise(rop, ctx);
    }
}

