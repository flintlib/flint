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

    if (padic_val(op1) > padic_val(op2))
    {
        /* u1 p^v1 - u2 p^v2 = p^v2 (u1 p^{v1-v2} - u2) */
        fmpz_t pow;
        int alloc;

        _padic_ctx_pow_ui(pow, &alloc, padic_val(op1) - padic_val(op2), ctx);
        fmpz_neg(padic_unit(rop), padic_unit(op2));
        fmpz_addmul(padic_unit(rop), padic_unit(op1), pow);
        if (alloc)
            fmpz_clear(pow);

        padic_val(rop) = padic_val(op2);

        _padic_reduce_unit(rop, ctx);
    }
    else if (padic_val(op1) < padic_val(op2))
    {
        /* u1 p^v1 - u2 p^v2 = p^v1 (u1 - u2 p^{v2-v1}) */
        fmpz_t pow;
        int alloc;

        _padic_ctx_pow_ui(pow, &alloc, padic_val(op2) - padic_val(op1), ctx);
        fmpz_set(padic_unit(rop), padic_unit(op1));
        fmpz_submul(padic_unit(rop), padic_unit(op2), pow);
        if (alloc)
            fmpz_clear(pow);

        padic_val(rop) = padic_val(op1);

        _padic_reduce_unit(rop, ctx);
    }
    else
    {
        fmpz_sub(padic_unit(rop), padic_unit(op1), padic_unit(op2));
        padic_val(rop) = padic_val(op1);

        padic_normalise(rop, ctx);
    }
}

