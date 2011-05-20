#include "padic.h"

void padic_teichmuller(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    fmpz_t u, x, ppow;

    /* If op is not a p-adic integer, raise an abort signal. */
    if (padic_val(op) < 0)
    {
        printf("ERROR (padic_teichmuller).  op is not a p-adic integer.\n");
        abort();
    }

    /* If op is divisible by p, return zero. */
    if (padic_val(op) > 0)
    {
        padic_zero(rop, ctx);
        return;
    }

    fmpz_init(x);
    fmpz_init(u);
    fmpz_init(ppow);
    
    /* Let ppow = p^N, set x = op modulo p^N */
    fmpz_pow_ui(ppow, ctx->p, ctx->N);
    fmpz_mod(x, padic_unit(op), ppow);
    
    /* Let u be the inverse of 1-p modulo p^N */
    fmpz_sub(u, ppow, ctx->p);
    fmpz_add_ui(u, u, 1);
    fmpz_invmod(u, u, ppow);
    
    /* Let rop = x + u * (x^p - x) modulo p^N */
    fmpz_powm(rop, x, ctx->p, ppow);
    fmpz_sub(padic_unit(rop), padic_unit(rop), x);
    fmpz_mul(padic_unit(rop), u, padic_unit(rop));
    fmpz_add(padic_unit(rop), x, padic_unit(rop));
    fmpz_mod(padic_unit(rop), padic_unit(rop), ppow);
    
    /* Repeat this until rop == x modulo p^N */
    while (!fmpz_equal(padic_unit(rop), x))
    {
        fmpz_swap(x, padic_unit(rop));
        fmpz_powm(padic_unit(rop), x, ctx->p, ppow);
        fmpz_sub(padic_unit(rop), padic_unit(rop), x);
        fmpz_mul(padic_unit(rop), u, padic_unit(rop));
        fmpz_add(padic_unit(rop), x, padic_unit(rop));
        fmpz_mod(padic_unit(rop), padic_unit(rop), ppow);
    }
    
    padic_val(rop) = 0;
    
    fmpz_clear(x);
    fmpz_clear(u);
    fmpz_clear(ppow);
}

