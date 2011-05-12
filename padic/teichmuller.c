#include "padic.h"

void padic_teichmuller(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    fmpz_t u, x, ppow;

    /* If op is not a p-adic integer, raise an abort signal. */
    if (op[1] < 0)
    {
        printf("ERROR (padic_teichmuller).  op is not a p-adic integer.\n");
        abort();
    }

    /* If op is divisible by p, return zero. */
    if (op[1] > 0)
    {
        padic_zero(rop, ctx);
        return;
    }

    fmpz_init(x);
    fmpz_init(u);
    fmpz_init(ppow);
    
    /* Let ppow = p^N, set x = op modulo p^N */
    fmpz_pow_ui(ppow, ctx->p, ctx->N);
    fmpz_mod(x, op, ppow);
    
    /* Let u be the inverse of 1-p modulo p^N */
    fmpz_sub(u, ppow, ctx->p);
    fmpz_add_ui(u, u, 1);
    fmpz_invmod(u, u, ppow);
    
    /* Let rop = x + u * (x^p - x) modulo p^N */
    fmpz_powm(rop, x, ctx->p, ppow);
    fmpz_sub(rop, rop, x);
    fmpz_mul(rop, u, rop);
    fmpz_add(rop, x, rop);
    fmpz_mod(rop, rop, ppow);
    
    /* Repeat this until rop == x modulo p^N */
    while (!fmpz_equal(rop, x))
    {
        fmpz_swap(x, rop);
        fmpz_powm(rop, x, ctx->p, ppow);
        fmpz_sub(rop, rop, x);
        fmpz_mul(rop, u, rop);
        fmpz_add(rop, x, rop);
        fmpz_mod(rop, rop, ppow);
    }
    
    rop[1] = 0;
    
    fmpz_clear(x);
    fmpz_clear(u);
    fmpz_clear(ppow);
}

