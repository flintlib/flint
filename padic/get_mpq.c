#include "padic.h"

void padic_get_mpq(mpq_t rop, const padic_t op, const padic_ctx_t ctx)
{
    fmpq_t t;

    fmpq_init(t);
    padic_get_fmpq(t, op, ctx);
    fmpq_get_mpq(rop, t);
    fmpq_clear(t);
}

