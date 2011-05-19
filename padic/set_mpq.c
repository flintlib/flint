#include "padic.h"

void padic_set_mpq(padic_t rop, const mpq_t op, const padic_ctx_t ctx)
{
    fmpq_t t;

    fmpq_init(t);
    fmpq_set_mpq(t, op);
    padic_set_fmpq(rop, t, ctx);
    fmpq_clear(t);
}

