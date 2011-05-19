#include "padic.h"

void padic_val_fac(fmpz_t rop, const fmpz_t op, const fmpz_t p)
{
    fmpz_t t, q, pow;

    if (fmpz_sgn(op) <= 0)
    {
        printf("Exception (padic_val_fac).  op is non-positive.\n");
        abort();
    }

    fmpz_init(t);
    fmpz_init(q);
    fmpz_init(pow);
    fmpz_set_ui(pow, 1);

    do 
    {
        fmpz_mul(pow, pow, p);
        fmpz_fdiv_q(q, op, pow);
        fmpz_add(t, t, q);
    }
    while (!fmpz_is_zero(q));

    fmpz_swap(rop, t);
    fmpz_clear(t);
    fmpz_clear(q);
    fmpz_clear(pow);
}

