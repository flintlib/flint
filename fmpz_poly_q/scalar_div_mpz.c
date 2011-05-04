#include "fmpz_poly_q.h"

void fmpz_poly_q_scalar_div_mpz(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const mpz_t x)
{
    fmpz_t y;

    if (mpz_sgn(x) == 0)
    {
        printf("ERROR (fmpz_poly_q_scalar_div_mpz).  Division by zero.\n");
        abort();
    }

    fmpz_init(y);
    fmpz_set_mpz(y, x);

    fmpz_poly_set(rop->num, op->num);
    fmpz_poly_scalar_mul_fmpz(rop->den, op->den, y);
    fmpz_poly_q_canonicalise(rop);

    fmpz_clear(y);
}
