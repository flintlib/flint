#include "fmpz_poly_q.h"

void fmpz_poly_q_scalar_mul_mpz(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const mpz_t x)
{
    fmpz_t y;

    fmpz_init(y);
    fmpz_set_mpz(y, x);

    fmpz_poly_scalar_mul_fmpz(rop->num, op->num, y);
    fmpz_poly_set(rop->den, op->den);
    fmpz_poly_q_canonicalise(rop);

    fmpz_clear(y);
}
