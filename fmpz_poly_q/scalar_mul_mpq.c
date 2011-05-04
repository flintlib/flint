#include "fmpz_poly_q.h"

void fmpz_poly_q_scalar_mul_mpq(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const mpq_t x)
{
    fmpz_t num, den;

    fmpz_init(num);
    fmpz_init(den);
    fmpz_set_mpz(num, mpq_numref(x));
    fmpz_set_mpz(den, mpq_denref(x));

    fmpz_poly_scalar_mul_fmpz(rop->num, op->num, num);
    fmpz_poly_scalar_mul_fmpz(rop->den, op->den, den);
    fmpz_poly_q_canonicalise(rop);

    fmpz_clear(num);
    fmpz_clear(den);
}
