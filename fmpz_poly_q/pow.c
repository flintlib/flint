#include "fmpz_poly_q.h"

void fmpz_poly_q_pow(fmpz_poly_q_t rop, const fmpz_poly_q_t op, unsigned long exp)
{
    if (exp == 0)
    {
        fmpz_poly_q_one(rop);
    }
    else
    {
        fmpz_poly_pow(rop->num, op->num, exp);
        fmpz_poly_pow(rop->den, op->den, exp);
    }
}
