#include "fmpz_poly_q.h"

void fmpz_poly_q_set(fmpz_poly_q_t rop, const fmpz_poly_q_t op)
{
    if (rop != op)
    {
        fmpz_poly_set(rop->num, op->num);
        fmpz_poly_set(rop->den, op->den);
    }
}

