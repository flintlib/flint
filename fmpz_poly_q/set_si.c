#include "fmpz_poly_q.h"

void fmpz_poly_q_set_si(fmpz_poly_q_t rop, long op)
{
    fmpz_poly_set_si(rop->num, op);
    fmpz_poly_set_si(rop->den, 1);
}

