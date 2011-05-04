#include "fmpz_poly_q.h"

void fmpz_poly_q_addmul(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)
{
    fmpz_poly_q_t temp;

    fmpz_poly_q_init(temp);
    fmpz_poly_q_mul(temp, op1, op2);
    fmpz_poly_q_add(rop, rop, temp);
    fmpz_poly_q_clear(temp);
}
