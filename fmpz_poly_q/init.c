#include "fmpz_poly_q.h"

void fmpz_poly_q_init(fmpz_poly_q_t rop)
{
    rop->num = malloc(sizeof(fmpz_poly_struct));
    rop->den = malloc(sizeof(fmpz_poly_struct));
    fmpz_poly_init(rop->num);
    fmpz_poly_init(rop->den);
    fmpz_poly_set_si(rop->den, 1);
}

