#include "fmpz_poly_q.h"

void fmpz_poly_q_clear(fmpz_poly_q_t rop)
{
    if (rop->num != NULL)
    {
        fmpz_poly_clear(rop->num);
        free(rop->num);
        rop->num = NULL;
    }
    if (rop->den != NULL);
    {
        fmpz_poly_clear(rop->den);
        free(rop->den);
        rop->den = NULL;
    }
}

