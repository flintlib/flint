#include "fmpz_poly_q.h"

void fmpz_poly_q_swap(fmpz_poly_q_t op1, fmpz_poly_q_t op2)
{
    if (op1 != op2)
    {
        fmpz_poly_struct *t;

        t        = op1->num;
        op1->num = op2->num;
        op2->num = t;

        t        = op1->den;
        op1->den = op2->den;
        op2->den = t;
    }
}
