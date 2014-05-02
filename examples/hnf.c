#include <stdlib.h>
#include <stdio.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"

int main(int argc, char* argv[])
{
    fmpz_mat_t A, B;
    long i, j;

    fmpz_mat_init(A,2,3);
    fmpz_mat_init(B,2,3);
    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            fmpz_set_ui(fmpz_mat_entry(A, i, j), (i+1)*(i + j));

    fmpz_mat_hnf_classical(B,A);

    fmpz_mat_print_pretty(A);
    flint_printf("\n");
    fmpz_mat_print_pretty(B);

    fmpz_mat_clear(A);
    fmpz_mat_clear(B);

    return EXIT_SUCCESS;
}
