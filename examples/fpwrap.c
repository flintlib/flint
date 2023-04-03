/* This file is public domain. Author: Fredrik Johansson. */

#include "arb_fpwrap.h"

int main()
{
    double x, y;
    complex_double cx, cy;
    int flags = 0;    /* default options */

    x = 2.0;
    cx.real = 0.5;
    cx.imag = 123.0;

    arb_fpwrap_double_zeta(&y, x, flags);
    arb_fpwrap_cdouble_zeta(&cy, cx, flags);

    printf("zeta(%g) = %.16g\n", x, y);
    printf("zeta(%g + %gi) = %.16g + %.16gi\n", cx.real, cx.imag, cy.real, cy.imag);

    flint_cleanup();
    return 0;
}

