/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.1

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <string.h>
#include "test_helpers.h"
#include "fmpz.h"
#include "ecpp.h"

TEST_FUNCTION_START(ecpp_cert_str, state)
{
    slong iter;

    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        fmpz_t n;
        ecpp_cert_t cert, cert2;
        char * str;
        int format = ECPP_CERT_FORMAT_FLINT + n_randint(state, 2);

        fmpz_init(n);
        fmpz_randprime(n, state, 80 + n_randint(state, 120), 0);
        ecpp_cert_init(cert);
        ecpp_cert_init(cert2);

        if (ecpp_prove(cert, n) != 1)
        {
            flint_printf("FAIL: not proved\n");
            flint_abort();
        }
        str = ecpp_cert_get_str(cert, format);
        if (!ecpp_cert_set_str(cert2, str))
        {
            flint_printf("FAIL: cannot read back the certificate (format %d)\n%s\n", format, str);
            flint_abort();
        }
        if (cert2->num != cert->num || !ecpp_verify(cert2, n))
        {
            flint_printf("FAIL: certificate read back does not verify (format %d)\n%s\n", format, str);
            flint_abort();
        }
        /* the PARI format drops D; everything else round-trips */
        if (format == ECPP_CERT_FORMAT_FLINT)
        {
            slong i;
            for (i = 0; i < cert->num; i++)
                if (cert2->steps[i].D != cert->steps[i].D || !fmpz_equal(cert2->steps[i].b, cert->steps[i].b))
                {
                    flint_printf("FAIL: round trip\n");
                    flint_abort();
                }
        }
        flint_free(str);

        /* printing (to a temporary file) in both formats */
        {
            FILE * tmp = tmpfile();
            if (tmp != NULL)
            {
                ecpp_cert_fprint(tmp, cert, ECPP_CERT_FORMAT_FLINT);
                ecpp_cert_fprint(tmp, cert, ECPP_CERT_FORMAT_PARI);
                fclose(tmp);
            }
        }

        /* garbage is rejected */
        if (ecpp_cert_set_str(cert2, "[[12, 3, 4") || ecpp_cert_set_str(cert2, "ecpp certificate, 1 steps\n[0] n = x"))
        {
            flint_printf("FAIL: malformed string accepted\n");
            flint_abort();
        }

        ecpp_cert_clear(cert);
        ecpp_cert_clear(cert2);
        fmpz_clear(n);
    }

    TEST_FUNCTION_END(state);
}
