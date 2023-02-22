#ifndef EXTRACTING_SAMPLE
/*
    Copyright (C) 2013 Tom Bachmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/


// This file does not actually meant compile. It is used to check compiler error
// outputs. NB: File extension .cc to stop the build system from trying to compile
// this.

#include "fmpzxx.h"
#ifdef TEST_PADICXX_FORGET_EVAL
#include "padicxx.h"
#endif

struct newtype { };

using namespace flint;

int
main()
{
#endif // EXTRACTING_SAMPLE
#ifdef TEST_FMPZXX_INIT_WRONG
    {
        newtype n;
        fmpzxx a(n);
    }
#endif
#ifdef TEST_FMPZXX_INIT_2
    {
        fmpzxx a(3, 4);
    }
#endif
#ifdef TEST_FMPZXX_ASSIGN_WRONG
    {
        fmpzxx a;
        newtype n;
        a = n;
    }
#endif
#ifdef TEST_FMPZXX_CONVERT_WRONG
    {
        fmpzxx a;
        a.to<newtype>();
    }
#endif
#ifdef TEST_FMPZXX_REF_INIT_WRONG_1
    {
        const fmpzxx a;
        fmpzxx_ref ar(a);
    }
#endif
#ifdef TEST_FMPZXX_REF_INIT_WRONG_2
    {
        const fmpzxx a;
        fmpzxx_srcref asr(a);
        fmpzxx_ref ar(asr);
    }
#endif
#ifdef TEST_FMPZXX_SRCREF_ASSIGN
    {
        fmpzxx a;
        fmpzxx_srcref b(a);
        b = a;
    }
#endif
#ifdef TEST_FMPZXX_ARITH_WRONG
    {
        fmpzxx a;
        newtype n;
        a + n;
    }
#endif
#ifdef TEST_FMPZXX_ARITH_WRONG_DEEP
    {
        fmpzxx a;
        newtype n;
        a + (a*a + (a / n) + a)*a;
    }
#endif
#ifdef TEST_FMPZXX_ARITHFUNC_WRONG_NARGS
    {
        fmpzxx a;
        gcd(a);
    }
#endif
#ifdef TEST_FMPZXX_ARITHFUNC_WRONG_TYPE
    {
        newtype n;
        dlog(n);
    }
#endif
#ifdef TEST_FMPZXX_ARITHFUNC_WRONG_TYPE2
    {
        fmpzxx a;
        fac(a);
    }
#endif
#ifdef TEST_PADICXX_FORGET_EVAL
    {
        padicxx_ctx ctx(fmpzxx(5), 1, 2, PADIC_TERSE);
        padicxx a(ctx);
        (a + a).unit();
    }
#endif
#ifndef EXTRACTING_SAMPLE
    // TODO: assignment through tuple with wrong type
}
#endif
