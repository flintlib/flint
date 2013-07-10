#include <iostream>
#include <sstream>
#include <string>

#include "cxx/padicxx.h"
#include "cxx/test/helpers.h"

using namespace flint;

void
test_init()
{
    padicxx_ctx ctx(fmpzxx(5), 10, 20, PADIC_TERSE);
    tassert(ctx.get_p() == 5);

    padicxx a(ctx, 20);
    tassert(&a.get_ctx() == &ctx);
    tassert(a.prec() == 20);

    padicxx c(a);
    tassert(a == c);

    padicxx b(ctx, 30);
    tassert(&(a + b).estimate_ctx() == &ctx);
    tassert((a + b).prec() == 30);

    tassert((a + b).create_temporary()._prec() == 30);
    padicxx d((a + b).create_temporary());
    tassert(&d.get_ctx() == &ctx);
    tassert(d.prec() == 30);

    padicxx e(a + b);
    tassert(e.prec() == 30);
}

void
test_assignment()
{
    padicxx_ctx ctx(fmpzxx(5), 10, 20, PADIC_TERSE);
    padicxx a(ctx, 20), b(ctx, 20);
    fmpzxx c(17);

    a = 17; tassert(a != b);
    b = 17; tassert(a == b);
    b = 0; tassert(a != b);
    b = 17ul; tassert(a == b);
    b = 0; b = c; tassert(a == b);
    b = 0; b = fmpzxx_ref(c); tassert(a == b);
    b = 0; b = fmpzxx_cref(c); tassert(a == b);
    // TODO fmpq
}

void
test_conversion()
{
    padicxx_ctx ctx(fmpzxx(5), 10, 20, PADIC_TERSE);
    padicxx_ctx ctx2(fmpzxx(5), 10, 20, PADIC_VAL_UNIT);

    padicxx a(ctx), b(ctx2);
    a = 15; b = 15;
    tassert(a.to_string() == "15");
    tassert(b.to_string() == "3*5");

    tassert(a.to<fmpzxx>() == 15);
    // TODO fmpq
}

void
test_arithmetic()
{
    // TODO
}

int
main()
{
    std::cout << "padicxx....";

    test_init();
    test_assignment();
    test_conversion();
    test_arithmetic();

    std::cout << "PASS" << std::endl;
    return 0;
}
