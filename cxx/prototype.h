/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Tom Bachmann

******************************************************************************/

#ifndef CXX_PROTOTYPE_H
#define CXX_PROTOTYPE_H

#include <cstdlib>

// TODO
// * non-explicit ctos
// * contexts

#include "cxx/expression.h"
#include "cxx/stdmath.h"

#include "fmpz.h"

namespace flint {
// "concrete" expression classes

template<class Operation, class Data>
class mpz_expression
    : public expression<derived_wrapper<mpz_expression>, Operation, Data>
{
public:
    mpz_expression() {}
    template<class T>
    explicit mpz_expression(const T& t) : mpz_expression::expression(t) {}

    template<class T>
    mpz_expression& operator=(const T& t)
    {
        this->set(t);
        return *this;
    }

protected:
    explicit mpz_expression(const Data& d) : mpz_expression::expression(d) {}

    template<class D, class O, class Da>
    friend class expression;
};
typedef mpz_expression<operations::immediate, fmpz_t> mpz;

namespace rules {
template<>
struct empty_initialization<mpz>
{
    static void doit(mpz& v)
    {
        fmpz_init(v._data());
    }
};

template<>
struct destruction<mpz>
{
    static void doit(mpz& v)
    {
        fmpz_clear(v._data());
    }
};

FLINT_DEFINE_DOIT(initialization, mpz, mpz,
        fmpz_init_set(to._data(), from._data()))

FLINT_DEFINE_DOIT_COND(initialization, mpz, traits::is_unsigned_integer<T>,
        fmpz_init_set_ui(to._data(), from))

FLINT_DEFINE_DOIT_COND(initialization, mpz, traits::is_signed_integer<T>,
        fmpz_init(to._data());fmpz_set_si(to._data(), from))

// TODO should expression automatically deduce this from the assignment
// implementation (or the other way round)?
template<int n>
struct initialization<mpz, char[n]>
{
    static void doit(mpz& target, const char* source)
    {
        fmpz_init(target._data());
        // TODO what about different bases
        fmpz_set_str(target._data(), const_cast<char*>(source), 10);
    }
};

FLINT_DEFINE_DOIT(assignment, mpz, mpz, fmpz_set(to._data(), from._data()))

FLINT_DEFINE_DOIT_COND(assignment, mpz, traits::is_unsigned_integer<T>,
        fmpz_set_ui(to._data(), from))

FLINT_DEFINE_DOIT_COND(assignment, mpz, traits::is_signed_integer<T>,
        fmpz_set_si(to._data(), from))

template<int n>
struct assignment<mpz, char[n]>
{
    static void doit(mpz& target, const char* source)
    {
        fmpz_set_str(target._data(), const_cast<char*>(source), 10);
    }
};

template<>
struct cmp<mpz, mpz>
{
    static int get(const mpz& l, const mpz& r)
    {
        return fmpz_cmp(l._data(), r._data());
    }
};

template<class T>
struct cmp<mpz, T,
    typename mp::enable_if<traits::is_signed_integer<T> >::type>
{
    static int get(const mpz& v, const T& t)
    {
        return fmpz_cmp_si(v._data(), t);
    }
};

template<class T>
struct cmp<mpz, T,
    typename mp::enable_if<traits::is_unsigned_integer<T> >::type>
{
    static int get(const mpz& v, const T& t)
    {
        return fmpz_cmp_ui(v._data(), t);
    }
};

template<>
struct to_string<mpz>
{
    static std::string get(const mpz& v, int base)
    {
        char* str = fmpz_get_str(0, base, v._data());
        std::string res(str);
        std::free(str);
        return res;
    }
};

FLINT_DEFINE_GET(conversion, slong, mpz, fmpz_get_si(from._data()))
FLINT_DEFINE_GET(conversion, ulong, mpz, fmpz_get_ui(from._data()))
FLINT_DEFINE_GET(conversion, double, mpz, fmpz_get_d(from._data()))

FLINT_DEFINE_BINARY_EXPR(plus, mpz,
        fmpz_add(to._data(), e1._data(), e2._data()))

FLINT_DEFINE_CBINARY_EXPR_COND(plus, mpz, traits::is_unsigned_integer<T>,
        fmpz_add_ui(to._data(), e1._data(), e2))

FLINT_DEFINE_BINARY_EXPR(times, mpz,
        fmpz_mul(to._data(), e1._data(), e2._data()))

FLINT_DEFINE_CBINARY_EXPR_COND(times, mpz, traits::is_unsigned_integer<T>,
        fmpz_mul_ui(to._data(), e1._data(), e2))

FLINT_DEFINE_CBINARY_EXPR_COND(times, mpz, traits::is_signed_integer<T>,
        fmpz_mul_si(to._data(), e1._data(), e2))

FLINT_DEFINE_BINARY_EXPR(minus, mpz,
        fmpz_sub(to._data(), e1._data(), e2._data()))

FLINT_DEFINE_BINARY_EXPR_COND(minus, mpz, traits::is_unsigned_integer<T>,
        fmpz_sub_ui(to._data(), e1._data(), e2))

FLINT_DEFINE_BINARY_EXPR(divided_by, mpz,
        fmpz_fdiv_q(to._data(), e1._data(), e2._data()))

FLINT_DEFINE_BINARY_EXPR_COND(divided_by, mpz, traits::is_unsigned_integer<T>,
        fmpz_fdiv_q_ui(to._data(), e1._data(), e2))

FLINT_DEFINE_BINARY_EXPR_COND(divided_by, mpz, traits::is_signed_integer<T>,
        fmpz_fdiv_q_si(to._data(), e1._data(), e2))

// TODO this interpretation of mod is not the same as for builtin types!
FLINT_DEFINE_BINARY_EXPR(modulo, mpz,
        fmpz_mod(to._data(), e1._data(), e2._data()))

FLINT_DEFINE_BINARY_EXPR_COND(modulo, mpz, traits::is_unsigned_integer<T>,
        fmpz_mod_ui(to._data(), e1._data(), e2))

FLINT_DEFINE_UNARY_EXPR(negate, mpz, fmpz_neg(to._data(), from._data()))

// standard math functions (c/f stdmath.h)
FLINT_DEFINE_BINARY_EXPR_COND(pow, mpz, traits::is_unsigned_integer<T>,
        fmpz_pow_ui(to._data(), e1._data(), e2))
FLINT_DEFINE_BINARY_EXPR_COND(root, mpz, traits::fits_into_slong<T>,
        fmpz_root(to._data(), e1._data(), e2))
FLINT_DEFINE_UNARY_EXPR(sqrt, mpz, fmpz_sqrt(to._data(), from._data()))

// TODO ternary rules
} // rules


///////////////////////////////////////////////////////////////////////////
// FUNCTIONS
///////////////////////////////////////////////////////////////////////////
namespace traits {
template<class T, class Enable = void>
struct is_mpz
    : mp::equal_types<typename T::evaluated_t, mpz> { };
template<class T>
struct is_mpz<T, typename mp::disable_if<traits::is_expression<T> >::type>
    : false_ { };
} // traits
namespace mp {
template<class Out, class T1, class T2 = void>
struct enable_all_mpz
    : mp::enable_if<mp::and_<traits::is_mpz<T1>, traits::is_mpz<T2> >, Out> { };
template<class Out, class T>
struct enable_all_mpz<Out, T, void>
    : mp::enable_if<traits::is_mpz<T>, Out> { };
} // mp

// These functions evaluate immediately and do not yield mpzs

template<class T1, class T2>
inline typename mp::enable_all_mpz<bool, T1, T2>::type
divisible(const T1& t1, const T2& t2)
{
    return fmpz_divisible(t1.evaluate()._data(), t2.evaluate()._data());
}
template<class T1, class T2>
inline typename mp::enable_if<mp::and_<
    traits::is_mpz<T1>, traits::fits_into_slong<T2> >, bool>::type
divisible(const T1& t1, const T2& t2)
{
    return fmpz_divisible_si(t1.evaluate()._data(), t2);
}

// TODO should this be lazy?
inline mpz fac(ulong n)
{
    mpz ret;
    fmpz_fac_ui(ret._data(), n);
    return ret;
}

// These functions are evaluated lazily
FLINT_DEFINE_BINOP(rfac)
namespace rules {
FLINT_DEFINE_BINARY_EXPR_COND(rfac, mpz, traits::is_unsigned_integer<T>,
        fmpz_rfac_ui(to._data(), e1._data(), e2))

// TODO many more functions
} // rules

} // flint

#endif
