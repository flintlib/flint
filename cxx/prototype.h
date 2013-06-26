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
template<>
struct initialization<mpz, mpz>
{
    static void doit(mpz& target, const mpz& source)
    {
        fmpz_init_set(target._data(), source._data());
    }
};

template<class T>
struct initialization<mpz, T,
    typename mp::enable_if<traits::is_unsigned_integer<T> >::type>
{
    static void doit(mpz& target, T source)
    {
        fmpz_init_set_ui(target._data(), source);
    }
};

template<class T>
struct initialization<mpz, T,
    typename mp::enable_if<traits::is_signed_integer<T> >::type>
{
    static void doit(mpz& target, T source)
    {
        fmpz_init(target._data());
        fmpz_set_si(target._data(), source);
    }
};

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

template<>
struct assignment<mpz, mpz>
{
    static void doit(mpz& target, const mpz& source)
    {
        fmpz_set(target._data(), source._data());
    }
};

template<class T>
struct assignment<mpz, T,
    typename mp::enable_if<traits::is_unsigned_integer<T> >::type>
{
    static void doit(mpz& target, T source)
    {
        fmpz_set_ui(target._data(), source);
    }
};

template<class T>
struct assignment<mpz, T,
    typename mp::enable_if<traits::is_signed_integer<T> >::type>
{
    static void doit(mpz& target, T source)
    {
        fmpz_set_si(target._data(), source);
    }
};

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

template<>
struct conversion<slong, mpz>
{
    static slong get(const mpz& from)
    {
        return fmpz_get_si(from._data());
    }
};

template<>
struct conversion<ulong, mpz>
{
    static slong get(const mpz& from)
    {
        return fmpz_get_ui(from._data());
    }
};

template<>
struct conversion<double, mpz>
{
    static double get(const mpz& from)
    {
        return fmpz_get_d(from._data());
    }
};

template<>
struct commutative_binary_expression<mpz, operations::plus, mpz>
{
    typedef mpz return_t;
    static void doit(mpz& to, const mpz& e1, const mpz& e2)
    {
        fmpz_add(to._data(), e1._data(), e2._data());
    }
};
} // rules
} // flint

#endif
