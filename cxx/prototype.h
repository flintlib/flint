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
struct initialization<mpz, mpz>
{
    static void doit(mpz& target, const mpz& source)
    {
        fmpz_init_set(target._data(), source._data());
    }
};

template<class T>
struct initialization<mpz, T, typename mp::enable_if<traits::is_unsigned_integer<T> >::type>
{
    static void doit(mpz& target, T source)
    {
        fmpz_init_set_ui(target._data(), source);
    }
};

template<class T>
struct initialization<mpz, T, typename mp::enable_if<traits::is_signed_integer<T> >::type>
{
    static void doit(mpz& target, T source)
    {
        fmpz_init(target._data());
        fmpz_set_si(target._data(), source);
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
struct assignment<mpz, T, typename mp::enable_if<traits::is_unsigned_integer<T> >::type>
{
    static void doit(mpz& target, T source)
    {
        fmpz_set_ui(target._data(), source);
    }
};

template<class T>
struct assignment<mpz, T, typename mp::enable_if<traits::is_signed_integer<T> >::type>
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
struct print<mpz>
{
    static void doit(const mpz& v, std::ostream& o)
    {
        char* str = fmpz_get_str(0, 10, v._data());
        o << str;
        std::free(str);
    }
};

template<bool result_is_temporary>
struct evaluation<mpz_expression<operations::plus, tuple<const mpz&, tuple<const mpz&, empty_tuple> > >, result_is_temporary, 2>
{
    typedef mpz return_t;
    typedef empty_tuple temporaries_t;
    static void doit(const mpz_expression<operations::plus, tuple<const mpz&, tuple<const mpz&, empty_tuple> > >& input, temporaries_t temps, return_t* output)
    {
        fmpz_add(output->_data(), input._data().first()._data(), input._data().second()._data());
    }
};

// TODO make more generic
template<bool result_is_temporary, class Op, class Data1, class Data2>
struct evaluation<
    mpz_expression<Op, tuple<Data1, tuple<Data2, empty_tuple> > >,
    result_is_temporary, 1,
    typename mp::enable_if<mp::and_<mp::not_<traits::is_immediate<Data1> >,
                                    mp::not_<traits::is_immediate<Data1> > > >::type>
{
    typedef mpz return_t;
    typedef typename mp::find_evaluation<Data1, true>::type ev1_t;
    typedef typename mp::find_evaluation<Data2, true>::type ev2_t;
    typedef mp::concat_tuple<tuple<mpz*, typename ev1_t::temporaries_t>, tuple<mpz*, typename ev2_t::temporaries_t> > concater;
    typedef typename concater::type temporaries_t;

    static void doit(const mpz_expression<Op, tuple<Data1, tuple<Data2, empty_tuple> > >& input, temporaries_t temps, return_t* output)
    {
        tuple<mpz*, typename ev1_t::temporaries_t> temps1 = concater::get_first(temps);
        tuple<mpz*, typename ev2_t::temporaries_t> temps2 = concater::get_second(temps);
        ev1_t::doit(input._data().first(), temps1.tail, temps1.head);
        ev2_t::doit(input._data().second(), temps2.tail, temps2.head);
        *output = *temps1.head + *temps2.head;
    }
};
} // rules
} // flint

#endif
