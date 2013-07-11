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

#ifndef CXX_FMPQXX_H
#define CXX_FMPQXX_H

#include <cstdlib>

#include "fmpq.h"

#include "cxx/expression.h"
#include "cxx/fmpzxx.h"

// TODO exhibit this as a specialisation of a generic fraction<fmpzxx>

namespace flint {
template<class Operation, class Data>
class fmpqxx_expression
    : public expression<derived_wrapper<fmpqxx_expression>, Operation, Data>
{
public:
    typedef expression<derived_wrapper< ::flint::fmpqxx_expression>,
              Operation, Data> base_t;

    fmpqxx_expression() {}

    template<class T>
    explicit fmpqxx_expression(const T& t)
        : base_t(t) {}

    template<class T, class U>
    fmpqxx_expression(const T& t, const U& u)
        : base_t(t, u) {}

    template<class T>
    fmpqxx_expression& operator=(const T& t)
    {
        this->set(t);
        return *this;
    }

    // These only make sense with immediates
    fmpq_t& _fmpq() {return this->_data().f;}
    const fmpq_t& _fmpq() const {return this->_data().f;}
    fmpzxx_ref num() {return fmpzxx_ref::make(fmpq_numref(_fmpq()));}
    fmpzxx_cref num() const {return fmpzxx_cref::make(fmpq_numref(_fmpq()));}
    fmpzxx_ref den() {return fmpzxx_ref::make(fmpq_denref(_fmpq()));}
    fmpzxx_cref den() const {return fmpzxx_cref::make(fmpq_denref(_fmpq()));}
    void canonicalise() {fmpq_canonicalise(_fmpq());}
    bool is_canonical() const {return fmpq_is_canonical(_fmpq());}

protected:
    explicit fmpqxx_expression(const Data& d) : base_t(d) {}

    template<class D, class O, class Da>
    friend class expression;
};

namespace detail {
struct fmpq_data
{
    fmpq_t f;

    fmpq_data() {fmpq_init(f);}
    ~fmpq_data() {fmpq_clear(f);}

    fmpq_data(const fmpq_data& o)
    {
        fmpq_init(f);
        fmpq_set(f, o.f);
    }

    fmpq_data(fmpzxx_cref num, fmpzxx_cref den)
    {
        fmpq_init(f);
        fmpq_set_fmpz_frac(f, num._fmpz(), den._fmpz());
    }

    template<class T, class U>
    fmpq_data(T num, U den,
            typename mp::enable_if<traits::fits_into_slong<T> >::type* = 0,
            typename mp::enable_if<traits::is_unsigned_integer<U> >::type* = 0)
    {
        fmpq_init(f);
        fmpq_set_si(f, num, den);
    }
};
} // detail

typedef fmpqxx_expression<operations::immediate, detail::fmpq_data> fmpqxx;

namespace rules {
FLINT_DEFINE_DOIT(assignment, fmpqxx, fmpqxx,
        fmpq_set(to._fmpq(), from._fmpq()))
FLINT_DEFINE_DOIT_COND(assignment, fmpqxx, traits::fits_into_slong<T>, 
        fmpq_set_si(to._fmpq(), from, 1))
// TODO mpq, mpfr?

template<>
struct to_string<fmpqxx>
{
    static std::string get(const fmpqxx& p, int base)
    {
        char* str = fmpq_get_str(0, base, p._fmpq());
        std::string res(str);
        std::free(str);
        return res;
    }
};

FLINT_DEFINE_GET2(cmp, int, fmpqxx, fmpqxx, fmpq_cmp(e1._fmpq(), e2._fmpq()))

FLINT_DEFINE_CBINARY_EXPR(plus, fmpqxx,
        fmpq_add(to._fmpq(), e1._fmpq(), e2._fmpq()))
FLINT_DEFINE_CBINARY_EXPR(minus, fmpqxx,
        fmpq_sub(to._fmpq(), e1._fmpq(), e2._fmpq()))
FLINT_DEFINE_CBINARY_EXPR(times, fmpqxx,
        fmpq_mul(to._fmpq(), e1._fmpq(), e2._fmpq()))
FLINT_DEFINE_CBINARY_EXPR(divided_by, fmpqxx,
        fmpq_div(to._fmpq(), e1._fmpq(), e2._fmpq()))
FLINT_DEFINE_CBINARY_EXPR_COND(times, fmpqxx, fmpzxx_traits::is_source<T>,
        fmpq_mul_fmpz(to._fmpq(), e1._fmpq(), e2._fmpz()))
FLINT_DEFINE_BINARY_EXPR_COND(divided_by, fmpqxx, fmpzxx_traits::is_source<T>,
        fmpq_div_fmpz(to._fmpq(), e1._fmpq(), e2._fmpz()))

FLINT_DEFINE_UNARY_EXPR(negate, fmpqxx, fmpq_neg(to._fmpq(), from._fmpq()))

// TODO addmul

// TODO functions
} // rules
} // flint

#endif
