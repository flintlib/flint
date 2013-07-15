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
#include "cxx/flint_classes.h"
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

    FLINTXX_DEFINE_BASICS(fmpqxx_expression)
    FLINTXX_DEFINE_CTORS(fmpqxx_expression)
    FLINTXX_DEFINE_C_REF(fmpqxx_expression, fmpq, _fmpq)

    // These only make sense with immediates
    fmpzxx_ref num() {return fmpzxx_ref::make(fmpq_numref(_fmpq()));}
    fmpzxx_srcref num() const {return fmpzxx_srcref::make(fmpq_numref(_fmpq()));}
    fmpzxx_ref den() {return fmpzxx_ref::make(fmpq_denref(_fmpq()));}
    fmpzxx_srcref den() const {return fmpzxx_srcref::make(fmpq_denref(_fmpq()));}
    void canonicalise() {fmpq_canonicalise(_fmpq());}
    bool is_canonical() const {return fmpq_is_canonical(_fmpq());}
};

namespace detail {
struct fmpq_data;
}

typedef fmpqxx_expression<operations::immediate, detail::fmpq_data> fmpqxx;
typedef fmpqxx_expression<operations::immediate,
            flint_classes::ref_data<fmpqxx, fmpq> > fmpqxx_ref;
typedef fmpqxx_expression<operations::immediate,
            flint_classes::srcref_data<fmpqxx, fmpqxx_ref, fmpq> > fmpqxx_srcref;

namespace detail {
struct fmpq_data
{
    fmpq_t inner;
    typedef fmpq_t& data_ref_t;
    typedef const fmpq_t& data_srcref_t;

    fmpq_data() {fmpq_init(inner);}
    ~fmpq_data() {fmpq_clear(inner);}

    fmpq_data(const fmpq_data& o)
    {
        fmpq_init(inner);
        fmpq_set(inner, o.inner);
    }

    fmpq_data(fmpzxx_srcref num, fmpzxx_srcref den)
    {
        fmpq_init(inner);
        fmpq_set_fmpz_frac(inner, num._fmpz(), den._fmpz());
    }

    fmpq_data(fmpqxx_srcref r)
    {
        fmpq_init(inner);
        fmpq_set(inner, r._fmpq());
    }

    template<class T, class U>
    fmpq_data(T num, U den,
            typename mp::enable_if<traits::fits_into_slong<T> >::type* = 0,
            typename mp::enable_if<traits::is_unsigned_integer<U> >::type* = 0)
    {
        fmpq_init(inner);
        fmpq_set_si(inner, num, den);
    }
};
} // detail

namespace rules {
#define FMPQXX_COND_S FLINTXX_COND_S(fmpqxx)
#define FMPQXX_COND_T FLINTXX_COND_T(fmpqxx)

FLINT_DEFINE_DOIT_COND2(assignment, FMPQXX_COND_T, FMPQXX_COND_S,
        fmpq_set(to._fmpq(), from._fmpq()))
FLINT_DEFINE_DOIT_COND2(assignment, FMPQXX_COND_T, traits::fits_into_slong, 
        fmpq_set_si(to._fmpq(), from, 1))
// TODO mpq, mpfr?

FLINTXX_DEFINE_TO_STR(fmpqxx, fmpq_get_str(0,  base, from._fmpq()))
FLINTXX_DEFINE_CMP(fmpqxx, fmpq_cmp(e1._fmpq(), e2._fmpq()))

FLINT_DEFINE_CBINARY_EXPR_COND2(plus, fmpqxx, FMPQXX_COND_S, FMPQXX_COND_S,
        fmpq_add(to._fmpq(), e1._fmpq(), e2._fmpq()))
FLINT_DEFINE_CBINARY_EXPR_COND2(minus, fmpqxx, FMPQXX_COND_S, FMPQXX_COND_S,
        fmpq_sub(to._fmpq(), e1._fmpq(), e2._fmpq()))
FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpqxx, FMPQXX_COND_S, FMPQXX_COND_S,
        fmpq_mul(to._fmpq(), e1._fmpq(), e2._fmpq()))
FLINT_DEFINE_CBINARY_EXPR_COND2(divided_by, fmpqxx, FMPQXX_COND_S,
        FMPQXX_COND_S, fmpq_div(to._fmpq(), e1._fmpq(), e2._fmpq()))
FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpqxx, FMPQXX_COND_S, FMPZXX_COND_S,
        fmpq_mul_fmpz(to._fmpq(), e1._fmpq(), e2._fmpz()))
FLINT_DEFINE_BINARY_EXPR_COND2(divided_by, fmpqxx, FMPQXX_COND_S, FMPZXX_COND_S,
        fmpq_div_fmpz(to._fmpq(), e1._fmpq(), e2._fmpz()))

FLINT_DEFINE_UNARY_EXPR_COND(negate, fmpqxx, FMPQXX_COND_S,
        fmpq_neg(to._fmpq(), from._fmpq()))

// TODO addmul

// TODO functions
} // rules
} // flint

#endif
