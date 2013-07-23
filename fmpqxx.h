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

#include "flintxx/expression.h"
#include "flintxx/flint_classes.h"
#include "flintxx/flint_exception.h"
#include "flintxx/frandxx.h"
#include "fmpzxx.h"

// TODO exhibit this as a specialisation of a generic fraction<fmpzxx>
// TODO summation

namespace flint {
// function "declarations"
FLINT_DEFINE_BINOP(fmpqxx_reconstruct)
FLINT_DEFINE_UNOP(fmpqxx_next_minimal)
FLINT_DEFINE_UNOP(fmpqxx_next_signed_minimal)
FLINT_DEFINE_UNOP(fmpqxx_next_calkin_wilf)
FLINT_DEFINE_UNOP(fmpqxx_next_signed_calkin_wilf)

FLINT_DEFINE_UNOP(height)

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

    // static methods which only make sense with fmpqxx
    FLINTXX_DEFINE_RANDFUNC(fmpq, randbits)
    FLINTXX_DEFINE_RANDFUNC(fmpq, randtest)
    FLINTXX_DEFINE_RANDFUNC(fmpq, randtest_not_zero)

    template<class Vec>
    static fmpqxx_expression from_cfrac(const Vec& v, slong n)
    {
        fmpqxx_expression res;
        res.set_cfrac(v, n);
        return res;
    }

    // TODO does this make more sense as standalone function?
    // TODO in any case we would like it to be lazy (but too many args..)
    template<class Fmpz1, class Fmpz2, class Fmpz3, class Fmpz4>
    static typename mp::enable_all_fmpzxx<
        fmpqxx_expression, Fmpz1, Fmpz2, Fmpz3, Fmpz4>::type reconstruct(
                const Fmpz1& a, const Fmpz2& m, const Fmpz3& N, const Fmpz4& D)
    {
        fmpqxx_expression res;
        // TODO should this throw a different exception type?
        execution_check(fmpq_reconstruct_fmpz_2(res._fmpq(),
                    a.evaluate()._fmpz(), m.evaluate()._fmpz(),
                    N.evaluate()._fmpz(), D.evaluate()._fmpz()),
                "rational reconstruction (v2)", "fmpq");
        return res;
    }

    template<class Fmpz1, class Fmpz2>
    static FLINT_BINOP_ENABLE_RETTYPE(fmpqxx_reconstruct, Fmpz1, Fmpz2)
    reconstruct(const Fmpz1& a, const Fmpz2& m)
    {
        return fmpqxx_reconstruct(a, m);
    }

    // These only make sense with immediates
    fmpzxx_ref num() {return fmpzxx_ref::make(fmpq_numref(_fmpq()));}
    fmpzxx_srcref num() const {return fmpzxx_srcref::make(fmpq_numref(_fmpq()));}
    fmpzxx_ref den() {return fmpzxx_ref::make(fmpq_denref(_fmpq()));}
    fmpzxx_srcref den() const {return fmpzxx_srcref::make(fmpq_denref(_fmpq()));}
    void canonicalise() {fmpq_canonicalise(_fmpq());}
    bool is_canonical() const {return fmpq_is_canonical(_fmpq());}

    template<class Vec>
    void set_cfrac(const Vec& v, slong n)
    {
        fmpq_set_cfrac(this->_fmpq(), v._array(), n);
    }

    // These cause evaluation
    bool is_zero() const {return fmpq_is_zero(this->evaluate()._fmpq());}
    bool is_one() const {return fmpq_is_one(this->evaluate()._fmpq());}
    // TODO make this only work on immediates?
    slong cfrac_bound() const {return fmpq_cfrac_bound(this->evaluate()._fmpq());}

    FLINTXX_DEFINE_MEMBER_UNOP(next_minimal, fmpqxx_next_minimal)
    FLINTXX_DEFINE_MEMBER_UNOP(next_signed_minimal, fmpqxx_next_signed_minimal)
    FLINTXX_DEFINE_MEMBER_UNOP(next_calkin_wilf, fmpqxx_next_calkin_wilf)
    FLINTXX_DEFINE_MEMBER_UNOP(next_signed_calkin_wilf,
            fmpqxx_next_signed_calkin_wilf)
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

// TODO macroize?
namespace traits {
template<class T> struct is_fmpqxx : mp::or_<
     traits::is_T_expr<T, fmpqxx>,
     flint_classes::is_source<fmpqxx, T> > { };
} // traits

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
FLINTXX_DEFINE_SWAP(fmpqxx, fmpq_swap(e1._fmpq(), e2._fmpq()))

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

FLINT_DEFINE_BINARY_EXPR_COND2(modulo, fmpzxx, FMPQXX_COND_S, FMPZXX_COND_S,
        execution_check(fmpq_mod_fmpz(to._fmpz(), e1._fmpq(), e2._fmpz()),
            "modular inversion", "fmpq"))

// TODO macroize?
namespace rdetail {
template<class Fmpq1, class Fmpq2, class T>
void fmpqxx_shift(Fmpq1& to, const Fmpq2& from, T howmuch)
{
    if(howmuch < 0)
        fmpq_div_2exp(to._fmpq(), from._fmpq(), -howmuch);
    else
        fmpq_mul_2exp(to._fmpq(), from._fmpq(), howmuch);
}
} // rdetail
FLINT_DEFINE_BINARY_EXPR_COND2(shift, fmpqxx,
        FMPQXX_COND_S, traits::is_integer,
        rdetail::fmpqxx_shift(to, e1, e2))
} // rules

FLINTXX_DEFINE_TERNARY(fmpqxx,
        fmpq_addmul(to._fmpq(), e1._fmpq(), e2._fmpq()),
        fmpq_submul(to._fmpq(), e1._fmpq(), e2._fmpq()),
        FLINTXX_UNADORNED_MAKETYPES)

// immediate functions
// TODO maybe as a member function?
template<class Fmpq>
inline typename mp::enable_if<traits::is_fmpqxx<Fmpq>, mp_bitcnt_t>::type
height_bits(const Fmpq& f)
{
    return fmpq_height_bits(f.evaluate()._fmpq());
}

// TODO maybe as a member function?
template<class Fmpq1, class Fmpq2, class Vec>
inline typename mp::enable_if<mp::and_<
        traits::is_fmpqxx<Fmpq2>,
        FMPQXX_COND_T<Fmpq1> >,
    int>::type
get_cfrac(Vec& v, Fmpq1& rem, const Fmpq2& x)
{
    return fmpq_get_cfrac(v._array(), rem._fmpq(), x.evaluate()._fmpq(),
            v.size());
}
// TODO also set_cfrac? c/f fmpqxx::set_cfrac ...

// TODO maybe as a member function?
template<class Fmpq>
inline typename mp::enable_if<traits::is_fmpqxx<Fmpq>, int>::type
sgn(const Fmpq& f)
{
    return fmpq_sgn(f.evaluate()._fmpq());
}

namespace rules {
FLINT_DEFINE_UNARY_EXPR_COND(abs_op, fmpqxx, FMPQXX_COND_S,
        fmpq_abs(to._fmpq(), from._fmpq()))
FLINT_DEFINE_UNARY_EXPR_COND(height_op, fmpzxx, FMPQXX_COND_S,
        fmpq_height(to._fmpz(), from._fmpq()))
FLINT_DEFINE_UNARY_EXPR_COND(inv_op, fmpqxx, FMPQXX_COND_S,
        fmpq_inv(to._fmpq(), from._fmpq()))
FLINT_DEFINE_UNARY_EXPR_COND(fmpqxx_next_minimal_op, fmpqxx, FMPQXX_COND_S,
        fmpq_next_minimal(to._fmpq(), from._fmpq()))
FLINT_DEFINE_UNARY_EXPR_COND(fmpqxx_next_signed_minimal_op, fmpqxx, FMPQXX_COND_S,
        fmpq_next_signed_minimal(to._fmpq(), from._fmpq()))
FLINT_DEFINE_UNARY_EXPR_COND(fmpqxx_next_calkin_wilf_op, fmpqxx, FMPQXX_COND_S,
        fmpq_next_calkin_wilf(to._fmpq(), from._fmpq()))
FLINT_DEFINE_UNARY_EXPR_COND(fmpqxx_next_signed_calkin_wilf_op, fmpqxx, FMPQXX_COND_S,
        fmpq_next_signed_calkin_wilf(to._fmpq(), from._fmpq()))

// TODO should this throw a different exception type?
FLINT_DEFINE_BINARY_EXPR_COND2(fmpqxx_reconstruct_op, fmpqxx,
        FMPZXX_COND_S, FMPZXX_COND_S,
        execution_check(fmpq_reconstruct_fmpz(
                    to._fmpq(), e1._fmpz(), e2._fmpz()),
                "rational reconstruction", "fmpq"))

FLINT_DEFINE_BINARY_EXPR_COND2(pow_op, fmpqxx,
        FMPQXX_COND_S, traits::fits_into_slong,
        fmpq_pow_si(to._fmpq(), e1._fmpq(), e2))
}

} // flint

#if 0
// fmpq_vecxx

#include "flintxx/vector.h"

namespace flint {
namespace detail {
struct fmpq_vector_data
{
    long size;
    fmpq* array;

    fmpq_vector_data(long n)
        : size(n), array(_fmpq_vec_init(n)) {}

    ~fmpq_vector_data() {_fmpq_vec_clear(array, size);}

    fmpq_vector_data(const fmpq_vector_data& o)
        : size(o.size), array(_fmpq_vec_init(o.size))
    {
        for(long i = 0;i < size;++i)
            fmpq_set(array + i, o.array + i);
    }

    fmpqxx_ref at(long i) {return fmpqxx_ref::make(array + i);}
    fmpqxx_srcref at(long i) const {return fmpqxx_srcref::make(array + i);}
};
} // detail

typedef vector_expression<
    detail::wrapped_vector_traits<fmpqxx, long, fmpqxx_ref, fmpqxx_srcref>,
    operations::immediate,
    detail::fmpq_vector_data> fmpq_vecxx;

template<>
struct enable_vector_rules<fmpq_vecxx> : mp::false_ { };
}
#endif

#endif
