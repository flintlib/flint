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
FLINT_DEFINE_FOURARY_HERE(fmpqxx_reconstruct) // four argument version
FLINT_DEFINE_UNOP(fmpqxx_next_minimal)
FLINT_DEFINE_UNOP(fmpqxx_next_signed_minimal)
FLINT_DEFINE_UNOP(fmpqxx_next_calkin_wilf)
FLINT_DEFINE_UNOP(fmpqxx_next_signed_calkin_wilf)
FLINT_DEFINE_UNOP(fmpqxx_num)
FLINT_DEFINE_UNOP(fmpqxx_den)

namespace detail {
template<class Fmpq>
struct fmpq_traits {
    typedef FLINT_UNOP_BUILD_RETTYPE(fmpqxx_num, fmpzxx, Fmpq) numreturn_t;
    typedef FLINT_UNOP_BUILD_RETTYPE(fmpqxx_den, fmpzxx, Fmpq) denreturn_t;
    typedef numreturn_t cnumreturn_t;
    typedef denreturn_t cdenreturn_t;
    static numreturn_t num(const Fmpq& f) {return fmpqxx_num(f);}
    static denreturn_t den(const Fmpq& f) {return fmpqxx_den(f);}
};
}

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
    template<class Fmpz1, class Fmpz2>
    static FLINT_BINOP_ENABLE_RETTYPE(fmpqxx_reconstruct, Fmpz1, Fmpz2)
    reconstruct(const Fmpz1& a, const Fmpz2& m)
    {
        return fmpqxx_reconstruct(a, m);
    }
    template<class Fmpz1, class Fmpz2, class Fmpz3, class Fmpz4>
    static FLINT_FOURARY_ENABLE_RETTYPE(fmpqxx_reconstruct,
            Fmpz1, Fmpz2, Fmpz3, Fmpz4)
    reconstruct(const Fmpz1& a, const Fmpz2& m, const Fmpz3& N, const Fmpz4& D)
    {
        return fmpqxx_reconstruct(a, m, N, D);
    }

    template<class F1, class F2>
    void set_frac(const F1& f1, const F2& f2)
    {
        num() = f1;
        den() = f2;
        canonicalise();
    }
    template<class F1, class F2>
    static fmpqxx_expression frac(const F1& f1, const F2& f2)
    {
        fmpqxx_expression res;
        res.set_frac(f1, f2);
        return res;
    }

    template<class T>
    void set_integer(const T& t)
    {
        num() = t;
        den() = 1u;
    }
    template<class T>
    static fmpqxx_expression integer(const T& t)
    {
        fmpqxx_expression res;
        res.set_integer(t);
        return res;
    }

    static fmpqxx_expression zero(){return fmpqxx_expression();}
    static fmpqxx_expression one()
    {
        fmpqxx_expression res;
        res.set_one();
        return res;
    }

    // These only make sense with immediates
    void canonicalise() {fmpq_canonicalise(_fmpq());}
    bool is_canonical() const {return fmpq_is_canonical(_fmpq());}
    void set_zero() {fmpq_zero(_fmpq());}
    void set_one() {fmpq_one(_fmpq());}

    template<class Vec>
    void set_cfrac(const Vec& v, slong n)
    {
        fmpq_set_cfrac(this->_fmpq(), v._array(), n);
    }

    // Numerator and denominator access
    typedef detail::fmpq_traits<fmpqxx_expression> traits_t;
    typename traits_t::numreturn_t num() {return traits_t::num(*this);}
    typename traits_t::cnumreturn_t num() const {return traits_t::num(*this);}
    typename traits_t::denreturn_t den() {return traits_t::den(*this);}
    typename traits_t::cdenreturn_t den() const
        {return traits_t::den(*this);}

    // These cause evaluation
    bool is_zero() const {return fmpq_is_zero(this->evaluate()._fmpq());}
    bool is_one() const {return fmpq_is_one(this->evaluate()._fmpq());}
    // TODO make this only work on immediates?
    slong cfrac_bound() const {return fmpq_cfrac_bound(this->evaluate()._fmpq());}
    int sgn() const {return fmpq_sgn(this->evaluate()._fmpq());}
    mp_bitcnt_t height_bits() const
        {return fmpq_height_bits(this->evaluate()._fmpq());}

    FLINTXX_DEFINE_MEMBER_UNOP_(next_minimal, fmpqxx_next_minimal)
    FLINTXX_DEFINE_MEMBER_UNOP_(next_signed_minimal, fmpqxx_next_signed_minimal)
    FLINTXX_DEFINE_MEMBER_UNOP_(next_calkin_wilf, fmpqxx_next_calkin_wilf)
    FLINTXX_DEFINE_MEMBER_UNOP_(next_signed_calkin_wilf,
            fmpqxx_next_signed_calkin_wilf)

    // forwarded member functions
    FLINTXX_DEFINE_MEMBER_UNOP(abs)
    FLINTXX_DEFINE_MEMBER_UNOP(inv)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(fmpzxx, height)
    FLINTXX_DEFINE_MEMBER_BINOP(pow)
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
template<>
struct fmpq_traits<fmpqxx_srcref>
{
    typedef fmpzxx_srcref numreturn_t;
    typedef fmpzxx_srcref cnumreturn_t;
    typedef fmpzxx_srcref denreturn_t;
    typedef fmpzxx_srcref cdenreturn_t;
    template<class T>
    static cnumreturn_t num(T f)
        {return cnumreturn_t::make(fmpq_numref(f._fmpq()));}
    template<class T>
    static cnumreturn_t den(T f)
        {return cnumreturn_t::make(fmpq_denref(f._fmpq()));}
};
template<>
struct fmpq_traits<fmpqxx_ref>
{
    typedef fmpzxx_ref numreturn_t;
    typedef fmpzxx_ref denreturn_t;
    typedef fmpzxx_ref cnumreturn_t;
    typedef fmpzxx_ref cdenreturn_t;
    template<class T>
    static cnumreturn_t num(T f)
        {return cnumreturn_t::make(fmpq_numref(f._fmpq()));}
    template<class T>
    static cnumreturn_t den(T f)
        {return cnumreturn_t::make(fmpq_denref(f._fmpq()));}
};
template<> struct fmpq_traits<fmpqxx>
{
    typedef fmpzxx_ref numreturn_t;
    typedef fmpzxx_ref denreturn_t;
    typedef fmpzxx_srcref cnumreturn_t;
    typedef fmpzxx_srcref cdenreturn_t;
    template<class T>
    static cnumreturn_t num(const T& f)
        {return cnumreturn_t::make(fmpq_numref(f._fmpq()));}
    template<class T>
    static cnumreturn_t den(const T& f)
        {return cnumreturn_t::make(fmpq_denref(f._fmpq()));}
    template<class T>
    static numreturn_t num(T& f)
        {return numreturn_t::make(fmpq_numref(f._fmpq()));}
    template<class T>
    static numreturn_t den(T& f)
        {return numreturn_t::make(fmpq_denref(f._fmpq()));}
};

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

    fmpq_data(fmpqxx_srcref r)
    {
        fmpq_init(inner);
        fmpq_set(inner, r._fmpq());
    }

    fmpq_data(fmpzxx_srcref num, fmpzxx_srcref den)
    {
        fmpq_init(inner);
        fmpq_set_fmpz_frac(inner, num._fmpz(), den._fmpz());
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

FLINT_DEFINE_PRINT_COND(FMPQXX_COND_S, (fmpq_fprint(to, from._fmpq()), 1))

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
template<class Fmpq>
inline typename mp::enable_if<traits::is_fmpqxx<Fmpq>, mp_bitcnt_t>::type
height_bits(const Fmpq& f)
{
    return f.height_bits();
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

template<class Fmpq>
inline typename mp::enable_if<traits::is_fmpqxx<Fmpq>, int>::type
sgn(const Fmpq& f)
{
    return f.sgn();
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

FLINT_DEFINE_UNARY_EXPR_COND(fmpqxx_num_op, fmpzxx, FMPQXX_COND_S,
        fmpz_set(to._fmpz(), fmpq_numref(from._fmpq())))
FLINT_DEFINE_UNARY_EXPR_COND(fmpqxx_den_op, fmpzxx, FMPQXX_COND_S,
        fmpz_set(to._fmpz(), fmpq_denref(from._fmpq())))

// TODO should this throw a different exception type?
FLINT_DEFINE_BINARY_EXPR_COND2(fmpqxx_reconstruct_op, fmpqxx,
        FMPZXX_COND_S, FMPZXX_COND_S,
        execution_check(fmpq_reconstruct_fmpz(
                    to._fmpq(), e1._fmpz(), e2._fmpz()),
                "rational reconstruction", "fmpq"))
FLINT_DEFINE_FOURARY_EXPR_COND4(fmpqxx_reconstruct_op, fmpqxx,
        FMPZXX_COND_S, FMPZXX_COND_S, FMPZXX_COND_S, FMPZXX_COND_S,
        execution_check(fmpq_reconstruct_fmpz_2(
                    to._fmpq(), e1._fmpz(), e2._fmpz(), e3._fmpz(), e4._fmpz()),
                "rational reconstruction (v2)", "fmpq"))


FLINT_DEFINE_BINARY_EXPR_COND2(pow_op, fmpqxx,
        FMPQXX_COND_S, traits::fits_into_slong,
        fmpq_pow_si(to._fmpq(), e1._fmpq(), e2))
}

} // flint

// fmpq_vecxx

#include "flintxx/vector.h"

namespace flint {
namespace detail {
struct fmpq_vector_data
{
    slong size;
    fmpq* array;

    fmpq_vector_data(slong n)
        : size(n), array(_fmpq_vec_init(n)) {}

    ~fmpq_vector_data() {_fmpq_vec_clear(array, size);}

    fmpq_vector_data(const fmpq_vector_data& o)
        : size(o.size), array(_fmpq_vec_init(o.size))
    {
        for(slong i = 0;i < size;++i)
            fmpq_set(array + i, o.array + i);
    }

    fmpqxx_ref at(slong i) {return fmpqxx_ref::make(array + i);}
    fmpqxx_srcref at(slong i) const {return fmpqxx_srcref::make(array + i);}
};
} // detail

typedef vector_expression<
    detail::wrapped_vector_traits<fmpqxx, slong, fmpqxx_ref, fmpqxx_srcref, fmpq>,
    operations::immediate,
    detail::fmpq_vector_data> fmpq_vecxx;

template<>
struct enable_vector_rules<fmpq_vecxx> : mp::false_ { };

namespace detail {
inline bool fmpq_vec_equal(const fmpq* v1, const fmpq* v2, slong n)
{
    for(slong i = 0;i < n;++i)
        if(!fmpq_equal(v1+i, v2+i))
            return false;
    return true;
}
}
namespace rules {
// TODO hack to make code look like references are implemented
template<class T> struct FMPQ_VECXX_COND_S : mp::equal_types<T, fmpq_vecxx> { };
#define FMPQ_VECXX_COND_T FMPQ_VECXX_COND_S

// TODO references
FLINT_DEFINE_GET(equals, bool, fmpq_vecxx,
        e1.size() == e2.size()
        && detail::fmpq_vec_equal(e1._data().array, e2._data().array, e1.size()))
} // rules
} // flint

#endif
