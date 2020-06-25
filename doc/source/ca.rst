.. ca:

**ca.h** -- exact real and complex numbers
===============================================================================

A :type:`ca_t` represents a real or complex number in a form suitable
for exact field arithmetic or comparison.
Exceptionally, a :type:`ca_t` may represent a special nonnumerical value,
such as an infinity.

Introduction: numbers
-------------------------------------------------------------------------------

A *Calcium number* is a real or complex number
represented as an element of a formal field `K = \mathbb{Q}(a_1, \ldots, a_n)`
where the symbols `a_k` denote fixed algebraic or transcendental numbers.
For example, `e^{-2 \pi} - 3 i` may be represented as
`(1 - 3 a_2^2 a_1) / a_2^2` in the field `\mathbb{Q}(a_1,a_2)` with
`a_1 = i, a_2 = e^{\pi}`.
Extension elements `a_k` can be algebraic numbers represented
in canonical form by :type:`qqbar_t` instances (example: `i`),
symbolic constants (example: `\pi`),
or symbolic functions applied to :type:`ca_t` arguments
(example: `e^{\pi} = \exp(x)`, where `x = 1 \cdot \pi`
is a :type:`ca_t` representing an element of `\mathbb{Q}(\pi)`).
The user does not need to construct formal fields explicitly:
operations on Calcium numbers generate and cache fields automatically as needed
to express the results.

This representation is not canonical (in general). There can be
many representations for the same number, depending on
the choice of formal field *K*. Even within a fixed field *K*,
a number can have many representations if there are algebraic
relations between the extension elements.
Two numbers *x* and *y* can be tested for inequality using numerical
evaluation; to test for equality, it may be necessary to
eliminate redundant extension elements.
One of the central goals of Calcium will be to
implement heuristics for such elimination.

Introduction: special values
-------------------------------------------------------------------------------

In order to provide a closed arithmetic system and express limiting
cases of operators and special functions, a :type:`ca_t` can hold
any of the following special values besides ordinary numbers:

* *Unsigned infinity*, a formal object `{\tilde \infty}` representing
  the value of `1 / 0`. More generally, this is the
  value of meromorphic functions at poles.

* *Signed infinity*, a formal object `a \cdot \infty` where the sign `a`
  is a Calcium number with `|a| = 1`.
  The most common values are `+\infty, -\infty, +i \infty, -i \infty`.
  Signed infinities are used to denote directional limits and logarithmic
  singularities (for example, `\log(0) = -\infty`).

* *Undefined*, a formal object representing the value of indeterminate
  forms such as `0 / 0` and essential singularities such as 
  `\exp(\tilde \infty)`, where a number or infinity would not make sense
  as an answer.

* *Unknown*, a meta-value used to signal that the actual desired value
  could not be computed, either because Calcium does not (yet) have a data
  structure or algorithm for that case, or because doing so would be
  unreasonably expensive. This occurs, for example, if Calcium performs
  a division and is unable to decide whether the result is a number,
  unsigned infinity or undefined (because testing for zero fails).
  Wrappers may want to check output variables for
  *Unknown* and throw an exception (e.g. *NotImplementedError* in Python).

The distinction between *Calcium numbers* (which must represent
elements of `\mathbb{C}`) and the different kinds of nonnumerical values
(infinities, Undefined or Unknown) is essential. Nonnumerical values may
not be used as field extension elements `a_k`, and the denominator of
a formal field element must always represent a nonzero complex number.
Accordingly, for any given Calcium value *x* that is not *Unknown*,
it is exactly known whether *x* represents A) a number, B) unsigned infinity,
C) a signed infinity, or D) Undefined.


Number objects
-------------------------------------------------------------------------------

For all types, a *type_t* is defined as an array of length one of type
*type_struct*, permitting a *type_t* to be passed by reference.

.. type:: ca_struct

.. type:: ca_t

    A :type:`ca_t` contains an index to a field *K*, and
    data representing an element *x* of *K*.
    The data is either an inline rational number (:type:`fmpq_t`), an
    inline Antic number field element (:type:`nf_elem_t`)
    when *K* is an absolute algebraic number field `\mathbb{Q}(a)`,
    or a pointer to a heap-allocated :type:`fmpz_mpoly_q_t` representing
    an element of a generic field `\mathbb{Q}(a_1,\ldots,a_n)`.
    Special values are encoded using magic bits in the field index.

.. type:: ca_ptr

   Alias for ``ca_struct *``, used for vectors of numbers.

.. type:: ca_srcptr

   Alias for ``const ca_struct *``, used for vectors of numbers
   when passed as constant input to functions.


Context objects
-------------------------------------------------------------------------------

.. type:: ca_ctx_struct

.. type:: ca_ctx_t

    A :type:`ca_ctx_t` context object holds a cache of fields *K* and
    constituent extension elements `a_k`.
    The field index in an individual :type:`ca_t` instance represents
    a shallow reference to the object defining the field *K* within the
    context object, so creating many elements of the same field is cheap.

    Since context objects are mutable (and may be mutated even when
    performing read-only operations on :type:`ca_t` instances), they must not
    be accessed simultaneously by different threads: in multithreaded
    environments, the user must use a separate context object for each thread.

.. function:: void ca_ctx_init(ca_ctx_t ctx)

    Initializes the context object *ctx* for use.
    Any evaluation options stored in the context object
    are set to default values.

.. function:: void ca_ctx_clear(ca_ctx_t ctx)

    Clears the context object *ctx*, freeing any memory allocated internally.
    This function should only be called after all :type:`ca_t` instances
    referring to this context have been cleared.

Memory management for numbers
-------------------------------------------------------------------------------

.. function:: void ca_init(ca_t x, ca_ctx_t ctx)

    Initializes the variable *x* for use, associating it with the
    context object *ctx*. The value of *x* is set to the rational number 0.

.. function:: void ca_clear(ca_t x, ca_ctx_t ctx)

    Clears the variable *x*.

.. function:: void ca_swap(ca_t x, ca_t y, ca_ctx_t ctx)

    Efficiently swaps the variables *x* and *y*.

.. function:: ca_ptr ca_vec_init(slong n, ca_ctx_t ctx)

    Returns a pointer to a heap-allocated array of *n* initialized
    :type:`ca_struct` entries.

.. function:: void ca_vec_clear(ca_ptr v, slong n, ca_ctx_t ctx)

    Clears an array of *n* initialized :type:`ca_struct` entries.
    This also frees the array *v*.

Input and output
-------------------------------------------------------------------------------

.. function:: void ca_ctx_print(const ca_ctx_t ctx)

    Prints a description of the context *ctx* to standard output.
    This will give a complete listing of the cached fields in *ctx*.

.. function:: void ca_print(const ca_t x, const ca_ctx_t ctx)

    Prints a description of *x* to standard output.

Assignment and specific values
-------------------------------------------------------------------------------

.. function:: void ca_set(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to a copy of *x*.

.. function:: void ca_zero(ca_t res, ca_ctx_t ctx)
              void ca_one(ca_t res, ca_ctx_t ctx)
              void ca_neg_one(ca_t res, ca_ctx_t ctx)

    Sets *res* to the integer 0, 1 or -1. This creates a canonical representation
    of this number as an element of the trivial field `\mathbb{Q}`.

.. function:: void ca_set_si(ca_t res, slong v, ca_ctx_t ctx)
              void ca_set_ui(ca_t res, ulong v, ca_ctx_t ctx)
              void ca_set_fmpz(ca_t res, const fmpz_t v, ca_ctx_t ctx)
              void ca_set_fmpq(ca_t res, const fmpq_t v, ca_ctx_t ctx)

    Sets *res* to the integer or rational number *v*. This creates a canonical
    representation of this number as an element of the trivial field
    `\mathbb{Q}`.

.. function:: void ca_i(ca_t res, ca_ctx_t ctx)
              void ca_neg_i(ca_t res, ca_ctx_t ctx)

    Sets *res* to the imaginary unit `i = \sqrt{-1}`, or its negation `-i`.
    This creates a canonical representation of `i` as the generator of the
    algebraic number field `\mathbb{Q}(i)`.

.. function:: void ca_pi(ca_t res, ca_ctx_t ctx)

    Sets *res* to the constant `\pi`. This creates an element
    of the transcendental number field `\mathbb{Q}(\pi)`.

.. function:: void ca_pi_i(ca_t res, ca_ctx_t ctx)

    Sets *res* to the constant `\pi i`. This creates an element of the
    composite field `\mathbb{Q}(i,\pi)` rather than representing `\pi i`
    (or even `2 \pi i`, which for some purposes would be more elegant)
    as an atomic quantity.

.. function:: void ca_euler(ca_t res, ca_ctx_t ctx)

    Sets *res* to Euler's constant `\gamma`. This creates an element
    of the (transcendental?) number field `\mathbb{Q}(\gamma)`.

.. function:: void ca_unknown(ca_t res, ca_ctx_t ctx)

    Sets *res* to the meta-value *Unknown*.

.. function:: void ca_undefined(ca_t res, ca_ctx_t ctx)

    Sets *res* to *Undefined*.

.. function:: void ca_uinf(ca_t res, ca_ctx_t ctx)

    Sets *res* to unsigned infinity `{\tilde \infty}`.

.. function:: void ca_pos_inf(ca_t res, ca_ctx_t ctx)
              void ca_neg_inf(ca_t res, ca_ctx_t ctx)
              void ca_pos_i_inf(ca_t res, ca_ctx_t ctx)
              void ca_neg_i_inf(ca_t res, ca_ctx_t ctx)

    Sets *res* to the signed infinity `+\infty`, `-\infty`, `+i \infty` or `-i \infty`.

Conversion of algebraic numbers
-------------------------------------------------------------------------------

.. function:: void ca_set_qqbar(ca_t res, const qqbar_t x, ca_ctx_t ctx)

    Sets *res* to the algebraic number *x*.

    If *x* is rational, *res* is set to the canonical representation as
    an element in the trivial field `\mathbb{Q}`.

    If *x* is irrational, this function always sets *res* to an element of
    a univariate number field `\mathbb{Q}(a)`. It will not, for example,
    identify `\sqrt{2} + \sqrt{3}`
    as an element of `\mathbb{Q}(\sqrt{2}, \sqrt{3})`. However, it may
    attempt to find a simpler number field than that generated by *x*
    itself. For example:

    * If *x* is quadratic, it will be expressed as an element of
      `\mathbb{Q}(\sqrt{N})` where *N* has no small repeated factors
      (obtained by performing a smooth factorisation of the discriminant).

    * TODO: if possible, coerce *x* to a low-degree cyclotomic field.

.. function:: int ca_get_fmpz(fmpz_t res, const ca_t x, ca_ctx_t ctx)
              int ca_get_fmpq(fmpz_t res, const ca_t x, ca_ctx_t ctx)
              int ca_get_qqbar(qqbar_t res, const ca_t x, ca_ctx_t ctx)

    Attempts to evaluate *x* to an explicit integer, rational or
    algebraic number. If successful, sets *res* to this number and
    returns 1. If unsuccessful, returns 0.

    The conversion certainly fails if *x* does not represent an integer,
    rational or algebraic number (respectively), but can also fail if *x*
    is too expensive to compute under
    the current evaluation limits.
    In particular, the evaluation will be aborted if an intermediate
    algebraic number (or more precisely, the resultant polynomial prior
    to factorization) exceeds ``CA_OPT_QQBAR_DEG_LIMIT``
    or the coefficients exceed some multiple of ``CA_OPT_PREC_LIMIT``.
    Note that evaluation may hit those limits even if the minimal polynomial
    for *x* itself is small. The conversion can also fail if no algorithm
    has been implemented for the functions appearing in the construction
    of *x*.


Random generation
-------------------------------------------------------------------------------

.. function:: void ca_randtest_rational(ca_t res, flint_rand_t state, slong bits, ca_ctx_t ctx)

    Sets *res* to a random rational number with numerator and denominator
    up to *bits* bits in size.

.. function:: void ca_randtest(ca_t res, flint_rand_t state, slong depth, slong bits, ca_ctx_t ctx)

    Sets *res* to a random number generated by evaluating a random expression.
    The algorithm randomly selects between generating a "simple" number
    (a random rational number or quadratic field element with coefficients
    up to *bits* in size, or a random builtin constant),
    or if *depth* is nonzero, applying a random arithmetic operation or
    function to operands produced through recursive calls with
    *depth* - 1. The output is guaranteed to be a number, not a special value.

.. function:: void ca_randtest_special(ca_t res, flint_rand_t state, slong depth, slong bits, ca_ctx_t ctx)

    Randomly generates either a special value or a number.

Representation properties
-------------------------------------------------------------------------------

The following functions deal with the representation of a :type:`ca_t` and
hence can always be decided quickly and unambiguously. The return value
for predicates is 0 for false and 1 for true.

.. function:: int ca_is_unknown(const ca_t x, ca_ctx_t ctx)

    Returns 1 if *x* is Unknown, and 0 otherwise.

.. function:: int ca_equal_repr(const ca_t x, const ca_t y, ca_ctx_t ctx)

    Returns whether *x* and *y* have identical representation. For field
    elements, this checks if *x* and *y* belong to the same formal field
    (with generators having identical representation) and are represented by
    the same rational function within that field.

    For special values, this tests
    equality of the special values, with *Unknown* handled as if it were
    a value rather than a meta-value: that is, *Unknown* = *Unknown* gives 1,
    and *Unknown* = *y* gives 0 for any other kind of value *y*.
    If neither *x* nor *y* is *Unknown*, then representation equality
    implies that *x* and *y* describe to the same mathematical value, but if
    either operand is *Unknown*, the result is meaningless for
    mathematical comparison.

.. function:: int ca_cmp_repr(const ca_t x, const ca_t y, ca_ctx_t ctx)

    Compares the representations of *x* and *y* in a canonical sort order,
    returning -1, 0 or 1. This only performs a lexicographic comparison
    of the representations of *x* and *y*; the return value does not say
    anything meaningful about the numbers represented by *x* and *y*.


Value predicates
-------------------------------------------------------------------------------

The following predicates check a mathematical property which might
not be effectively decidable. The result is a :type:`truth_t` to allow
representing an unknown outcome.

.. function:: truth_t ca_check_is_number(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is a number. The result is ``T_TRUE`` is *x* is
    a field element (and hence a complex number), ``T_FALSE`` if *x* is
    an infinity or *Undefined*, and ``T_UNKNOWN`` if *x* is *Unknown*.

.. function:: truth_t ca_check_is_zero(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_one(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_neg_one(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_i(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_neg_i(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is equal to the number `0`, `1`, `-1`, `i`, or `-i`.

.. function:: truth_t ca_check_is_algebraic(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_rational(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_integer(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is respectively an algebraic number, a rational number,
    or an integer.

.. function:: truth_t ca_check_is_real(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is a real number. Warning: this returns ``T_FALSE`` if *x* is an
    infinity with real sign.

.. function:: truth_t ca_check_is_imaginary(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is an imaginary number. Warning: this returns ``T_FALSE`` if
    *x* is an infinity with imaginary sign.

.. function:: truth_t ca_check_is_nonreal(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is a nonreal complex number.

.. function:: truth_t ca_check_is_undefined(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is the special value *Undefined*.

.. function:: truth_t ca_check_is_infinity(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is any infinity (unsigned or signed).

.. function:: truth_t ca_check_is_uinf(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is unsigned infinity `{\tilde \infty}`.

.. function:: truth_t ca_check_is_signed_inf(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is any signed infinity.

.. function:: truth_t ca_check_is_pos_inf(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_neg_inf(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_pos_i_inf(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_neg_i_inf(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is equal to the signed infinity `+\infty`, `-\infty`,
    `+i \infty`, `-i \infty`, respectively.

Comparisons
-------------------------------------------------------------------------------

.. function:: truth_t ca_check_equal(const ca_t x, const ca_t y, ca_ctx_t ctx)

    Tests `x = y` as a mathematical equality.
    The result is ``T_UNKNOWN`` if either operand is *Unknown*.
    The result may also be ``T_UNKNOWN`` if *x* and *y* are numerically
    indistinguishable and cannot be proved equal or unequal by
    an exact computation.

.. function:: truth_t ca_check_lt(const ca_t x, const ca_t y, ca_ctx_t ctx)
              truth_t ca_check_le(const ca_t x, const ca_t y, ca_ctx_t ctx)
              truth_t ca_check_gt(const ca_t x, const ca_t y, ca_ctx_t ctx)
              truth_t ca_check_ge(const ca_t x, const ca_t y, ca_ctx_t ctx)

    Compares *x* and *y*, implementing the respective operations
    `x < y`, `x \le y`, `x > y`, `x \ge y`.
    Only real numbers and `-\infty` and `+\infty` are considered comparable.
    The result is ``T_FALSE`` (not ``T_UNKNOWN``) if either operand is not
    comparable (being a nonreal complex number, unsigned infinity, or
    undefined).

Field structure operations
-------------------------------------------------------------------------------

.. function:: void ca_merge_fields(ca_t resx, ca_t resy, const ca_t x, const ca_t y, ca_ctx_t ctx)

    Sets *resx* and *resy* to copies of *x* and *y* coerced to a common field.
    Both *x* and *y* must be field elements (not special values).

    In the present implementation, this simply merges the lists of generators,
    avoiding duplication. In the future, it will be able to eliminate
    generators satisfying algebraic relations.

.. function:: void ca_condense_field(ca_t res, ca_ctx_t ctx)

    Attempts to demote the value of *res* to a trivial subfield of its
    current field by removing unused generators. In particular, this demotes
    any obviously rational value to the trivial field `\mathbb{Q}`.

    This function is applied automatically in most operations
    (arithmetic operations, etc.).

Arithmetic
-------------------------------------------------------------------------------

.. function:: void ca_neg(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the negation of *x*.
    For numbers, this operation amounts to a direct negation within the
    formal field.
    For a signed infinity `c \infty`, negation gives `(-c) \infty`; all other
    special values are unchanged.

.. function:: void ca_add_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
              void ca_add_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx)
              void ca_add_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx)
              void ca_add_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx)
              void ca_add(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx)

    Sets *res* to the sum of *x* and *y*. For special values, the following
    rules apply (`c \infty` denotes a signed infinity, `|c| = 1`):

    * `c \infty + d \infty = c \infty` if `c = d`

    * `c \infty + d \infty = \text{Undefined}` if `c \ne d`

    * `\tilde \infty + c \infty = \tilde \infty + \tilde \infty = \text{Undefined}`

    * `c \infty + z = c \infty` if `z \in \mathbb{C}`

    * `\tilde \infty + z = \tilde \infty` if `z \in \mathbb{C}`

    * `z + \text{Undefined} = \text{Undefined}` for any value *z* (including *Unknown*)

    In any other case involving special values, or if the specific case cannot
    be distinguished, the result is *Unknown*.

.. function:: void ca_sub_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
              void ca_sub_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx)
              void ca_sub_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx)
              void ca_sub_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx)
              void ca_sub(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx)

    Sets *res* to the difference of *x* and *y*.
    This is equivalent to computing `x + (-y)`.

.. function:: void ca_mul_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
              void ca_mul_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx)
              void ca_mul_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx)
              void ca_mul_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx)
              void ca_mul(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx)

    Sets *res* to the product of *x* and *y*. For special values, the following
    rules apply (`c \infty` denotes a signed infinity, `|c| = 1`):

    * `c \infty \cdot d \infty = c d \infty`

    * `c \infty \cdot \tilde \infty = \tilde \infty`

    * `\tilde \infty \cdot \tilde \infty = \tilde \infty`

    * `c \infty \cdot z = \operatorname{sgn}(z) c \infty` if `z \in \mathbb{C} \setminus \{0\}`

    * `c \infty \cdot 0 = \text{Undefined}`

    * `\tilde \infty \cdot 0 = \text{Undefined}`

    * `z \cdot  \text{Undefined} = \text{Undefined}` for any value *z* (including *Unknown*)

    In any other case involving special values, or if the specific case cannot
    be distinguished, the result is *Unknown*.

.. function:: void ca_inv(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the multiplicative inverse of *x*. In a univariate algebraic
    number field, this always produces a rational denominator, but the
    denominator might not be rationalized in a multivariate
    field.
    For special values
    and zero, the following rules apply:

    * `1 / (c \infty) = 1 / \tilde \infty = 0`

    * `1 / 0 = \tilde \infty`

    * `1 / \text{Undefined} = \text{Undefined}`

    * `1 / \text{Unknown} = \text{Unknown}`

    If it cannot be determined whether *x* is zero or nonzero,
    the result is *Unknown*.

.. function:: void ca_fmpq_div(ca_t res, const fmpq_t x, const ca_t y, ca_ctx_t ctx)
              void ca_fmpz_div(ca_t res, const fmpz_t x, const ca_t y, ca_ctx_t ctx)
              void ca_ui_div(ca_t res, ulong x, const ca_t y, ca_ctx_t ctx)
              void ca_si_div(ca_t res, slong x, const ca_t y, ca_ctx_t ctx)
              void ca_div_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
              void ca_div_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx)
              void ca_div_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx)
              void ca_div_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx)
              void ca_div(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx)

    Sets *res* to the quotient of *x* and *y*. This is equivalent
    to computing `x \cdot (1 / y)`. For special values and division
    by zero, this implies the following rules
    (`c \infty` denotes a signed infinity, `|c| = 1`):

    * `(c \infty) / (d \infty) = (c \infty) / \tilde \infty = \tilde \infty / (c \infty) = \tilde \infty / \tilde \infty = \text{Undefined}`

    * `c \infty / z = (c / \operatorname{sgn}(z)) \infty` if `z \in \mathbb{C} \setminus \{0\}`

    * `c \infty / 0 = \tilde \infty / 0 = \tilde \infty`

    * `z / (c \infty) = z / \tilde \infty = 0` if `z \in \mathbb{C}`

    * `z / 0 = \tilde \infty` if `z \in \mathbb{C} \setminus \{0\}`

    * `0 / 0 = \text{Undefined}`

    * `z / \text{Undefined} = \text{Undefined}` for any value *z* (including *Unknown*)

    * `\text{Undefined} / z = \text{Undefined}` for any value *z* (including *Unknown*)

    In any other case involving special values, or if the specific case cannot
    be distinguished, the result is *Unknown*.

Roots
-------------------------------------------------------------------------------

.. function:: void ca_sqrt(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the square root of *x*.

    For special values, the following definitions apply:

    * `\sqrt{c \infty} = \sqrt{c} \infty`

    * `\sqrt{\tilde \infty} = \tilde \infty`.

    * Both *Undefined* and *Unknown* map to themselves.

    This function will attempt to simplify its argument through an exact
    computation. It may in particular attempt to simplify `\sqrt{x}` to
    a single element in `\overline{\mathbb{Q}}`.

    In the generic case, this function outputs an element of the formal
    field `\mathbb{Q}(\sqrt{x})`.

Complex parts
-------------------------------------------------------------------------------

.. function:: void ca_abs(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the absolute value of *x*.

    For special values, the following definitions apply:

    * `|c \infty| = |\tilde \infty| = +\infty`.

    * Both *Undefined* and *Unknown* map to themselves.

    This function will attempt to simplify its argument through an exact
    computation. It may in particular attempt to simplify `|x|` to
    a single element in `\overline{\mathbb{Q}}`.

    In the generic case, this function outputs an element of the formal
    field `\mathbb{Q}(|x|)`.

.. function:: void ca_sgn(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the sign of *x*, defined by

    .. math ::

        \operatorname{sgn}(x) = \begin{cases} 0 & x = 0 \\ \frac{x}{|x|} & x \ne 0 \end{cases}

    for numbers. For special values, the following definitions apply:

    * `\operatorname{sgn}(c \infty) = c`.

    * `\operatorname{sgn}(\tilde \infty) = \operatorname{Undefined}`.

    * Both *Undefined* and *Unknown* map to themselves.

    This function will attempt to simplify its argument through an exact
    computation. It may in particular attempt to simplify `\operatorname{sgn}(x)` to
    a single element in `\overline{\mathbb{Q}}`.

    In the generic case, this function outputs an element of the formal
    field `\mathbb{Q}(\operatorname{sgn}(x))`.

Elementary functions
-------------------------------------------------------------------------------

.. function:: void ca_exp(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the exponential function of *x*.

    For special values, the following definitions apply:

    * `e^{+\infty} = +\infty`

    * `e^{c \infty} = \tilde \infty` if `0 < \operatorname{Re}(c) < 1`.

    * `e^{c \infty} = 0` if `\operatorname{Re}(c) < 0`.

    * `e^{c \infty} = \text{Undefined}` if `\operatorname{Re}(c) = 0`.

    * `e^{\tilde \infty} = \text{Undefined}`.

    * Both *Undefined* and *Unknown* map to themselves.

    In the generic case, this function outputs an element of the formal
    field `\mathbb{Q}(e^x)`.
    Presently, no simplifications are performed apart from `e^0 = 1`.

.. function:: void ca_log(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the natural logarithm of *x*.

    For special values and at the origin, the following definitions apply:

    * For any infinity, `\log(c\infty) = \log(\tilde \infty) = +\infty`.

    * `\log(0) = -\infty`. The result is *Unknown* if deciding `x = 0` fails.

    * Both *Undefined* and *Unknown* map to themselves.

    In the generic case, this function outputs an element of the formal
    field `\mathbb{Q}(\log(x))`.
    Presently, no simplifications are performed apart from `\log(1) = 0`.

Numerical evaluation
-------------------------------------------------------------------------------

.. function:: void ca_get_acb_raw(acb_t res, const ca_t x, slong prec, ca_ctx_t ctx)

    Sets *res* to an enclosure of the numerical value of *x*.
    A working precision of *prec* bits is used internally for the evaluation,
    without adaptive refinement.
    If *x* is any special value, *res* is set to *acb_indeterminate*.


Context options
-------------------------------------------------------------------------------

The *options* member of a :type:`ca_ctx_t` object is an array of *slong*
values controlling simplification behavior and various other settings.
The values of the array at the following indices can be changed by the user
(example: ``ctx->options[CA_OPT_PREC_LIMIT] = 65536``).

.. macro:: CA_OPT_VERBOSE

    Whether to print debug information. Default value: 0.

.. macro:: CA_OPT_PREC_LIMIT

    Maximum precision to use internally for numerical evaluation with Arb.
    This parameter affects the possibility to prove inequalities
    and find simplifications between related extension elements.
    This is not a strict limit; some calculations may use higher precision
    when there is a good reason to do so.
    Default value: 4096.

.. macro:: CA_OPT_QQBAR_DEG_LIMIT

    Maximum degree of :type:`qqbar_t` elements allowed internally during
    simplification of algebraic numbers. This limit may be exceeded
    when the user provides explicit :type:`qqbar_t` input of higher degree.
    Default value: 120.

.. macro:: CA_OPT_LOW_PREC

    Numerical precision to use for fast checks (typically, before attempting
    more expensive operations). Default value: 64.

Internal representation
-------------------------------------------------------------------------------

.. macro:: CA_FMPQ(x)

.. macro:: CA_FMPQ_NUMREF(x)

.. macro:: CA_FMPQ_DENREF(x)

    Assuming that *x* holds an element of the trivial field `\mathbb{Q}`,
    this macro returns a pointer which can be used as an :type:`fmpq_t`,
    or respectively to the numerator or denominator as an :type:`fmpz_t`.

.. macro:: CA_MPOLY_Q(x)

    Assuming that *x* holds a generic field element as data,
    this macro returns a pointer which can be used as
    an :type:`fmpz_mpoly_q_t`.

.. macro:: CA_NF_ELEM(x)

    Assuming that *x* holds an Antic number field element as data,
    this macro returns a pointer which can be used as
    an :type:`nf_elem_t`.


.. function:: void _ca_make_field_element(ca_t x, slong new_index, ca_ctx_t ctx)

    Changes the internal representation of *x* to that of an element
    of the field with index *new_index* in the context object *ctx*.
    This may destroy the value of *x*.

.. function:: void _ca_make_fmpq(ca_t x, ca_ctx_t ctx)

    Changes the internal representation of *x* to that of an element of
    the trivial field `\mathbb{Q}`. This may destroy the value of *x*.

Extension field objects
-------------------------------------------------------------------------------

The following methods are intended for internal use.
The user should normally only manipulate :type:`ca_t` instances,
leaving the construction of field objects to the context object.

.. type:: ca_field_struct

.. type:: ca_field_t

    Represents a formal field.

.. function:: void ca_field_init_qq(ca_field_t K)

    Initializes *K* to represent the trivial field `\mathbb{Q}`.

.. function:: void ca_field_init_nf(ca_field_t K, const qqbar_t x)

    Initializes *K* to represent the algebraic number field `\mathbb{Q}(x)`.

.. function:: void ca_field_init_const(ca_field_t K, ulong func)

    Initializes *K* to represent the field
    `\mathbb{Q}(x)` where *x* is a builtin constant defined by
    *func* (example: *func* = *CA_Pi* for `x = \pi`).

.. function:: void ca_field_init_fx(ca_field_t K, ulong func, const ca_t x, ca_ctx_t ctx)

    Initializes *K* to represent the field
    `\mathbb{Q}(a)` where `a = f(x)`, given a number *x* and a builtin
    univariate function *func* (example: *func* = *CA_Exp* for `e^x`).

.. function:: void ca_field_init_multi(ca_field_t K, slong len)

    Initializes *K* to represent a multivariate field
    `\mathbb{Q}(a_1, \ldots, a_n)` in *n*
    extension elements. The extension elements must subsequently be
    assigned one by one using :func:`ca_field_set_ext`.

.. function:: void ca_field_set_ext(ca_field_t K, slong i, slong x_index, ca_ctx_t ctx)

    Sets the extension element at position *i* (here indexed from 0) of *K*
    to the generator of the field with index *x_index* in *ctx*.
    (It is assumed that the generating field is a univariate field.)

    This only inserts a shallow reference: the field at index *x_index* must
    be kept alive until *K* has been cleared.

.. function:: void ca_field_clear(ca_field_t K)

    Clears the field *K*.

.. function:: void ca_field_print(const ca_field_t K, const ca_ctx_t ctx)

    Prints a description of the field *K* to standard output.

.. function:: int ca_field_cmp(const ca_field_t K1, const ca_field_t K2, ca_ctx_t ctx)

    Compares the field objects *K1* and *K2* in a canonical sort order,
    returning -1, 0 or 1. This only performs a lexicographic comparison
    of the representations of *K1* and *K2*; the return value does not say
    anything meaningful about the relative structures of *K1* and *K2*
    as mathematical fields.


.. raw:: latex

    \newpage

