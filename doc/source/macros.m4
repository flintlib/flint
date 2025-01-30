include(`config.m4')dnl
dnl############################################################################
dnl Change quotation marks to avoid conflicts
dnl############################################################################
changequote({{{,}}})dnl
dnl############################################################################
dnl helper stuff
dnl############################################################################
define({{{_neg_}}},{{{m4_assert_numargs(1)dnl
-$1}}})dnl
dnl
define({{{_add_}}},{{{m4_assert_numargs(2)dnl
$1 + $2}}})dnl
dnl
define({{{_sub_}}},{{{m4_assert_numargs(2)dnl
$1 - $2}}})dnl
dnl
define({{{_mul_}}},{{{m4_assert_numargs(2)dnl
$1 \cdot $2}}})dnl
dnl
define({{{_div_}}},{{{m4_assert_numargs(2)dnl
$1 / $2}}})dnl
dnl
define({{{_addmul_}}},{{{m4_assert_numargs(3)dnl
_add_($1, _mul_($2, $3))}}})dnl
dnl
define({{{_submul_}}},{{{m4_assert_numargs(3)dnl
_sub_($1, _mul_($2, $3))}}})dnl
dnl
define({{{_lt_}}},{{{m4_assert_numargs(2)dnl
$1 < $2}}})dnl
dnl
define({{{_gt_}}},{{{m4_assert_numargs(2)dnl
$1 > $2}}})dnl
dnl
define({{{_equal_}}},{{{m4_assert_numargs(2)dnl
$1 = $2}}})dnl
dnl############################################################################
dnl set
dnl############################################################################
define({{{desc_set}}},{{{m4_assert_numargs(2)
    Sets `$1` to `$2`.dnl
}}})dnl
define({{{desc_zero}}},{{{m4_assert_numargs(1)
    Sets `$1` to zero.dnl
}}})dnl
define({{{desc_one}}},{{{m4_assert_numargs(1)
    Sets `$1` to one.dnl
}}})dnl
dnl############################################################################
dnl negation, absolute value etc.
dnl############################################################################
define({{{desc_neg}}},{{{m4_assert_numargs(2)
    Sets `$1` to `_neg_($2)`.dnl
}}})dnl
define({{{desc_abs}}},{{{m4_assert_numargs(2)
    Sets `$1` to the absolute value of `$2`.dnl
}}})dnl
dnl############################################################################
dnl basic arithmetic operations
dnl############################################################################
define({{{desc_add}}},{{{m4_assert_numargs(3)
    Sets `$1` to `_add_($2, $3)`.dnl
}}})dnl
define({{{desc_sub}}},{{{m4_assert_numargs(3)
    Sets `$1` to `_sub_($2, $3)`.dnl
}}})dnl
define({{{desc_mul}}},{{{m4_assert_numargs(3)
    Sets `$1` to `_mul_($2, $3)`.dnl
}}})dnl
define({{{desc_divexact}}},{{{m4_assert_numargs(3)
    Sets `$1` to `_div_($2, $3)` under the assumption that the division is
    exact.  If `$3` is zero, an exception is raised.dnl
}}})dnl
dnl############################################################################
dnl extended basic arithmetic operations
dnl############################################################################
define({{{desc_addmul}}},{{{m4_assert_numargs(3)
    Sets `$1` to `_addmul_($1, $2, $3)`.dnl
}}})dnl
define({{{desc_submul}}},{{{m4_assert_numargs(3)
    Sets `$1` to `_submul_($1, $2, $3)`.dnl
}}})dnl
dnl############################################################################
dnl sqrt
dnl############################################################################
define({{{desc_sqrt_nonordered_ring}}},{{{m4_assert_numargs(2)
    If `$2` is a perfect square, sets `$1` to a square root of `$2`
    and returns nonzero.  Otherwise returns zero.dnl
}}})dnl
dnl############################################################################
dnl comparisons
dnl############################################################################
define({{{desc_cmp}}},{{{m4_assert_numargs(2)
    Returns a negative value if `_lt_($1, $2)`, positive value if
    `_gt_($1, $2)`, otherwise returns zero.dnl
}}})dnl
define({{{desc_equal}}},{{{m4_assert_numargs(2)
    Returns nonzero if `_equal_($1, $2)`, otherwise returns zero.
}}})dnl
define({{{desc_is_zero}}},{{{m4_assert_numargs(1)
    Returns nonzero if `_equal_($1, 0)`, otherwise returns zero.dnl
}}})dnl
define({{{desc_is_one}}},{{{m4_assert_numargs(1)
    Returns nonzero if `_equal_($1, 1)`, otherwise returns zero.dnl
}}})dnl
