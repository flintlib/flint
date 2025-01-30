dnl  Copyright (C) 1999-2006, 2011 Free Software Foundation, Inc.
dnl  Copyright (C) 2025 Albin Ahlbäck
dnl
dnl  This file is part of FLINT.
dnl
dnl  FLINT is free software: you can redistribute it and/or modify it under
dnl  the terms of the GNU Lesser General Public License (LGPL) as published
dnl  by the Free Software Foundation; either version 3 of the License, or
dnl  (at your option) any later version.  See <https://www.gnu.org/licenses/>.
dnl
dnl  --------------------------------------------------------------------------
dnl  Basic error handling things.
dnl
dnl  Usage: m4_dollarhash_1_if_noparen_p
dnl
dnl  Expand to 1 if a call "foo" gives $# set to 1 (as opposed to 0 like GNU
dnl  and SysV m4 give).
define(m4_dollarhash_1_if_noparen_test,`$#')dnl
define(m4_dollarhash_1_if_noparen_p,
eval(m4_dollarhash_1_if_noparen_test==1))dnl
undefine(`m4_dollarhash_1_if_noparen_test')dnl
dnl  Usage: m4wrap_prepend(string)
dnl
dnl  Prepend the given string to what will be expanded under m4wrap at the
dnl  end of input.
dnl
dnl  This macro exists to work around variations in m4wrap() behaviour in
dnl  the various m4s (notes at the start of this file).  Don't use m4wrap()
dnl  directly since it will interfere with this scheme.
define(m4wrap_prepend,dnl
m4_assert_numargs(1)dnl
`define(`m4wrap_string',`$1'defn(`m4wrap_string'))')dnl
define(m4wrap_string,`')dnl
define(m4wrap_works_p,dnl
`ifelse(M4WRAP_SPURIOUS,yes,0,1)')dnl
ifelse(m4wrap_works_p,1,dnl
`m4wrap(`m4wrap_string')')dnl
dnl  Usage: m4_file_and_line
dnl
dnl  Expand to the current file and line number, if the GNU m4 extensions
dnl  __file__ and __line__ are available.
dnl
dnl  In GNU m4 1.4 at the end of input when m4wrap text is expanded,
dnl  __file__ is NONE and __line__ is 0, which is not a helpful thing to
dnl  print.  If m4_file_seen() has been called to note the last file seen,
dnl  then that file at a big line number is used, otherwise "end of input"
dnl  is used (although "end of input" won't parse as an error message).
define(m4_file_and_line,dnl
`ifdef(`__file__',dnl
`ifelse(__file__`'__line__,`NONE0',dnl
`ifdef(`m4_file_seen_last',`m4_file_seen_last: 999999: ',`end of input: ')',dnl
`__file__: __line__: ')')')dnl
dnl  Usage: m4_errprint_commas(arg,...)
dnl
dnl  The same as errprint(), but commas are printed between arguments
dnl  instead of spaces.
define(m4_errprint_commas,dnl
`errprint(`$1')dnl
ifelse(eval($#>1),1,`errprint(`,')m4_errprint_commas(shift($@))')')dnl
dnl  Usage: m4_error(args...)
dnl         m4_warning(args...)
dnl
dnl  Print an error message, using m4_errprint_commas, prefixed with the
dnl  current filename and line number (if available).  m4_error sets up to
dnl  give an error exit at the end of processing, m4_warning just prints.
dnl  These macros are the recommended way to print errors.
dnl
dnl  The arguments here should be quoted in the usual way to prevent them
dnl  being expanded when the macro call is read.  (m4_error takes care not
dnl  to do any further expansion.)
dnl
dnl  For example,
dnl
dnl         m4_error(`some error message
dnl         ')
dnl
dnl  which prints
dnl
dnl         foo.asm:123: some error message
dnl
dnl  or if __file__ and __line__ aren't available
dnl
dnl         some error message
dnl
dnl  The "file:line:" format is a basic style, used by gcc and GNU m4, so
dnl  emacs and other editors will recognise it in their normal error message
dnl  parsing.
define(m4_warning,dnl
`m4_errprint_commas(m4_file_and_line`'$@)')dnl
dnl
define(m4_error,dnl
`define(`m4_error_occurred',1)m4_warning($@)dnl
ifelse(m4wrap_works_p,0,`m4exit(1)')')dnl
dnl
define(`m4_error_occurred',0)dnl
dnl
dnl  This m4wrap_prepend() is first, so it'll be executed last.
m4wrap_prepend(dnl
`ifelse(m4_error_occurred,1,dnl
`m4_error(`Errors occurred during m4 processing
')m4exit(1)')')dnl
dnl
dnl  Usage: m4_assert_numargs(num)
dnl
dnl  Put this unquoted on a line on its own at the start of a macro
dnl  definition to add some code to check that num many arguments get passed
dnl  to the macro.  For example,
dnl
dnl         define(foo,
dnl         m4_assert_numargs(2)
dnl         `something `$1' and `$2' blah blah')
dnl
dnl  Then a call like foo(one,two,three) will provoke an error like
dnl
dnl         file:10: foo expected 2 arguments, got 3 arguments
dnl
dnl  Here are some calls and how many arguments they're interpreted as passing.
dnl
dnl         foo(abc,def)  2
dnl         foo(xyz)      1
dnl         foo()         0
dnl         foo          -1
dnl
dnl  The -1 for no parentheses at all means a macro that's meant to be used
dnl  that way can be checked with m4_assert_numargs(-1).  For example,
dnl
dnl         define(SPECIAL_SUFFIX,
dnl         m4_assert_numargs(-1)
dnl         `ifdef(`FOO',`_foo',`_bar')')
dnl
dnl  But as an alternative see also deflit() below where parenthesized
dnl  expressions following a macro are passed through to the output.
dnl
dnl  Note that in BSD m4 there's no way to differentiate calls "foo" and
dnl  "foo()", so in BSD m4 the distinction between the two isn't enforced.
dnl  (In GNU and SysV m4 it can be checked, and is.)
dnl
dnl
dnl  m4_assert_numargs is able to check its own arguments by calling
dnl  assert_numargs_internal directly.
dnl
dnl  m4_doublequote($`'0) expands to ``$0'', whereas ``$`'0'' would expand
dnl  to `$`'0' and do the wrong thing, and likewise for $1.  The same is
dnl  done in other assert macros.
dnl
dnl  $`#' leaves $# in the new macro being defined, and stops # being
dnl  interpreted as a comment character.
dnl
dnl  `dnl ' means an explicit dnl isn't necessary when m4_assert_numargs is
dnl  used.  The space means that if there is a dnl it'll still work.
dnl
dnl  Usage: m4_doublequote(x) expands to ``x''
define(m4_doublequote,dnl
`m4_assert_numargs_internal(`$0',1,$#,len(`$1'))``$1''')dnl
dnl
define(m4_assert_numargs,dnl
`m4_assert_numargs_internal(`$0',1,$#,len(`$1'))dnl
`m4_assert_numargs_internal'(m4_doublequote($`'0),$1,$`#',`len'(m4_doublequote($`'1)))`dnl '')dnl
dnl
dnl  Called: m4_assert_numargs_internal(`macroname',wantargs,$#,len(`$1'))
define(m4_assert_numargs_internal,dnl
`m4_assert_numargs_internal_check(`$1',`$2',m4_numargs_count(`$3',`$4'))')dnl
dnl
dnl  Called: m4_assert_numargs_internal_check(`macroname',wantargs,gotargs)
dnl
dnl  If m4_dollarhash_1_if_noparen_p (BSD m4) then gotargs can be 0 when it
dnl  should be -1.  If wantargs is -1 but gotargs is 0 and the two can't be
dnl  distinguished then it's allowed to pass.
dnl
define(m4_assert_numargs_internal_check,dnl
`ifelse(eval($2 == $3
             || ($2==-1 && $3==0 && m4_dollarhash_1_if_noparen_p)),0,dnl
`m4_error(`$1 expected 'm4_Narguments(`$2')`, got 'm4_Narguments(`$3')dnl
)')')dnl
dnl
dnl  Called: m4_numargs_count($#,len(`$1'))
dnl  If $#==0 then -1 args, if $#==1 but len(`$1')==0 then 0 args, otherwise
dnl  $# args.
define(m4_numargs_count,dnl
`ifelse($1,0, -1,dnl
`ifelse(eval($1==1 && $2-0==0),1, 0, $1)')')dnl
dnl
dnl  Usage: m4_Narguments(N)
dnl  "$1 argument" or "$1 arguments" with the plural according to $1.
define(m4_Narguments,dnl
`$1 argument`'ifelse(`$1',1,,s)')dnl
