dnl Define stuff here
define(`__function_prefix',`.. function::')dnl
dnl
define(`addition',dnl
m4_assert_numargs(1)dnl
`__function_prefix $1_add($1_t r, $1_t a, $1_t b)

poopy doopy'dnl
)dnl
dnl
dnl Change quotes to something that is not used anywhere
changequote({{{,}}})dnl
