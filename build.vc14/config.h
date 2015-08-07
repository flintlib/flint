
/* Templates for defines setup by configure.

Copyright 2000, 2001, 2002 Free Software Foundation, Inc.

This file is part of the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MP Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#if _MSC_VER >= 1600 && !defined( HAVE_STDINT_H ) 
#  define HAVE_STDINT_H 1
#endif

#if _MSC_VER >= 1800 
#  define HAVE_INTTYPES_H 1
#endif

#define HAVE_LITTLE_ENDIAN 1

/* Define if you have the `alarm' function. */
#undef HAVE_ALARM

/* Define if alloca() works (via gmp-impl.h). */
#define HAVE_ALLOCA			1

/* Define if you have <alloca.h> and it should be used (not on Ultrix). */
#undef HAVE_ALLOCA_H

/* Define if the compiler accepts gcc style __attribute__ ((const)) */
#undef HAVE_ATTRIBUTE_CONST

/* Define if the compiler accepts gcc style __attribute__ ((malloc)) */
#undef HAVE_ATTRIBUTE_MALLOC

/* Define if the compiler accepts gcc style __attribute__ ((mode (XX))) */
#undef HAVE_ATTRIBUTE_MODE

/* Define if the compiler accepts gcc style __attribute__ ((noreturn)) */
#undef HAVE_ATTRIBUTE_NORETURN

/* Define if tests/libtests has calling conventions checking for the CPU */
#undef HAVE_CALLING_CONVENTIONS

/* Define if you have the `clock' function. */
#define HAVE_CLOCK			1

/* Define if you have the `clock_gettime' function. */
#undef HAVE_CLOCK_GETTIME

/* Define if you have the `cputime' function. */
#undef HAVE_CPUTIME

/* Define to 1 if you have the declaration of `fgetc', and to 0 if you don't.
   */
#define HAVE_DECL_FGETC		1

/* Define to 1 if you have the declaration of `fscanf', and to 0 if you don't.
   */
#define HAVE_DECL_FSCANF	1

/* Define to 1 if you have the declaration of `optarg', and to 0 if you don't.
   */
#define HAVE_DECL_OPTARG	0

/* Define to 1 if you have the declaration of `ungetc', and to 0 if you don't.
   */
#define HAVE_DECL_UNGETC	1

/* Define to 1 if you have the declaration of `vfprintf', and to 0 if you
   don't. */
#define HAVE_DECL_VFPRINTF	1

/* Define if denormalized floats work. */
#define HAVE_DENORMS		1

/* Define if you have the <dlfcn.h> header file. */
#undef HAVE_DLFCN_H

/* Define one (and only one) of the following for the format of a `double'.
   If your format is not among these choices, or you don't know what it is,
   then leave all of them undefined.
   "IEEE_LITTLE_SWAPPED" means little endian, but with the two 4-byte halves
   swapped, as used by ARM CPUs in little endian mode.  */
#undef HAVE_DOUBLE_IEEE_BIG_ENDIAN
#define HAVE_DOUBLE_IEEE_LITTLE_ENDIAN	1
#undef HAVE_DOUBLE_IEEE_LITTLE_SWAPPED
#undef HAVE_DOUBLE_VAX_D
#undef HAVE_DOUBLE_VAX_G
#undef HAVE_DOUBLE_CRAY_CFP

/* Define if you have the <fcntl.h> header file. */
#define HAVE_FCNTL_H		1

/* Define if you have the <fpu_control.h> header file. */
#undef HAVE_FPU_CONTROL_H

/* Define if you have the `getpagesize' function. */
#undef HAVE_GETPAGESIZE

/* Define if you have the `getrusage' function. */
#undef HAVE_GETRUSAGE

/* Define if you have the `gettimeofday' function. */
#undef HAVE_GETTIMEOFDAY

/* Define if 0/0, 1/0, -1/0 and sqrt(-1) work to generate NaN/infinities. */
#define HAVE_INFS			1

/* Define if the system has the type `intmax_t'. */
#undef HAVE_INTMAX_T

/* Define one (just one) of the following for the endiannes of `mp_limb_t'.
   If the endianness is not a simple big or little, or you don't know what
   it is, then leave both of these undefined. */
#undef HAVE_LIMB_BIG_ENDIAN
#define HAVE_LIMB_LITTLE_ENDIAN		1

#define HAVE_STD__LOCALE 1

/* Define if you have the `localeconv' function. */
#define HAVE_LOCALECONV		1

/* Define if you have the <locale.h> header file. */
#define HAVE_LOCALE_H		1

/* now required by MPFR */
#define HAVE_STRUCT_LCONV_DECIMAL_POINT 1
#define HAVE_STRUCT_LCONV_THOUSANDS_SEP 1

/* Define if the system has the type `long double'. */
#define HAVE_LONG_DOUBLE	1

/* Define if the system has the type `long long'. */
#define HAVE_LONG_LONG		1

/* Define if you have the `lrand48' function. */
#undef HAVE_LRAND48

/* Define if you have the <memory.h> header file. */
#define HAVE_MEMORY_H		1

/* Define if you have the `memset' function. */
#define HAVE_MEMSET			1

/* Define if you have the `mmap' function. */
#undef HAVE_MMAP

/* Define if you have the `mprotect' function. */
#undef HAVE_MPROTECT

/* Define if you have the `obstack_vprintf' function. */
#undef HAVE_OBSTACK_VPRINTF

/* Define if you have the `popen' function. */
#undef HAVE_POPEN

/* Define if you have the `processor_info' function. */
#undef HAVE_PROCESSOR_INFO

/* Define if the system has the type `ptrdiff_t'. */
#define HAVE_PTRDIFF_T		1

/* Define if the system has the type `quad_t'. */
#undef HAVE_QUAD_T

#define HAVE_RAISE			1

/* Define if you have the `read_real_time' function. */
#undef HAVE_READ_REAL_TIME

#define HAVE_SIGNAL         1
#define HAVE_SIGNAL_H       1

/* Define if you have the `sigaction' function. */
#undef HAVE_SIGACTION

/* Define if you have the `sigaltstack' function. */
#undef HAVE_SIGALTSTACK

/* Define if you have the `sigstack' function. */
#undef HAVE_SIGSTACK

/* Tune directory speed_cyclecounter, undef=none, 1=32bits, 2=64bits) */
#define HAVE_SPEED_CYCLECOUNTER	2

/* Define if the system has the type `stack_t'. */
#undef HAVE_STACK_T

/* Define if <stdarg.h> exists and works */
#define HAVE_STDARG			1

/* Define if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H		1

/* Define if you have the `strcasecmp' function. */
#undef HAVE_STRCASECMP

/* Define if you have the `strchr' function. */
#define HAVE_STRCHR			1

/* Define if cpp supports the ANSI # stringizing operator. */
#define HAVE_STRINGIZE		1

/* Define if you have the <strings.h> header file. */
#undef HAVE_STRINGS_H

/* Define if you have the <string.h> header file. */
#define HAVE_STRING_H		1

/* Define if you have the `strnlen' function. */
#define HAVE_STRNLEN        1

/* Define if you have the `strtoul' function. */
#define HAVE_STRTOUL		1

/* Define if you have the `sysconf' function. */
#undef HAVE_SYSCONF

/* Define if you have the `sysctl' function. */
#undef HAVE_SYSCTL

/* Define if you have the `sysctlbyname' function. */
#undef HAVE_SYSCTLBYNAME

/* Define if you have the `syssgi' function. */
#undef HAVE_SYSSGI

/* Define if you have the <sys/mman.h> header file. */
#undef HAVE_SYS_MMAN_H

/* Define if you have the <sys/param.h> header file. */
#undef HAVE_SYS_PARAM_H

/* Define if you have the <sys/processor.h> header file. */
#undef HAVE_SYS_PROCESSOR_H

/* Define if you have the <sys/resource.h> header file. */
#undef HAVE_SYS_RESOURCE_H

/* Define if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H		1

/* Define if you have the <sys/sysctl.h> header file. */
#undef HAVE_SYS_SYSCTL_H

/* Define if you have the <sys/syssgi.h> header file. */
#undef HAVE_SYS_SYSSGI_H

/* Define if you have the <sys/systemcfg.h> header file. */
#undef HAVE_SYS_SYSTEMCFG_H

/* Define if you have the <sys/times.h> header file. */
#undef HAVE_SYS_TIMES_H

/* Define if you have the <sys/time.h> header file. */
#undef HAVE_SYS_TIME_H

/* Define if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H	1

/* Define if you have the `times' function. */
#undef HAVE_TIMES

/* Define if you have the <unistd.h> header file. */
#undef HAVE_UNISTD_H

/* Define if you have vsnprintf and it works properly. */
#undef HAVE_VSNPRINTF

/* Assembler local label prefix */
#undef LSYM_PREFIX

/* Define if you have the `fesetround' function via the <fenv.h> header file.
   */
#undef MPFR_HAVE_FESETROUND

#define HAVE_SSTREAM 1

/* Define if compiler has function prototypes */
#define PROTOTYPES			1

/* Define as the return type of signal handlers (`int' or `void'). */
#define RETSIGTYPE	void

/* The size of a `unsigned long', as computed by sizeof. */
#define SIZEOF_UNSIGNED_LONG	4

/* Define if sscanf requires writable inputs */
#undef SSCANF_WRITABLE_INPUT

/* Define if you have the ANSI C header files. */
#define STDC_HEADERS		1

/* ./configure --enable-assert option, to enable some ASSERT()s */
#undef WANT_ASSERT

/* Define if your processor stores words with the most significant byte first
   (like Motorola and SPARC, unlike Intel and VAX). */
#undef WORDS_BIGENDIAN

#define HAVE_PTHREAD 1

/* Define as `__inline' if that's what the C compiler calls it, or to nothing
   if it is not supported. */
#ifndef __cplusplus
#define inline	__inline
#endif

#ifdef HAVE_STDINT_H
#define HAVE_INTMAX_T        1
#define HAVE_UINTMAX_T       1
#define HAVE_PTRDIFF_T       1
#define HAVE_UINT_LEAST32_T  1
#define SIZEOF_UINTMAX_T	 8
#endif
#define NPRINTF_J            1
#define NPRINTF_T            1

#ifdef _MSC_VER
#define access _access
#define strcasecmp _stricmp
#define strncasecmp	_strnicmp
#define alloca _alloca
#define HAVE_STRCASECMP		1
#define HAVE_STRNCASECMP	1
#define MSC_C_(x) #x  
#define MSC_CC_(x)  MSC_C_(x)
#define MSC_VERSION "Microsoft C++ (Version " MSC_CC_(_MSC_FULL_VER) ")"

#if defined (MSC_BUILD_DLL)
#define FLINT_DLL __declspec(dllexport)
#elif defined(MSC_USE_DLL)
#define FLINT_DLL __declspec(dllimport)
#else
#define FLINT_DLL
#endif
#endif
