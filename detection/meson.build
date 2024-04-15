# -----------------------------------------------------------------------------
#   FLINT configuration
# -----------------------------------------------------------------------------

# Data set on cfg_data is used to generate flint-config.h
# Data set on cfg_data_internal (below) is also used to generate config.h
cfg_data = configuration_data()

# -----------------------------------------------------------------------------
#   Project options
# -----------------------------------------------------------------------------

gmp_internals_opt = get_option('gmp_internals')
pthread_opt = get_option('pthread')

# -----------------------------------------------------------------------------
#   fmpz memory management
# -----------------------------------------------------------------------------

enable_gc = get_option('gc').enabled()
enable_reentrant = get_option('reentrant').enabled()
enable_thread_safe = get_option('thread_safe').enabled()

if enable_gc
  if enable_thread_safe
    error('Garbage collection is not thread safe')
  else
    fmpz_c_in = 'fmpz_gc.c'
  endif
else
  if enable_reentrant
    fmpz_c_in = 'fmpz_reentrant.c'
  else
    fmpz_c_in = 'fmpz_single.c'
  endif
endif

# -----------------------------------------------------------------------------
#   libm
# -----------------------------------------------------------------------------

if not cc.has_function('atan2', dependencies: m_dep)
  error('libm is required')
endif

# -----------------------------------------------------------------------------
#   GMP
# -----------------------------------------------------------------------------

gmp_long_long_limb = cc.compiles('''
#include <gmp.h>
#if !defined(_LONG_LONG_LIMB)
# error mp_limb_t != unsigned long long int
#endif
''')

# We need to check with cc.links rather than cc.has_function etc because we
# are checking if we can link symbols not in the headers.

gmp_test_code = '''
#include <gmp.h>

@0@

int main() {
  @1@
  return 0;
}
'''

gmp_base_tests = {
  'mpz_init': ['', 'mpz_t x; __gmpz_init(x);'],
  'mpn_mul_basecase': [
    'void __gmpn_mul_basecase(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);',
    'mp_limb_t a[2], b[2], p[4]; __gmpn_mul_basecase(p, a, (mp_size_t) 1, b, (mp_size_t) 1);',
  ],
}

foreach name, snippets : gmp_base_tests
  _code = gmp_test_code.format(snippets[0], snippets[1])
  if not cc.links(_code, dependencies: gmp_dep, name: name + ' links')
    message('This code did not link with -lgmp: ', _code)
    error('GMP does not have ' + name)
  endif
endforeach

if gmp_internals_opt.allowed()

  gmp_internals_required = {
    'mpn_gcd_11': ['', 'mp_limb_t res, a=0, b=0; res = __gmpn_gcd_11(a, b);'],
    'mpn_div_q': [
      'void __gmpn_div_q(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t, mp_ptr);',
      'mp_limb_t q[2], n[2], d[2], t[2]; __gmpn_div_q(q, t, (mp_size_t) 1, n, (mp_size_t) 1, d);'
    ],
  }

  foreach name, snippets : gmp_internals_required
    _code = gmp_test_code.format(snippets[0], snippets[1])
    if not cc.links(_code, dependencies: gmp_dep, name: name + ' links')
      message('This code did not link with -lgmp: ', _code)
      # Maybe disable gmp_internals rather than error?
      error('GMP does not have ' + name)
    endif
  endforeach

  gmp_internals_optional = {
    'mpn_add_n_sub_n': [
      'void __gmpn_add_n_sub_n(mp_ptr, mp_ptr, mp_ptr, mp_size_t);',
      'mp_limb_t a[2], b[2], c[2], d[2]; __gmpn_add_n_sub_n(a, b, c, (mp_size_t) 1);',
    ],
    'mpn_add_nc': [
      'mp_limb_t __gmpn_add_nc(mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);',
      'mp_limb_t a[2], b[2]; __gmpn_add_nc(a, b, (mp_size_t) 1, (mp_limb_t) 0);',
    ],
    'mpn_addlsh1_n': [
      'mp_limb_t __gmpn_addlsh1_n(mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);',
      'mp_limb_t a[2], b[2]; __gmpn_addlsh1_n(a, b, (mp_size_t) 1, (mp_limb_t) 0);',
    ],
    'mpn_addlsh1_n_ip1': [
      'mp_limb_t __gmpn_addlsh1_n_ip1(mp_ptr, mp_srcptr, mp_size_t);',
      'mp_limb_t a[2], b[2]; __gmpn_addlsh1_n_ip1(a, b, (mp_size_t) 1);',
    ],
    'mpn_addmul_2': [
      'void __gmpn_addmul_2(mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);',
      'mp_limb_t a[2], b[2]; __gmpn_addmul_2(a, b, (mp_size_t) 1, (mp_limb_t) 0);',
    ],
    'mpn_modexact_1_odd': [
      'mp_limb_t __gmpn_modexact_1_odd(mp_srcptr, mp_size_t, mp_limb_t);',
      'mp_limb_t a[2]; __gmpn_modexact_1_odd(a, (mp_size_t) 1, (mp_limb_t) 0);',
    ],
    'mpn_rsh1add_n': [
      'mp_limb_t __gmpn_rsh1add_n(mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);',
      'mp_limb_t a[2], b[2]; __gmpn_rsh1add_n(a, b, (mp_size_t) 1, (mp_limb_t) 0);',
    ],
    'mpn_rsh1sub_n': [
      'mp_limb_t __gmpn_rsh1sub_n(mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);',
      'mp_limb_t a[2], b[2]; __gmpn_rsh1sub_n(a, b, (mp_size_t) 1, (mp_limb_t) 0);',
    ],
    'mpn_sub_nc': [
      'mp_limb_t __gmpn_sub_nc(mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);',
      'mp_limb_t a[2], b[2]; __gmpn_sub_nc(a, b, (mp_size_t) 1, (mp_limb_t) 0);',
    ],
  }

  foreach name, snippets : gmp_internals_optional
    _code = gmp_test_code.format(snippets[0], snippets[1])
    if cc.links(_code, dependencies: gmp_dep, name: name + ' links')
      cfg_data.set('FLINT_HAVE_NATIVE_' + name, 1,
        description: 'Define if GMP has ' + name)
    endif
  endforeach

endif

# -----------------------------------------------------------------------------
#    fft_small
# -----------------------------------------------------------------------------

fft_small_opt = get_option('fft_small')

# check for headers used by fft_small
# XXX: Need to actually compile a test program to check for AVX support
has_avx = cc.has_header('immintrin.h')
has_neon = cc.has_header('arm_neon.h')

fft_small_supported = has_avx or has_neon

# Error if fft_small enabled but not possible, otherwise disable if auto
fft_small_opt = fft_small_opt.require(fft_small_supported,
    error_message: 'host CPU does not support AVX/NEON. Set fft_small to disabled or auto')

# If fft_small was auto we enable it if possible
fft_small_opt = fft_small_opt.enable_auto_if(fft_small_supported)
have_fft_small = fft_small_opt.enabled()

message('FFT_SMALL: ', have_fft_small ? 'enabled' : 'disabled')

# -----------------------------------------------------------------------------
#   assembly
# -----------------------------------------------------------------------------

assembly_opt = get_option('assembly')

# Might be better to check for known good operating systems rather than
# exclude known bad ones.
assembly_opt = assembly_opt.require(host_machine.system() != 'windows',
    error_message: 'No assembly support for Windows')

# Checking for assembly support needs to be more complicated than this. ARM
# needs to be v8 and not sure exactly what to check for x86_64. For now
# the checks below would attempt to build assembly on some systems that do not
# support it. The feature can be disabled with -Dassembly=disabled to build on
# those systems.
if host_machine.cpu_family() == 'aarch64'
  assembly_cpu_supported = true
  assembly_dir = 'mpn_extras/arm64'
  FLINT_HAVE_ASSEMBLY = 'armv8'
elif host_machine.cpu_family() == 'x86_64'
  assembly_cpu_supported = true
  assembly_dir = 'mpn_extras/x86_64/broadwell'
  FLINT_HAVE_ASSEMBLY = 'x86_64_adx'
else
  assembly_cpu_supported = false
  FLINT_HAVE_ASSEMBLY = ''
endif

assembly_opt = assembly_opt.require(assembly_cpu_supported,
    error_message: 'No assembly support for host CPU')

assembly_opt = assembly_opt.enable_auto_if(true)
have_assembly = assembly_opt.enabled()

message('ASSEMBLY: ', have_assembly ? 'enabled' : 'disabled')

# XXX: Need better detection of ASM language support
if have_assembly
  LSYM_PREFIX = 'L'
endif

# -----------------------------------------------------------------------------
#   mparam
# -----------------------------------------------------------------------------

flint_mparam = 'generic'

if host_machine.system() != 'windows'
  if host_machine.cpu_family() == 'aarch64'
    if host_machine.system() == 'darwin'
      flint_mparam = 'arm64/applem1'
    else
      flint_mparam = 'arm64'
    endif
  elif host_machine.cpu_family() == 'x86_64'
    # XXX: Need better CPU detection
    if have_assembly
      flint_mparam = 'x86_64/broadwell'
    else
      flint_mparam = 'generic'
    endif
  endif
endif

# -----------------------------------------------------------------------------
#  external configuration for flint_config.h
# -----------------------------------------------------------------------------

# cfg_data.set('FLINT_BIG_ENDIAN', 1,
#   description: 'Define if system is big endian.')

if have_assembly and FLINT_HAVE_ASSEMBLY == 'armv8'
  cfg_data.set('FLINT_HAVE_ASSEMBLY_armv8', 1,
    description: 'Define if Arm v8 assembly is available')
elif have_assembly and FLINT_HAVE_ASSEMBLY == 'x86_64_adx'
  cfg_data.set('FLINT_HAVE_ASSEMBLY_x86_64_adx', 1,
    description: 'Define if x86_64 ADX assembly is available')
endif

if have_fft_small
  cfg_data.set('FLINT_HAVE_FFT_SMALL', 1,
    description: 'Define to use the fft_small module')
endif

# cfg_data.set('FLINT_KNOW_STRONG_ORDER', 1,
#   description: 'Define if system is strongly ordered')
# cfg_data.set('FLINT_REENTRANT', 0,
#   description: 'Define to enable reentrant.')
cfg_data.set('FLINT_UNROLL_LOOPS', 1,
  description: 'Define to locally unroll some loops')
# cfg_data.set('FLINT_USES_BLAS', 0,
#   description: 'Define to enable BLAS.')
# cfg_data.set('FLINT_USES_CPUSET', 0,
#   description: 'Define if system has cpu_set_t')
# cfg_data.set('FLINT_USES_GC', 0,
#   description: 'Define to enable the Boehm-Demers-Weise garbage collector.')

if pthread_opt.allowed()
  cfg_data.set('FLINT_USES_PTHREAD', 1,
    description: 'Define to enable the use of pthread.')
endif

if enable_thread_safe
  cfg_data.set('FLINT_USES_TLS', 1,
    description: 'Define to enable thread-local storage.')
endif

# cfg_data.set('FLINT_WANT_WANT_ASSERT', 1,
#   description: 'Define to enable use of asserts.')

if gmp_internals_opt.allowed()
  cfg_data.set('FLINT_WANT_GMP_INTERNALS', 1,
    description: 'Define to enable use of GMP internals.')
endif

cfg_data.set('FLINT_WANT_PRETTY_TESTS', 1,
  description: 'Define to enable pretty printing for tests.')

# -----------------------------------------------------------------------------
#  internal configuration for config.h
# -----------------------------------------------------------------------------

cfg_data_internal = cfg_data

#cfg_data_internal.set('AC_APPLE_UNIVERSAL_BUILD', 1)

# XXX: Change this when coverage options are added
#cfg_data_internal.set('FLINT_COVERAGE', 1)

if cc.has_header_symbol('stdlib.h', 'aligned_alloc')
  cfg_data_internal.set('HAVE_ALIGNED_ALLOC', 1,
    description: 'Define to 1 if you have the \'aligned_alloc\' function.')
endif

ac_headers_check = [
  'alloca.h',
  'arm_neon.h',
  'dlfcn.h',
  'errno.h',
  'fenv.h',
  'float.h',
  'immintrin.h', # XXX: has_header is not equivalent to autotools check...
  'inttypes.h',
  'malloc.h',
  'math.h',
  'pthread_np.h',
  'stdarg.h',
  'stdint.h',
  'stdio.h',
  'stdlib.h',
  'string.h',
  'strings.h',
  'sys/param.h',
  'sys/stat.h',
  'sys/types.h',
  'unistd.h',
  'windows.h',
]

foreach h : ac_headers_check
  if cc.has_header(h)
    cfg_data_internal.set('HAVE_' + h.underscorify().to_upper(), 1,
      description: f'Define to 1 if you have the <@h@> header file.')
  endif
endforeach

# XXX: Do we need to check for _alligned_malloc?

cfg_data_internal.set_quoted('LSYM_PREFIX', 'L',
  description: 'Assembler local label prefix')
cfg_data_internal.set_quoted('LT_OBJDIR', '.libs/',
  description: 'Define to the sub-directory where libtool stores uninstalled libraries.')
cfg_data_internal.set_quoted('PACKAGE_BUGREPORT',
  'https://github.com/flintlib/flint/issues/',
  description: 'Define to the address where bug reports for this package should be sent.')
cfg_data_internal.set_quoted('PACKAGE_NAME', 'FLINT',
  description: 'Define to the full name of this package.')
cfg_data_internal.set_quoted('PACKAGE_STRING', 'FLINT ' + FLINT_VERSION_FULL,
  description: 'Define to the full name and version of this package.')
cfg_data_internal.set_quoted('PACKAGE_TARNAME', 'flint',
  description: 'Define to the one symbol short name of this package.')
cfg_data_internal.set_quoted('PACKAGE_URL', 'https://flintlib.org/',
  description: 'Define to the home page for this package.')
cfg_data_internal.set_quoted('PACKAGE_VERSION', FLINT_VERSION_FULL,
  description: 'Define to the version of this package.')

cfg_data_internal.set('STDC_HEADERS', 1,
  description: '''Define to 1 if all of the C89 standard headers exist (not just the ones
   required in a freestanding environment). This macro is provided for
   backward compatibility; new code need not use it.''')

# XXX: set WORDS_BIGENDIAN or undef inline here?
