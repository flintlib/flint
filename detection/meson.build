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

if have_fft_small
    cfg_data.set('FLINT_HAVE_FFT_SMALL', 1,
        description: 'Define to use the fft_small module')
endif

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
  # FLINT_HAVE_ASSEMBLY_armv8
elif host_machine.cpu_family() == 'x86_64'
  assembly_cpu_supported = true
  assembly_dir = 'mpn_extras/x86_64/broadwell'
  # FLINT_HAVE_ASSEMBLY_x86_64_adx
else
  assembly_cpu_supported = false
endif

assembly_opt = assembly_opt.require(assembly_cpu_supported,
    error_message: 'No assembly support for host CPU')

assembly_opt = assembly_opt.enable_auto_if(true)
have_assembly = assembly_opt.enabled()

message('ASSEMBLY: ', have_assembly ? 'enabled' : 'disabled')
