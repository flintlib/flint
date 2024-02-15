dnl
dnl Copyright (C) 2023 Albin Ahlbäck
dnl
dnl This file is part of FLINT.
dnl
dnl FLINT is free software: you can redistribute it and/or modify it under
dnl the terms of the GNU Lesser General Public License (LGPL) as published
dnl by the Free Software Foundation; either version 3 of the License, or
dnl (at your option) any later version.  See <https://www.gnu.org/licenses/>.
dnl

include(`config.m4')

dnl TODO:
dnl * Store in memory to extend to larger cases.
dnl * Minimize 32-bit instructions.
dnl * Why is not adc-chains used more? Replace more adcx chains.
dnl * Fix latencies.

define(`rp', `%rdi')
define(`ap', `%rsi')

define(`s0', `%rcx')
define(`s1', `%r8')
define(`s2', `%r9')
define(`s3', `%r10')
define(`s4', `%r11')
define(`s5', `%rbx')
define(`s6', `%rbp')
define(`s7', `%r12')
define(`s8', `%r13')
define(`s9', `%r14')
define(`s10', `%r15')

define(`sx', `%rax')

	TEXT

	ALIGN(16)
PROLOGUE(flint_mpn_sqr_1)
	mov	0*8(ap), %rdx
	mulx	%rdx, %rcx, %rax
	mov	%rcx, 0*8(rp)
	mov	%rax, 1*8(rp)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqr_2)
	mov	0*8(ap), %rdx
	xor	%r9d, %r9d
	mulx	1*8(ap), %rcx, %r8
	add	%rcx, %rcx
	adc	%r8, %r8
	adc	%r9, %r9
	mulx	%rdx, %r10, %rax
	mov	%r10, 0*8(rp)
	add	%rax, %rcx
	mov	1*8(ap), %rdx
	mov	%rcx, 1*8(rp)
	mulx	%rdx, %r10, %rax
	adc	%r10, %r8
	adc	%r9, %rax
	mov	%r8, 2*8(rp)
	mov	%rax, 3*8(rp)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqr_3)
	mov	0*8(ap), %rdx
	xor	%r11d, %r11d
	mulx	1*8(ap), %rcx, %rax
	mulx	2*8(ap), %r8, %r9
	mov	1*8(ap), %rdx
	add	%rax, %r8
	mulx	2*8(ap), %rax, %r10
	adc	%rax, %r9
	adc	%r11, %r10
	mov	0*8(ap), %rdx
	add	%rcx, %rcx
	adc	%r8, %r8
	adc	%r9, %r9
	adc	%r10, %r10
	adc	%r11, %r11
	mulx	%rdx, %rdx, %rax
	mov	%rdx, 0*8(rp)
	add	%rax, %rcx
	mov	1*8(ap), %rdx
	mov	%rcx, 1*8(rp)
	mulx	%rdx, %rdx, %rax
	adc	%rdx, %r8
	adc	%rax, %r9
	mov	2*8(ap), %rdx
	mov	%r8, 2*8(rp)
	mov	%r9, 3*8(rp)
	mulx	%rdx, %rdx, %rax
	adc	%rdx, %r10
	adc	%r11, %rax
	mov	%r10, 4*8(rp)
	mov	%rax, 5*8(rp)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqr_4)
	mov	0*8(ap), %rdx
	push	%rbx
	push	%rbp
	xor	%ebp, %ebp
	mulx	1*8(ap), %rcx, %r11
	mulx	2*8(ap), %r8, %rbx
	mulx	3*8(ap), %r9, %r10
	mov	1*8(ap), %rdx
	adox	%r11, %r8
	adox	%rbx, %r9
	mulx	2*8(ap), %r11, %rbx
	adcx	%r11, %r9
	adcx	%rbx, %r10
	mulx	3*8(ap), %rbx, %r11
	mov	2*8(ap), %rdx
	adox	%rbx, %r10
	adox	%rbp, %r11
	mulx	3*8(ap), %rdx, %rbx
	adcx	%rdx, %r11
	adc	%rbp, %rbx
	mov	0*8(ap), %rdx
	add	%rcx, %rcx
	adc	%r8, %r8
	adc	%r9, %r9
	adc	%r10, %r10
	adc	%r11, %r11
	adc	%rbx, %rbx
	setc	%bpl
	mulx	%rdx, %rdx, %rax
	mov	%rdx, 0*8(rp)
	add	%rax, %rcx
	mov	1*8(ap), %rdx
	mov	%rcx, 1*8(rp)
	mulx	%rdx, %rdx, %rax
	adc	%rdx, %r8
	adc	%rax, %r9
	mov	2*8(ap), %rdx
	mov	%r8, 2*8(rp)
	mov	%r9, 3*8(rp)
	mulx	%rdx, %rcx, %rax
	adc	%rcx, %r10
	adc	%rax, %r11
	mov	3*8(ap), %rdx
	mov	%r10, 4*8(rp)
	mov	%r11, 5*8(rp)
	mulx	%rdx, %rcx, %rax
	adc	%rcx, %rbx
	adc	%rbp, %rax
	mov	%rbx, 6*8(rp)
	mov	%rax, 7*8(rp)
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqr_5)
	mov	0*8(ap), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	xor	%r13d, %r13d
	mulx	1*8(ap), %rcx, %rbx
	mulx	2*8(ap), %r8, %rbp
	mulx	3*8(ap), %r9, %r12
	mulx	4*8(ap), %r10, %r11
	adcx	%rbx, %r8
	adcx	%rbp, %r9
	adcx	%r12, %r10
	adcx	%r13, %r11
	mov	1*8(ap), %rdx
	mulx	2*8(ap), %r12, %rbx
	adcx	%r12, %r9
	adcx	%rbx, %r10
	mulx	3*8(ap), %r12, %rbx
	adox	%r12, %r10
	adox	%rbx, %r11
	mulx	4*8(ap), %r12, %rbx
	adcx	%r12, %r11
	adcx	%r13, %rbx
	mov	2*8(ap), %rdx
	mulx	3*8(ap), %r12, %rbp
	adcx	%r12, %r11
	adcx	%rbp, %rbx
	mulx	4*8(ap), %r12, %rbp
	adox	%r12, %rbx
	adox	%r13, %rbp
	mov	3*8(ap), %rdx
	mulx	4*8(ap), %rdx, %r12
	adcx	%rdx, %rbp
	adc	%r13, %r12
	mov	0*8(ap), %rdx
	add	%rcx, %rcx
	adc	%r8, %r8
	adc	%r9, %r9
	adc	%r10, %r10
	adc	%r11, %r11
	adc	%rbx, %rbx
	adc	%rbp, %rbp
	adc	%r12, %r12
	adc	%r13, %r13
	mulx	%rdx, %rdx, %rax
	mov	%rdx, 0*8(rp)
	add	%rax, %rcx
	mov	1*8(ap), %rdx
	mov	%rcx, 1*8(rp)
	mulx	%rdx, %rdx, %rax
	adc	%rdx, %r8
	adc	%rax, %r9
	mov	2*8(ap), %rdx
	mov	%r8, 2*8(rp)
	mov	%r9, 3*8(rp)
	mulx	%rdx, %rcx, %rax
	adc	%rcx, %r10
	adc	%rax, %r11
	mov	3*8(ap), %rdx
	mov	%r10, 4*8(rp)
	mov	%r11, 5*8(rp)
	mulx	%rdx, %rcx, %rax
	adc	%rcx, %rbx
	adc	%rax, %rbp
	mov	4*8(ap), %rdx
	mov	%rbx, 6*8(rp)
	mov	%rbp, 7*8(rp)
	mulx	%rdx, %rcx, %rax
	adc	%rcx, %r12
	adc	%r13, %rax
	mov	%r12, 8*8(rp)
	mov	%rax, 9*8(rp)
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqr_6)
	mov	0*8(ap), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15
	xor	%r15d, %r15d
	mulx	1*8(ap), %rcx, %rbp
	mulx	2*8(ap), %r8, %r12
	mulx	3*8(ap), %r9, %r13
	mulx	4*8(ap), %r10, %r14
	mulx	5*8(ap), %r11, %rbx
	adcx	%rbp, %r8
	adcx	%r12, %r9
	adcx	%r13, %r10
	adcx	%r14, %r11
	adcx	%r15, %rbx
	mov	1*8(ap), %rdx
	mulx	2*8(ap), %r14, %rbp
	adcx	%r14, %r9
	adcx	%rbp, %r10
	mulx	3*8(ap), %r14, %rbp
	adox	%r14, %r10
	adox	%rbp, %r11
	mulx	4*8(ap), %r14, %rbp
	adcx	%r14, %r11
	adcx	%rbp, %rbx
	mulx	5*8(ap), %r14, %rbp
	adox	%r14, %rbx
	adox	%r15, %rbp
	adcx	%r15, %rbp
	mov	2*8(ap), %rdx
	mulx	3*8(ap), %r14, %r12
	adcx	%r14, %r11
	adcx	%r12, %rbx
	mulx	4*8(ap), %r14, %r12
	adox	%r14, %rbx
	adox	%r12, %rbp
	mulx	5*8(ap), %r14, %r12
	adcx	%r14, %rbp
	adcx	%r15, %r12
	mov	3*8(ap), %rdx
	mulx	4*8(ap), %r14, %r13
	adcx	%r14, %rbp
	adcx	%r13, %r12
	mulx	5*8(ap), %r14, %r13
	adox	%r14, %r12
	adox	%r15, %r13
	mov	4*8(ap), %rdx
	mulx	5*8(ap), %rdx, %r14
	adcx	%rdx, %r13
	adc	%r15, %r14
	mov	0*8(ap), %rdx
	add	%rcx, %rcx
	adc	%r8, %r8
	adc	%r9, %r9
	adc	%r10, %r10
	adc	%r11, %r11
	adc	%rbx, %rbx
	adc	%rbp, %rbp
	adc	%r12, %r12
	adc	%r13, %r13
	adc	%r14, %r14
	setc	%r15b
	mulx	%rdx, %rdx, %rax
	mov	%rdx, 0*8(rp)
	mov	1*8(ap), %rdx
	add	%rax, %rcx
	mov	%rcx, 1*8(rp)
	mulx	%rdx, %rdx, %rax
	adc	%rdx, %r8
	adc	%rax, %r9
	mov	2*8(ap), %rdx
	mov	%r8, 2*8(rp)
	mov	%r9, 3*8(rp)
	mulx	%rdx, %rcx, %rax
	adc	%rcx, %r10
	adc	%rax, %r11
	mov	3*8(ap), %rdx
	mov	%r10, 4*8(rp)
	mov	%r11, 5*8(rp)
	mulx	%rdx, %rcx, %rax
	adc	%rcx, %rbx
	adc	%rax, %rbp
	mov	4*8(ap), %rdx
	mov	%rbx, 6*8(rp)
	mov	%rbp, 7*8(rp)
	mulx	%rdx, %rcx, %rax
	adc	%rcx, %r12
	adc	%rax, %r13
	mov	5*8(ap), %rdx
	mov	%r12, 8*8(rp)
	mov	%r13, 9*8(rp)
	mulx	%rdx, %rcx, %rax
	adc	%rcx, %r14
	adc	%r15, %rax
	mov	%r14, 10*8(rp)
	mov	%rax, 11*8(rp)
	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqr_7)
	mov	0*8(ap), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15
	xor	%r10d, %r10d
	mulx	1*8(ap), %rcx, %r11
	mulx	2*8(ap), %r8, %rbx
	mulx	3*8(ap), %r9, %r14
	adcx	%r11, %r8
	adcx	%rbx, %r9
	mulx	4*8(ap), %rax, %r12
	mulx	5*8(ap), %r11, %r13
	mulx	6*8(ap), %rbx, %rbp
	adcx	%r14, %rax
	adcx	%r12, %r11
	mov	1*8(ap), %rdx
	adcx	%r13, %rbx
	adcx	%r10, %rbp
	mulx	2*8(ap), %r14, %r12
	mulx	3*8(ap), %r15, %r13
	adcx	%r14, %r9
	adcx	%r12, %rax
	adox	%r15, %rax
	adox	%r13, %r11
	mulx	4*8(ap), %r14, %r12
	mulx	5*8(ap), %r15, %r13
	adcx	%r14, %r11
	adcx	%r12, %rbx
	mulx	6*8(ap), %r14, %r12
	adox	%r15, %rbx
	adox	%r13, %rbp
	mov	2*8(ap), %rdx
	adcx	%r14, %rbp
	adox	%r10, %r12
	adcx	%r10, %r12
	mulx	3*8(ap), %r15, %r13
	adox	%r15, %r11
	adox	%r13, %rbx
	mulx	4*8(ap), %r15, %r13
	adcx	%r15, %rbx
	adcx	%r13, %rbp
	mulx	5*8(ap), %r15, %r13
	adox	%r15, %rbp
	adox	%r13, %r12
	mulx	6*8(ap), %r15, %r13
	adcx	%r15, %r12
	adox	%r10, %r13
	mov	3*8(ap), %rdx
	adcx	%r10, %r13
	mulx	4*8(ap), %r15, %r14
	adcx	%r15, %rbp
	adcx	%r14, %r12
	mulx	5*8(ap), %r15, %r14
	adox	%r15, %r12
	adox	%r14, %r13
	mulx	6*8(ap), %r15, %r14
	mov	4*8(ap), %rdx
	adcx	%r15, %r13
	adcx	%r10, %r14
	mulx	5*8(ap), %r10, %r15
	adcx	%r10, %r13
	adcx	%r15, %r14
	mulx	6*8(ap), %r10, %r15
	mov	$0, %edx
	adox	%r10, %r14
	adox	%rdx, %r15
	mov	5*8(ap), %rdx
	mulx	6*8(ap), %rdx, %r10
	adcx	%rdx, %r15
	adc	$0, %r10
	test	%al, %al
	mov	0*8(ap), %rdx
	adcx	%rcx, %rcx
	adcx	%r8, %r8
	adcx	%r9, %r9
	adcx	%rax, %rax
	adcx	%r11, %r11
	adcx	%rbx, %rbx
	push	%rax
	adcx	%rbp, %rbp
	adcx	%r12, %r12
	adcx	%r13, %r13
	adcx	%r14, %r14
	adcx	%r15, %r15
	adcx	%r10, %r10
	mulx	%rdx, %rdx, %rax
	mov	%rdx, 0*8(rp)
	adox	%rax, %rcx
	mov	1*8(ap), %rdx
	mov	%rcx, 1*8(rp)
	mulx	%rdx, %rdx, %rax
	adox	%rdx, %r8
	adox	%rax, %r9
	pop	%rcx
	mov	2*8(ap), %rdx
	mov	%r8, 2*8(rp)
	mov	%r9, 3*8(rp)
	mulx	%rdx, %rdx, %rax
	adox	%rdx, %rcx
	adox	%rax, %r11
	mov	3*8(ap), %rdx
	mov	%rcx, 4*8(rp)
	mov	%r11, 5*8(rp)
	mulx	%rdx, %r8, %rax
	adox	%r8, %rbx
	adox	%rax, %rbp
	mov	4*8(ap), %rdx
	mov	%rbx, 6*8(rp)
	mov	%rbp, 7*8(rp)
	mulx	%rdx, %r8, %rax
	adox	%r8, %r12
	adox	%rax, %r13
	mov	5*8(ap), %rdx
	mov	%r12, 8*8(rp)
	mov	%r13, 9*8(rp)
	mulx	%rdx, %r8, %rax
	adox	%r8, %r14
	adox	%rax, %r15
	mov	$0, %ecx
	mov	6*8(ap), %rdx
	mov	%r14, 10*8(rp)
	mov	%r15, 11*8(rp)
	mulx	%rdx, %r8, %rax
	adox	%r8, %r10
	adox	%rcx, %rax
	adcx	%rcx, %rax
	mov	%r10, 12*8(rp)
	mov	%rax, 13*8(rp)
	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()
