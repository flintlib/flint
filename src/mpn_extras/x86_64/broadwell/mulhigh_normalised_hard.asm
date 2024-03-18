dnl
dnl Copyright (C) 2024 Albin Ahlbäck
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
dnl * Do stuff in between to avoid latency penalties.
dnl * Convert %rX registers using 32-bit operations to %rax, ..., %rbp
dnl   registers to save a few bytes.

	TEXT

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_normalised_1)
	mov	0*8(%rdx), %rdx
	mulx	0*8(%rsi), %rax, %rcx
	mov	$0, %rdx
	test	%rcx, %rcx
	setns	%dl
	js	L(1)
	add	%rax, %rax
	adc	%rcx, %rcx
L(1):
	mov	%rcx, 0*8(%rdi)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_normalised_2)
	mov	1*8(%rdx), %r11
	mov	0*8(%rdx), %rdx
	xor	%r10d, %r10d
	mulx	0*8(%rsi), %r9, %rax
	mulx	1*8(%rsi), %r9, %rcx
	adcx	%r9, %rax
	adcx	%r10, %rcx
	mov	%r11, %rdx
	mulx	0*8(%rsi), %r9, %r8
	adcx	%r9, %rax
	adcx	%r8, %rcx
	mulx	1*8(%rsi), %r9, %r8
	adox	%r9, %rcx
	adox	%r10, %r8
	adcx	%r10, %r8
	mov	$0, %rdx
	test	%r8, %r8
	setns	%dl
	js	L(2)
	add	%rax, %rax
	adc	%rcx, %rcx
	adc	%r8, %r8
L(2):
	mov	%rcx, 0*8(%rdi)
	mov	%r8, 1*8(%rdi)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_normalised_3)
	push	%rbx

	mov	%rdx, %rcx
	mov	0*8(%rdx), %rdx
	xor	%r9d, %r9d

	mulx	1*8(%rsi), %r8, %rax
	mulx	2*8(%rsi), %r8, %r10
	adcx	%r8, %rax
	adcx	%r9, %r10

	mov	1*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %rbx
	mulx	1*8(%rsi), %r8, %r11
	adcx	%rbx, %rax
	adox	%r8, %rax
	adcx	%r11, %r10
	mulx	2*8(%rsi), %r8, %r11
	adox	%r8, %r10
	adcx	%r9, %r11
	adox	%r9, %r11

	mov	2*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %rbx
	adcx	%r8, %rax
	adcx	%rbx, %r10
	mulx	1*8(%rsi), %r8, %rbx
	adox	%r8, %r10
	adox	%rbx, %r11
	mulx	2*8(%rsi), %r8, %rbx
	adcx	%r8, %r11
	adcx	%r9, %rbx
	adox	%r9, %rbx

	mov	$0, %rdx
	test	%rbx, %rbx
	setns	%dl
	js	L(3)
	add	%rax, %rax
	adc	%r10, %r10
	adc	%r11, %r11
	adc	%rbx, %rbx
L(3):
	mov	%r10, 0*8(%rdi)
	mov	%r11, 1*8(%rdi)
	mov	%rbx, 2*8(%rdi)

	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_normalised_4)
	push	%rbx
	push	%rbp

	mov	%rdx, %rcx
	mov	0*8(%rdx), %rdx
	xor	%r9d, %r9d

	mulx	2*8(%rsi), %r8, %rax
	mulx	3*8(%rsi), %r8, %r10
	adcx	%r8, %rax
	adcx	%r9, %r10

	mov	1*8(%rcx), %rdx
	mulx	1*8(%rsi), %r8, %rbx
	mulx	2*8(%rsi), %r8, %r11
	adcx	%rbx, %rax
	adox	%r8, %rax
	adcx	%r11, %r10
	mulx	3*8(%rsi), %r8, %r11
	adox	%r8, %r10
	adcx	%r9, %r11
	adox	%r9, %r11

	mov	2*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %rbp
	mulx	1*8(%rsi), %r8, %rbx
	adcx	%rbp, %rax
	adox	%r8, %rax
	adcx	%rbx, %r10
	mulx	2*8(%rsi), %r8, %rbx
	adox	%r8, %r10
	adcx	%rbx, %r11
	mulx	3*8(%rsi), %r8, %rbx
	adox	%r8, %r11
	adcx	%r9, %rbx
	adox	%r9, %rbx

	mov	3*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %rbp
	adcx	%r8, %rax
	adcx	%rbp, %r10
	mulx	1*8(%rsi), %r8, %rbp
	adox	%r8, %r10
	adox	%rbp, %r11
	mulx	2*8(%rsi), %r8, %rbp
	adcx	%r8, %r11
	adcx	%rbp, %rbx
	mulx	3*8(%rsi), %r8, %rbp
	adox	%r8, %rbx
	adcx	%r9, %rbp
	adox	%r9, %rbp

	mov	$0, %rdx
	test	%rbp, %rbp
	setns	%dl
	js	L(4)
	add	%rax, %rax
	adc	%r10, %r10
	adc	%r11, %r11
	adc	%rbx, %rbx
	adc	%rbp, %rbp
L(4):
	mov	%r10, 0*8(%rdi)
	mov	%r11, 1*8(%rdi)
	mov	%rbx, 2*8(%rdi)
	mov	%rbp, 3*8(%rdi)

	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_normalised_5)
	push	%rbx
	push	%rbp
	push	%r12

	mov	%rdx, %rcx
	mov	0*8(%rdx), %rdx
	xor	%r9d, %r9d

	mulx	3*8(%rsi), %r8, %rax
	mulx	4*8(%rsi), %r8, %r10
	adcx	%r8, %rax
	adcx	%r9, %r10

	mov	1*8(%rcx), %rdx
	mulx	2*8(%rsi), %r8, %rbx
	mulx	3*8(%rsi), %r8, %r11
	adcx	%rbx, %rax
	adox	%r8, %rax
	adcx	%r11, %r10
	mulx	4*8(%rsi), %r8, %r11
	adox	%r8, %r10
	adcx	%r9, %r11
	adox	%r9, %r11

	mov	2*8(%rcx), %rdx
	mulx	1*8(%rsi), %r8, %rbp
	mulx	2*8(%rsi), %r8, %rbx
	adcx	%rbp, %rax
	adox	%r8, %rax
	adcx	%rbx, %r10
	mulx	3*8(%rsi), %r8, %rbx
	adox	%r8, %r10
	adcx	%rbx, %r11
	mulx	4*8(%rsi), %r8, %rbx
	adox	%r8, %r11
	adcx	%r9, %rbx
	adox	%r9, %rbx

	mov	3*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r12
	mulx	1*8(%rsi), %r8, %rbp
	adcx	%r12, %rax
	adox	%r8, %rax
	adcx	%rbp, %r10
	mulx	2*8(%rsi), %r8, %rbp
	adox	%r8, %r10
	adcx	%rbp, %r11
	mulx	3*8(%rsi), %r8, %rbp
	adox	%r8, %r11
	adcx	%rbp, %rbx
	mulx	4*8(%rsi), %r8, %rbp
	adox	%r8, %rbx
	adcx	%r9, %rbp
	adox	%r9, %rbp

	mov	4*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r12
	adcx	%r8, %rax
	adcx	%r12, %r10
	mulx	1*8(%rsi), %r8, %r12
	adox	%r8, %r10
	adox	%r12, %r11
	mulx	2*8(%rsi), %r8, %r12
	adcx	%r8, %r11
	adcx	%r12, %rbx
	mulx	3*8(%rsi), %r8, %r12
	adox	%r8, %rbx
	adox	%r12, %rbp
	mulx	4*8(%rsi), %r8, %r12
	adcx	%r8, %rbp
	adcx	%r9, %r12
	adox	%r9, %r12

	mov	$0, %rdx
	test	%r12, %r12
	setns	%dl
	js	L(5)
	add	%rax, %rax
	adc	%r10, %r10
	adc	%r11, %r11
	adc	%rbx, %rbx
	adc	%rbp, %rbp
	adc	%r12, %r12
L(5):
	mov	%r10, 0*8(%rdi)
	mov	%r11, 1*8(%rdi)
	mov	%rbx, 2*8(%rdi)
	mov	%rbp, 3*8(%rdi)
	mov	%r12, 4*8(%rdi)

	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_normalised_6)
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	mov	%rdx, %rcx
	mov	0*8(%rdx), %rdx
	xor	%r9d, %r9d

	mulx	4*8(%rsi), %r8, %rax
	mulx	5*8(%rsi), %r8, %r10
	adcx	%r8, %rax
	adcx	%r9, %r10

	mov	1*8(%rcx), %rdx
	mulx	3*8(%rsi), %r8, %rbx
	mulx	4*8(%rsi), %r8, %r11
	adcx	%rbx, %rax
	adox	%r8, %rax
	adcx	%r11, %r10
	mulx	5*8(%rsi), %r8, %r11
	adox	%r8, %r10
	adcx	%r9, %r11
	adox	%r9, %r11

	mov	2*8(%rcx), %rdx
	mulx	2*8(%rsi), %r8, %rbp
	mulx	3*8(%rsi), %r8, %rbx
	adcx	%rbp, %rax
	adox	%r8, %rax
	adcx	%rbx, %r10
	mulx	4*8(%rsi), %r8, %rbx
	adox	%r8, %r10
	adcx	%rbx, %r11
	mulx	5*8(%rsi), %r8, %rbx
	adox	%r8, %r11
	adcx	%r9, %rbx
	adox	%r9, %rbx

	mov	3*8(%rcx), %rdx
	mulx	1*8(%rsi), %r8, %r12
	mulx	2*8(%rsi), %r8, %rbp
	adcx	%r12, %rax
	adox	%r8, %rax
	adcx	%rbp, %r10
	mulx	3*8(%rsi), %r8, %rbp
	adox	%r8, %r10
	adcx	%rbp, %r11
	mulx	4*8(%rsi), %r8, %rbp
	adox	%r8, %r11
	adcx	%rbp, %rbx
	mulx	5*8(%rsi), %r8, %rbp
	adox	%r8, %rbx
	adcx	%r9, %rbp
	adox	%r9, %rbp

	mov	4*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r13
	mulx	1*8(%rsi), %r8, %r12
	adcx	%r13, %rax
	adox	%r8, %rax
	adcx	%r12, %r10
	mulx	2*8(%rsi), %r8, %r12
	adox	%r8, %r10
	adcx	%r12, %r11
	mulx	3*8(%rsi), %r8, %r12
	adox	%r8, %r11
	adcx	%r12, %rbx
	mulx	4*8(%rsi), %r8, %r12
	adox	%r8, %rbx
	adcx	%r12, %rbp
	mulx	5*8(%rsi), %r8, %r12
	adox	%r8, %rbp
	adcx	%r9, %r12
	adox	%r9, %r12

	mov	5*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r13
	adcx	%r8, %rax
	adcx	%r13, %r10
	mulx	1*8(%rsi), %r8, %r13
	adox	%r8, %r10
	adox	%r13, %r11
	mulx	2*8(%rsi), %r8, %r13
	adcx	%r8, %r11
	adcx	%r13, %rbx
	mulx	3*8(%rsi), %r8, %r13
	adox	%r8, %rbx
	adox	%r13, %rbp
	mulx	4*8(%rsi), %r8, %r13
	adcx	%r8, %rbp
	adcx	%r13, %r12
	mulx	5*8(%rsi), %r8, %r13
	adox	%r8, %r12
	adcx	%r9, %r13
	adox	%r9, %r13

	mov	$0, %rdx
	test	%r13, %r13
	setns	%dl
	js	L(6)
	add	%rax, %rax
	adc	%r10, %r10
	adc	%r11, %r11
	adc	%rbx, %rbx
	adc	%rbp, %rbp
	adc	%r12, %r12
	adc	%r13, %r13
L(6):
	mov	%r10, 0*8(%rdi)
	mov	%r11, 1*8(%rdi)
	mov	%rbx, 2*8(%rdi)
	mov	%rbp, 3*8(%rdi)
	mov	%r12, 4*8(%rdi)
	mov	%r13, 5*8(%rdi)

	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_normalised_7)
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	mov	%rdx, %rcx
	mov	0*8(%rdx), %rdx
	xor	%r9d, %r9d

	mulx	5*8(%rsi), %r8, %rax
	mulx	6*8(%rsi), %r8, %r10
	adcx	%r8, %rax
	adcx	%r9, %r10

	mov	1*8(%rcx), %rdx
	mulx	4*8(%rsi), %r8, %rbx
	mulx	5*8(%rsi), %r8, %r11
	adcx	%rbx, %rax
	adox	%r8, %rax
	adcx	%r11, %r10
	mulx	6*8(%rsi), %r8, %r11
	adox	%r8, %r10
	adcx	%r9, %r11
	adox	%r9, %r11

	mov	2*8(%rcx), %rdx
	mulx	3*8(%rsi), %r8, %rbp
	mulx	4*8(%rsi), %r8, %rbx
	adcx	%rbp, %rax
	adox	%r8, %rax
	adcx	%rbx, %r10
	mulx	5*8(%rsi), %r8, %rbx
	adox	%r8, %r10
	adcx	%rbx, %r11
	mulx	6*8(%rsi), %r8, %rbx
	adox	%r8, %r11
	adcx	%r9, %rbx
	adox	%r9, %rbx

	mov	3*8(%rcx), %rdx
	mulx	2*8(%rsi), %r8, %r12
	mulx	3*8(%rsi), %r8, %rbp
	adcx	%r12, %rax
	adox	%r8, %rax
	adcx	%rbp, %r10
	mulx	4*8(%rsi), %r8, %rbp
	adox	%r8, %r10
	adcx	%rbp, %r11
	mulx	5*8(%rsi), %r8, %rbp
	adox	%r8, %r11
	adcx	%rbp, %rbx
	mulx	6*8(%rsi), %r8, %rbp
	adox	%r8, %rbx
	adcx	%r9, %rbp
	adox	%r9, %rbp

	mov	4*8(%rcx), %rdx
	mulx	1*8(%rsi), %r8, %r13
	mulx	2*8(%rsi), %r8, %r12
	adcx	%r13, %rax
	adox	%r8, %rax
	adcx	%r12, %r10
	mulx	3*8(%rsi), %r8, %r12
	adox	%r8, %r10
	adcx	%r12, %r11
	mulx	4*8(%rsi), %r8, %r12
	adox	%r8, %r11
	adcx	%r12, %rbx
	mulx	5*8(%rsi), %r8, %r12
	adox	%r8, %rbx
	adcx	%r12, %rbp
	mulx	6*8(%rsi), %r8, %r12
	adox	%r8, %rbp
	adcx	%r9, %r12
	adox	%r9, %r12

	mov	5*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r14
	mulx	1*8(%rsi), %r8, %r13
	adcx	%r14, %rax
	adox	%r8, %rax
	adcx	%r13, %r10
	mulx	2*8(%rsi), %r8, %r13
	adox	%r8, %r10
	adcx	%r13, %r11
	mulx	3*8(%rsi), %r8, %r13
	adox	%r8, %r11
	adcx	%r13, %rbx
	mulx	4*8(%rsi), %r8, %r13
	adox	%r8, %rbx
	adcx	%r13, %rbp
	mulx	5*8(%rsi), %r8, %r13
	adox	%r8, %rbp
	adcx	%r13, %r12
	mulx	6*8(%rsi), %r8, %r13
	adox	%r8, %r12
	adcx	%r9, %r13
	adox	%r9, %r13

	mov	6*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r14
	adcx	%r8, %rax
	adcx	%r14, %r10
	mulx	1*8(%rsi), %r8, %r14
	adox	%r8, %r10
	adox	%r14, %r11
	mulx	2*8(%rsi), %r8, %r14
	adcx	%r8, %r11
	adcx	%r14, %rbx
	mulx	3*8(%rsi), %r8, %r14
	adox	%r8, %rbx
	adox	%r14, %rbp
	mulx	4*8(%rsi), %r8, %r14
	adcx	%r8, %rbp
	adcx	%r14, %r12
	mulx	5*8(%rsi), %r8, %r14
	adox	%r8, %r12
	adox	%r14, %r13
	mulx	6*8(%rsi), %r8, %r14
	adcx	%r8, %r13
	adcx	%r9, %r14
	adox	%r9, %r14

	mov	$0, %rdx
	test	%r14, %r14
	setns	%dl
	js	L(7)
	add	%rax, %rax
	adc	%r10, %r10
	adc	%r11, %r11
	adc	%rbx, %rbx
	adc	%rbp, %rbp
	adc	%r12, %r12
	adc	%r13, %r13
	adc	%r14, %r14
L(7):
	mov	%r10, 0*8(%rdi)
	mov	%r11, 1*8(%rdi)
	mov	%rbx, 2*8(%rdi)
	mov	%rbp, 3*8(%rdi)
	mov	%r12, 4*8(%rdi)
	mov	%r13, 5*8(%rdi)
	mov	%r14, 6*8(%rdi)

	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_normalised_8)
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15

	mov	%rdx, %rcx
	mov	0*8(%rdx), %rdx
	xor	%r9d, %r9d

	mulx	6*8(%rsi), %r8, %rax
	mulx	7*8(%rsi), %r8, %r10
	adcx	%r8, %rax
	adcx	%r9, %r10

	mov	1*8(%rcx), %rdx
	mulx	5*8(%rsi), %r8, %rbx
	mulx	6*8(%rsi), %r8, %r11
	adcx	%rbx, %rax
	adox	%r8, %rax
	adcx	%r11, %r10
	mulx	7*8(%rsi), %r8, %r11
	adox	%r8, %r10
	adcx	%r9, %r11
	adox	%r9, %r11

	mov	2*8(%rcx), %rdx
	mulx	4*8(%rsi), %r8, %rbp
	mulx	5*8(%rsi), %r8, %rbx
	adcx	%rbp, %rax
	adox	%r8, %rax
	adcx	%rbx, %r10
	mulx	6*8(%rsi), %r8, %rbx
	adox	%r8, %r10
	adcx	%rbx, %r11
	mulx	7*8(%rsi), %r8, %rbx
	adox	%r8, %r11
	adcx	%r9, %rbx
	adox	%r9, %rbx

	mov	3*8(%rcx), %rdx
	mulx	3*8(%rsi), %r8, %r12
	mulx	4*8(%rsi), %r8, %rbp
	adcx	%r12, %rax
	adox	%r8, %rax
	adcx	%rbp, %r10
	mulx	5*8(%rsi), %r8, %rbp
	adox	%r8, %r10
	adcx	%rbp, %r11
	mulx	6*8(%rsi), %r8, %rbp
	adox	%r8, %r11
	adcx	%rbp, %rbx
	mulx	7*8(%rsi), %r8, %rbp
	adox	%r8, %rbx
	adcx	%r9, %rbp
	adox	%r9, %rbp

	mov	4*8(%rcx), %rdx
	mulx	2*8(%rsi), %r8, %r13
	mulx	3*8(%rsi), %r8, %r12
	adcx	%r13, %rax
	adox	%r8, %rax
	adcx	%r12, %r10
	mulx	4*8(%rsi), %r8, %r12
	adox	%r8, %r10
	adcx	%r12, %r11
	mulx	5*8(%rsi), %r8, %r12
	adox	%r8, %r11
	adcx	%r12, %rbx
	mulx	6*8(%rsi), %r8, %r12
	adox	%r8, %rbx
	adcx	%r12, %rbp
	mulx	7*8(%rsi), %r8, %r12
	adox	%r8, %rbp
	adcx	%r9, %r12
	adox	%r9, %r12

	mov	5*8(%rcx), %rdx
	mulx	1*8(%rsi), %r8, %r14
	mulx	2*8(%rsi), %r8, %r13
	adcx	%r14, %rax
	adox	%r8, %rax
	adcx	%r13, %r10
	mulx	3*8(%rsi), %r8, %r13
	adox	%r8, %r10
	adcx	%r13, %r11
	mulx	4*8(%rsi), %r8, %r13
	adox	%r8, %r11
	adcx	%r13, %rbx
	mulx	5*8(%rsi), %r8, %r13
	adox	%r8, %rbx
	adcx	%r13, %rbp
	mulx	6*8(%rsi), %r8, %r13
	adox	%r8, %rbp
	adcx	%r13, %r12
	mulx	7*8(%rsi), %r8, %r13
	adox	%r8, %r12
	adcx	%r9, %r13
	adox	%r9, %r13

	mov	6*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r15
	mulx	1*8(%rsi), %r8, %r14
	adcx	%r15, %rax
	adox	%r8, %rax
	adcx	%r14, %r10
	mulx	2*8(%rsi), %r8, %r14
	adox	%r8, %r10
	adcx	%r14, %r11
	mulx	3*8(%rsi), %r8, %r14
	adox	%r8, %r11
	adcx	%r14, %rbx
	mulx	4*8(%rsi), %r8, %r14
	adox	%r8, %rbx
	adcx	%r14, %rbp
	mulx	5*8(%rsi), %r8, %r14
	adox	%r8, %rbp
	adcx	%r14, %r12
	mulx	6*8(%rsi), %r8, %r14
	adox	%r8, %r12
	adcx	%r14, %r13
	mulx	7*8(%rsi), %r8, %r14
	adox	%r8, %r13
	adcx	%r9, %r14
	adox	%r9, %r14

	mov	7*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r15
	adcx	%r8, %rax
	adcx	%r15, %r10
	mulx	1*8(%rsi), %r8, %r15
	adox	%r8, %r10
	adox	%r15, %r11
	mulx	2*8(%rsi), %r8, %r15
	adcx	%r8, %r11
	adcx	%r15, %rbx
	mulx	3*8(%rsi), %r8, %r15
	adox	%r8, %rbx
	adox	%r15, %rbp
	mulx	4*8(%rsi), %r8, %r15
	adcx	%r8, %rbp
	adcx	%r15, %r12
	mulx	5*8(%rsi), %r8, %r15
	adox	%r8, %r12
	adox	%r15, %r13
	mulx	6*8(%rsi), %r8, %r15
	adcx	%r8, %r13
	adcx	%r15, %r14
	mulx	7*8(%rsi), %r8, %r15
	adox	%r8, %r14
	adcx	%r9, %r15
	adox	%r9, %r15

	mov	$0, %rdx
	test	%r15, %r15
	setns	%dl
	js	L(8)
	add	%rax, %rax
	adc	%r10, %r10
	adc	%r11, %r11
	adc	%rbx, %rbx
	adc	%rbp, %rbp
	adc	%r12, %r12
	adc	%r13, %r13
	adc	%r14, %r14
	adc	%r15, %r15
L(8):
	mov	%r10, 0*8(%rdi)
	mov	%r11, 1*8(%rdi)
	mov	%rbx, 2*8(%rdi)
	mov	%rbp, 3*8(%rdi)
	mov	%r12, 4*8(%rdi)
	mov	%r13, 5*8(%rdi)
	mov	%r14, 6*8(%rdi)
	mov	%r15, 7*8(%rdi)

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_normalised_9)
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15
	push	%rdi

	mov	%rdx, %rcx
	mov	0*8(%rdx), %rdx
	xor	%r9d, %r9d

	mulx	7*8(%rsi), %r8, %rax
	mulx	8*8(%rsi), %r8, %r10
	adcx	%r8, %rax
	adcx	%r9, %r10

	mov	1*8(%rcx), %rdx
	mulx	6*8(%rsi), %r8, %rbx
	mulx	7*8(%rsi), %r8, %r11
	adcx	%rbx, %rax
	adox	%r8, %rax
	adcx	%r11, %r10
	mulx	8*8(%rsi), %r8, %r11
	adox	%r8, %r10
	adcx	%r9, %r11
	adox	%r9, %r11

	mov	2*8(%rcx), %rdx
	mulx	5*8(%rsi), %r8, %rbp
	mulx	6*8(%rsi), %r8, %rbx
	adcx	%rbp, %rax
	adox	%r8, %rax
	adcx	%rbx, %r10
	mulx	7*8(%rsi), %r8, %rbx
	adox	%r8, %r10
	adcx	%rbx, %r11
	mulx	8*8(%rsi), %r8, %rbx
	adox	%r8, %r11
	adcx	%r9, %rbx
	adox	%r9, %rbx

	mov	3*8(%rcx), %rdx
	mulx	4*8(%rsi), %r8, %r12
	mulx	5*8(%rsi), %r8, %rbp
	adcx	%r12, %rax
	adox	%r8, %rax
	adcx	%rbp, %r10
	mulx	6*8(%rsi), %r8, %rbp
	adox	%r8, %r10
	adcx	%rbp, %r11
	mulx	7*8(%rsi), %r8, %rbp
	adox	%r8, %r11
	adcx	%rbp, %rbx
	mulx	8*8(%rsi), %r8, %rbp
	adox	%r8, %rbx
	adcx	%r9, %rbp
	adox	%r9, %rbp

	mov	4*8(%rcx), %rdx
	mulx	3*8(%rsi), %r8, %r13
	mulx	4*8(%rsi), %r8, %r12
	adcx	%r13, %rax
	adox	%r8, %rax
	adcx	%r12, %r10
	mulx	5*8(%rsi), %r8, %r12
	adox	%r8, %r10
	adcx	%r12, %r11
	mulx	6*8(%rsi), %r8, %r12
	adox	%r8, %r11
	adcx	%r12, %rbx
	mulx	7*8(%rsi), %r8, %r12
	adox	%r8, %rbx
	adcx	%r12, %rbp
	mulx	8*8(%rsi), %r8, %r12
	adox	%r8, %rbp
	adcx	%r9, %r12
	adox	%r9, %r12

	mov	5*8(%rcx), %rdx
	mulx	2*8(%rsi), %r8, %r14
	mulx	3*8(%rsi), %r8, %r13
	adcx	%r14, %rax
	adox	%r8, %rax
	adcx	%r13, %r10
	mulx	4*8(%rsi), %r8, %r13
	adox	%r8, %r10
	adcx	%r13, %r11
	mulx	5*8(%rsi), %r8, %r13
	adox	%r8, %r11
	adcx	%r13, %rbx
	mulx	6*8(%rsi), %r8, %r13
	adox	%r8, %rbx
	adcx	%r13, %rbp
	mulx	7*8(%rsi), %r8, %r13
	adox	%r8, %rbp
	adcx	%r13, %r12
	mulx	8*8(%rsi), %r8, %r13
	adox	%r8, %r12
	adcx	%r9, %r13
	adox	%r9, %r13

	mov	6*8(%rcx), %rdx
	mulx	1*8(%rsi), %r8, %r15
	mulx	2*8(%rsi), %r8, %r14
	adcx	%r15, %rax
	adox	%r8, %rax
	adcx	%r14, %r10
	mulx	3*8(%rsi), %r8, %r14
	adox	%r8, %r10
	adcx	%r14, %r11
	mulx	4*8(%rsi), %r8, %r14
	adox	%r8, %r11
	adcx	%r14, %rbx
	mulx	5*8(%rsi), %r8, %r14
	adox	%r8, %rbx
	adcx	%r14, %rbp
	mulx	6*8(%rsi), %r8, %r14
	adox	%r8, %rbp
	adcx	%r14, %r12
	mulx	7*8(%rsi), %r8, %r14
	adox	%r8, %r12
	adcx	%r14, %r13
	mulx	8*8(%rsi), %r8, %r14
	adox	%r8, %r13
	adcx	%r9, %r14
	adox	%r9, %r14

	mov	7*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %rdi
	mulx	1*8(%rsi), %r8, %r15
	adcx	%rdi, %rax
	adox	%r8, %rax
	adcx	%r15, %r10
	mulx	2*8(%rsi), %r8, %r15
	adox	%r8, %r10
	adcx	%r15, %r11
	mulx	3*8(%rsi), %r8, %r15
	adox	%r8, %r11
	adcx	%r15, %rbx
	mulx	4*8(%rsi), %r8, %r15
	adox	%r8, %rbx
	adcx	%r15, %rbp
	mulx	5*8(%rsi), %r8, %r15
	adox	%r8, %rbp
	adcx	%r15, %r12
	mulx	6*8(%rsi), %r8, %r15
	adox	%r8, %r12
	adcx	%r15, %r13
	mulx	7*8(%rsi), %r8, %r15
	adox	%r8, %r13
	adcx	%r15, %r14
	mulx	8*8(%rsi), %r8, %r15
	adox	%r8, %r14
	adcx	%r9, %r15
	adox	%r9, %r15

	mov	8*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %rdi
	adcx	%r8, %rax
	adcx	%rdi, %r10
	mulx	1*8(%rsi), %r8, %rdi
	adox	%r8, %r10
	adox	%rdi, %r11
	mulx	2*8(%rsi), %r8, %rdi
	adcx	%r8, %r11
	adcx	%rdi, %rbx
	mulx	3*8(%rsi), %r8, %rdi
	adox	%r8, %rbx
	adox	%rdi, %rbp
	mulx	4*8(%rsi), %r8, %rdi
	adcx	%r8, %rbp
	adcx	%rdi, %r12
	mulx	5*8(%rsi), %r8, %rdi
	adox	%r8, %r12
	adox	%rdi, %r13
	mulx	6*8(%rsi), %r8, %rdi
	adcx	%r8, %r13
	adcx	%rdi, %r14
	mulx	7*8(%rsi), %r8, %rdi
	adox	%r8, %r14
	adox	%rdi, %r15
	mulx	8*8(%rsi), %r8, %rdi
	adcx	%r8, %r15
	adcx	%r9, %rdi
	pop	%r8
	adox	%r9, %rdi

	mov	$0, %rdx
	test	%rdi, %rdi
	setns	%dl
	js	L(9)
	add	%rax, %rax
	adc	%r10, %r10
	adc	%r11, %r11
	adc	%rbx, %rbx
	adc	%rbp, %rbp
	adc	%r12, %r12
	adc	%r13, %r13
	adc	%r14, %r14
	adc	%r15, %r15
	adc	%rdi, %rdi
L(9):
	mov	%r10, 0*8(%r8)
	mov	%r11, 1*8(%r8)
	mov	%rbx, 2*8(%r8)
	mov	%rbp, 3*8(%r8)
	mov	%r12, 4*8(%r8)
	mov	%r13, 5*8(%r8)
	mov	%r14, 6*8(%r8)
	mov	%r15, 7*8(%r8)
	mov	%rdi, 8*8(%r8)

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()
