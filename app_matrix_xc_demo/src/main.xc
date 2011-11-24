/*
 * Parallel matrix operations XC demo
 *
 * Copyright (C) 2011 Steve Kerrison <github@stevekerrison.com>
 *
 * This software is freely distributable under a derivative of the
 * University of Illinois/NCSA Open Source License posted in
 * LICENSE.txt and at <http://github.xcore.com/>
 */
 
#include <stdio.h>
#include "matrix.h"

void perfcheck_sca(void);
void perfcheck_mul(void);

int main(void)
{
#if !defined(SKIP_BASIC_DEMO)
	MATRIX_CREATE(A,5,5,	1, 4, 3, 1, 5,
				1, 0, 3, 1, 3,
				4, 4, 4, 0, 2,
				3, 4, 3, 2, 2,
				2, 3, 3, 2, 2	);
		
	MATRIX_CREATE(B,5,5,	3, 3, 5, 4, 2,
				4, 4, 4, 3, 5,
				2, 1, 2, 2, 4,
				2, 2, 5, 2, 3,
				2, 1, 1, 4, 2	);
	MATRIX_CREATE(C,5,5,0);
	timer t;
	unsigned int tvA, tvB;
	int ops = 0;
	
	/* Time the operation in case we want to do some benchmarks :) */
	t :> tvA;
	/* In-place scalar multiplication */
	ops += matrix_sca_op(MUL,MATRIX_REF(A),11023,MATRIX_NULL(),0);
	ops += matrix_sca_op(MUL,MATRIX_REF(B),1798,MATRIX_NULL(),0);
	/* Matrix-multiplication into destination matrix */
	ops += matrix_mul(MATRIX_REF(A),MATRIX_REF(B),MATRIX_REF(C),0);
	t :> tvB;
	MATRIX_PRINT("A",A);
	MATRIX_PRINT("B",B);
	MATRIX_PRINT("C",C);
	printf("Ticks: %u, Ops: %u\n",tvB - tvA, ops);
#endif	
	/* That was a basic demo. Now let's compare performance */
	//perfcheck_sca();
	perfcheck_mul();
	return 0;
}


/* Some matrices used by the perfcheck_ functions, global to save stack */
#define MATRIX_MAX_PERF 48
MATRIX_CREATE(A,MATRIX_MAX_PERF,MATRIX_MAX_PERF,0);
MATRIX_CREATE(B,MATRIX_MAX_PERF,MATRIX_MAX_PERF,0);
MATRIX_CREATE(C,MATRIX_MAX_PERF,MATRIX_MAX_PERF,0);
MATRIX_CREATE(D,MATRIX_MAX_PERF,MATRIX_MAX_PERF,0);

void perfcheck_sca(void)
{
	timer t;
	unsigned int tvA, tvB;
	int e, lim;
	matrix_sca_op(RAND,MATRIX_REF(A),0,MATRIX_NULL(),0);
	matrix_sca_op(SET,MATRIX_REF(B),0,MATRIX_NULL(),0);
	matrix_sca_op(SET,MATRIX_REF(C),0,MATRIX_NULL(),0);
	printf("MATRIX_NTHREADS: %u\n",MATRIX_NTHREADS);
	for (int i = 1, ops = 0; i <= MATRIX_MAX_PERF; i++, ops = 0)
	{
		MATRIX_REDIM(A,i,i);
		printf("\nMatrix size: %dx%d\n",i,i);
		/* Scalar multiplication by 2 */
		lim = MATRIX_ROWS(A) * MATRIX_COLS(A);
		ops += lim;
		t :> tvA;
		for (e = 0; e < lim; e += 1)
		{
			B[e] = A[e] * 2;
		}
		t :> tvB;
		printf("By-hand scalar * 2 :: Ticks: %u, Ops: %u\n",tvB - tvA, ops);
		t :> tvA;
		ops = matrix_sca_op(MUL,MATRIX_REF(A),2,MATRIX_REF(C),0);
		t :> tvB;
		printf("sc_matrix scalar * 2 :: Ticks: %u, Ops: %u\n",tvB - tvA, ops);
		if (matrix_cmp(MATRIX_REF(B),MATRIX_REF(C)))
		{
			printf("ERROR: Results matrices do not match!\n");
			return;
		}
	}
	return;
}

void perfcheck_mul(void)
{
	timer t;
	unsigned int tvA, tvB;
	int e, lim;
	matrix_sca_op(RAND,MATRIX_REF(A),0,MATRIX_NULL(),0);
	matrix_sca_op(RAND,MATRIX_REF(B),0,MATRIX_NULL(),0);
	matrix_sca_op(SET,MATRIX_REF(C),0,MATRIX_NULL(),0);
	matrix_sca_op(SET,MATRIX_REF(D),0,MATRIX_NULL(),0);
	for (int i = 1, ops = 0; i <= MATRIX_MAX_PERF; i++, ops = 0)
	{
		MATRIX_REDIM(A,i,i);
		MATRIX_REDIM(B,i,i);
		printf("\nMatrix size: %dx%d\n",i,i);
		lim = MATRIX_ROWS(A) * MATRIX_COLS(A);
		ops += lim * (dimA[1] + 1);
		t :> tvA;
		for (int row = 0; row < dimA[0]; row += 1)
		{
			for (int col = 0; col < dimB[1]; col += 1)
			{
				int dst = row * dimB[1] + col;
				C[dst] = 0;
				for (e = 0; e < dimA[1]; e += 1)
				{
					C[dst] += A[row * dimA[1] + e] * B[dimB[1] * e + col];
				}
			}
		}
		t :> tvB;
		printf("By-hand A * B :: Ticks: %u, Ops: %u\n",tvB - tvA, ops);
		t :> tvA;
		ops = matrix_mul(MATRIX_REF(A),MATRIX_REF(B),MATRIX_REF(D),0);
		t :> tvB;
		printf("sc_matrix A * B :: Ticks: %u, Ops: %d\n",tvB - tvA, ops);
		if (matrix_cmp(MATRIX_REF(C),MATRIX_REF(D)))
		{
			printf("ERROR: Results matrices do not match!\n");
			return;
		}
	}
	return;
}
