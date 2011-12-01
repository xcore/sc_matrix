/*
 * Matrix manipulation operations
 *
 * Copyright (C) 2011 Steve Kerrison <github@stevekerrison.com>
 *
 * This software is freely distributable under a derivative of the
 * University of Illinois/NCSA Open Source License posted in
 * LICENSE.txt and at <http://github.xcore.com/>
 */

#include <stdio.h> //For matrix_print
#include "matrix.h"
#include "matrix_worker.h"
#include "xs1.h"

/* Make invoking the scalar and array threads a little bit shorter */
#define MATRIX_WORKER_SPAWN(name,blockSize,lastBlock, ...) 		\
par									\
{									\
	par (int t = 0; t < MATRIX_NTHREADS-1; t++)			\
	{								\
		name(__VA_ARGS__+(t * sizeof(int)), blockSize * t, 	\
			blockSize);					\
	}								\
	name(__VA_ARGS__ +((MATRIX_NTHREADS-1) * sizeof(int)),		\
		blockSize * (MATRIX_NTHREADS-1), lastBlock);		\
}

{int,int} matrix_calc_block(int size, int nthreads)
{
	int blockSize, blockRem, lastBlock;
	blockSize = size / nthreads;
	blockRem = size % blockSize;
	if (blockRem != 0)
	{
		if ((blockSize + 1)*(nthreads-1) < size)
		{
			blockSize += 1;
			
		}
		lastBlock = size - (blockSize * (nthreads-1));
	}
	else
	{
		lastBlock = blockSize;
	}
	return {blockSize,lastBlock};
}

int matrix_redim(short dims[4],short rows, short columns)
{
	if (dims[2] < rows || dims[3] < columns)
	{
		return -1; //Larger than initial size
	}
	dims[0] = rows;
	dims[1] = columns;
	return rows * columns;
}

int matrix_sca_op(enum matrix_ops op, int A[], short dimA[2], int S,
	int ? C[], short ? dimC[2], char nThreads)
{
	int retval[8] = {0,0,0,0,0,0,0,0}, i;
	int ptA, ptC, ptDimA, ptRetval;
	int srcSize = dimA[0] * dimA[1], blockSize, lastBlock;
	POINTER(ptA,A);
	POINTER(ptC,C);
	POINTER(ptDimA,dimA);
	POINTER(ptRetval,retval);
	/* First do some sanity checks... */
	if (isnull(C))
	{
		//No checks yet
	}
	else
	{
		if (dimC[1] < dimA[1])
		{
			return -4; //Not enough columns in destination matrix
		}
		if (dimC[0] < dimA[0])
		{
			return -5; //Not enough rows in destination matrix
		}
	}
	if (isnull(C) || op == SET)
	{
		ptC = ptA;
	}
	/* Early-out for small data */
	if (srcSize < MATRIX_NTHREADS * MATRIX_NTHREADS)
	{
		switch (op)
		{
		case ADD:
			matrix_sca_worker_add(ptA,S,ptC,ptRetval,0,srcSize);
			break;
		case SUB:
			matrix_sca_worker_sub(ptA,S,ptC,ptRetval,0,srcSize);
			break;
		case MUL:
			matrix_sca_worker_mul(ptA,S,ptC,ptRetval,0,srcSize);
			break;
		case DIV:
		case SDIV:
			matrix_sca_worker_div(ptA,S,ptC,ptRetval,0,srcSize);
			break;
		case UDIV:
			matrix_sca_worker_udiv(ptA,S,ptC,ptRetval,0,srcSize);
			break;
		case RAND: //Fall through to default
			matrix_sca_worker_rand(ptC,ptRetval,0,srcSize);
			break;
		case SET:
			matrix_sca_worker_set(S,ptC,ptRetval,0,srcSize);
			break;
		case SHR:
			matrix_sca_worker_shr(ptA,S,ptC,ptRetval,0,srcSize);
			break;
		case ASHR:
			matrix_sca_worker_ashr(ptA,S,ptC,ptRetval,0,srcSize);
			break;
		case SHL:
			matrix_sca_worker_shl(ptA,S,ptC,ptRetval,0,srcSize);
			break;
		case AND:
			matrix_sca_worker_and(ptA,S,ptC,ptRetval,0,srcSize);
			break;
		default:
			break;	
		}
		return retval[0];
	}
	/* More optimal distribution of workload */
	{blockSize,lastBlock} = matrix_calc_block(srcSize,MATRIX_NTHREADS);
	switch (op)
	{
	case ADD:
		MATRIX_WORKER_SPAWN(matrix_sca_worker_add,blockSize,lastBlock,ptA,S,ptC,ptRetval);
		break;
	case SUB:
		MATRIX_WORKER_SPAWN(matrix_sca_worker_sub,blockSize,lastBlock,ptA,S,ptC,ptRetval);
		break;
	case MUL:
		MATRIX_WORKER_SPAWN(matrix_sca_worker_mul,blockSize,lastBlock,ptA,S,ptC,ptRetval);
		break;
	case DIV:
	case SDIV:
		MATRIX_WORKER_SPAWN(matrix_sca_worker_div,blockSize,lastBlock,ptA,S,ptC,ptRetval);
		break;
	case UDIV:
		MATRIX_WORKER_SPAWN(matrix_sca_worker_udiv,blockSize,lastBlock,ptA,S,ptC,ptRetval);
		break;
	case RAND:
		MATRIX_WORKER_SPAWN(matrix_sca_worker_rand,blockSize,lastBlock,ptC,ptRetval);
		break;
	case SET:
		MATRIX_WORKER_SPAWN(matrix_sca_worker_set,blockSize,lastBlock,S,ptC,ptRetval);
		break;
	case SHR:
		MATRIX_WORKER_SPAWN(matrix_sca_worker_shr,blockSize,lastBlock,ptA,S,ptC,ptRetval);
		break;
	case ASHR:
		MATRIX_WORKER_SPAWN(matrix_sca_worker_ashr,blockSize,lastBlock,ptA,S,ptC,ptRetval);
		break;
	case SHL:
		MATRIX_WORKER_SPAWN(matrix_sca_worker_shl,blockSize,lastBlock,ptA,S,ptC,ptRetval);
		break;
	case AND:
		MATRIX_WORKER_SPAWN(matrix_sca_worker_and,blockSize,lastBlock,ptA,S,ptC,ptRetval);
		break;
	default:
		break;	
	}
	for (i = 1; i < 8; i++)
	{
		retval[0] += retval[i];
	}
	return retval[0];
}

int matrix_arr_op(enum matrix_ops op, int A[], short dimA[2], int B[], short dimB[2],
	int ? C[], short ? dimC[2], char nThreads)
{
	int retval[8] = {0,0,0,0,0,0,0,0}, i;
	int ptA, ptB, ptC, ptDimA, ptDimB, ptRetval;
	int srcSize = dimA[0] * dimA[1], blockSize, lastBlock;
	POINTER(ptA,A);
	POINTER(ptB,B);
	POINTER(ptC,C);
	POINTER(ptDimA,dimA);
	POINTER(ptDimB,dimB);
	POINTER(ptRetval,retval);
	/*int ptA = pointer_int(A), ptB = pointer_int(B), ptC = pointer_int(C),
		ptDimA = pointer_short(dimA), ptDimB = pointer_short(dimB),
		ptRetval = pointer_int(retval);*/
	/* First do some sanity checks... */
	if (dimA[0] != dimB[0] || dimA[1] != dimB[1]) return -2; //Invalid dimensions
	if (isnull(C))
	{
		//No checks yet
	}
	else
	{
		if (dimC[1] < dimB[1])
		{
			return -4; //Not enough columns in destination matrix
		}
		if (dimC[0] < dimA[0])
		{
			return -5; //Not enough rows in destination matrix
		}
	}
	if (isnull(C))
	{
		ptC = ptA;
	}
	/* Early-out for small data */
	if (srcSize < MATRIX_NTHREADS * MATRIX_NTHREADS)
	{
		switch (op)
		{
		case ADD:
			matrix_arr_worker_add(ptA,ptB,ptC,ptRetval,0,srcSize);
			break;
		case SUB:
			matrix_arr_worker_sub(ptA,ptB,ptC,ptRetval,0,srcSize);
			break;
		case MUL:
			matrix_arr_worker_mul(ptA,ptB,ptC,ptRetval,0,srcSize);
			break;
		case DIV:
		case SDIV:
			matrix_arr_worker_div(ptA,ptB,ptC,ptRetval,0,srcSize);
			break;
		case UDIV:
			matrix_arr_worker_udiv(ptA,ptB,ptC,ptRetval,0,srcSize);
			break;
		case AND:
			matrix_arr_worker_and(ptA,ptB,ptC,ptRetval,0,srcSize);
			break;
		case RAND: //Fall through to default
		case SET:  //Still falling...
		case SHR:  //It's a long way down...
		case ASHR:
		case SHL:
		default:
			break;	
		}
		return retval[0];
	}
	/* More optimal distribution of workload */
	{blockSize,lastBlock} = matrix_calc_block(srcSize,MATRIX_NTHREADS);
	switch (op)
	{
	case ADD:
		MATRIX_WORKER_SPAWN(matrix_arr_worker_add,blockSize,lastBlock,ptA,ptB,ptC,ptRetval);
		break;
	case SUB:
		MATRIX_WORKER_SPAWN(matrix_arr_worker_sub,blockSize,lastBlock,ptA,ptB,ptC,ptRetval);
		break;
	case MUL:
		MATRIX_WORKER_SPAWN(matrix_arr_worker_mul,blockSize,lastBlock,ptA,ptB,ptC,ptRetval);
		break;
	case DIV:
	case SDIV:
		MATRIX_WORKER_SPAWN(matrix_arr_worker_div,blockSize,lastBlock,ptA,ptB,ptC,ptRetval);
		break;
	case UDIV:
		MATRIX_WORKER_SPAWN(matrix_arr_worker_udiv,blockSize,lastBlock,ptA,ptB,ptC,ptRetval);
		break;
	case AND:
		MATRIX_WORKER_SPAWN(matrix_arr_worker_and,blockSize,lastBlock,ptA,ptB,ptC,ptRetval);
		break;
	case RAND: //Fall through to default
	case SET:  //Still falling...
	case SHR:  //It's a long way down...
	case ASHR:
	case SHL:
	default:
		break;	
	}
	for (i = 1; i < 8; i++)
	{
		retval[0] += retval[i];
	}
	return retval[0];
}

int matrix_mul(int A[], short dimA[2], int B[], short dimB[2],
	int ? C[], short ? dimC[2], char nThreads)
{
	int retval[8] = {0,0,0,0,0,0,0,0}, i;
	int ptA, ptB, ptC, ptDimA, ptDimB, ptRetval;
	int dstSize = dimA[0]* dimB[1], blockSize, lastBlock;
	POINTER(ptA,A);
	POINTER(ptB,B);
	POINTER(ptC,C);
	POINTER(ptDimA,dimA);
	POINTER(ptDimB,dimB);
	POINTER(ptRetval,retval);
	/* First do some sanity checks... */
	if (dimA[1] != dimB[0]) return -2; //Matrices cannot be multiplied
	if (isnull(C))
	{
		if (dimA[1] < dimB[1])
			return -3; //Inline result cannot fit in matrix A.
	}
	else
	{
		if (dimC[1] < dimB[1])
		{
			return -4; //Not enough columns in destination matrix
		}
		if (dimC[0] < dimA[0])
		{
			return -5; //Not enough rows in destination matrix
		}
	}
	if (isnull(C))
	{
		//FIXME - Use a thread-safe strategy for in-place results
		return -1; //In-place result not supported at the moment
	}
	
	/* Small matrix, use a single thread... */
	if (dstSize < MATRIX_NTHREADS * MATRIX_NTHREADS)
	{
		matrix_mul_worker(ptA,ptDimA,ptB,ptDimB,ptC,0,dimA[0] * dimB[1],ptRetval);
		return retval[0];
	}
	{blockSize,lastBlock} = matrix_calc_block(dstSize,MATRIX_NTHREADS);
	par
	{
		par (int t = 0; t < MATRIX_NTHREADS-1; t++)
		{
			matrix_mul_worker(ptA,ptDimA,ptB,ptDimB,ptC,
				blockSize * t, blockSize,
				ptRetval + t * sizeof(int));
		}
		matrix_mul_worker(ptA,ptDimA,ptB,ptDimB,ptC,
			blockSize * (MATRIX_NTHREADS-1), lastBlock,
			ptRetval + (MATRIX_NTHREADS-1) * sizeof(int));
	}
	for (i = 1; i < 8; i++)
	{
		retval[0] += retval[i];
	}
	return retval[0];
}

int matrix_cmp(int A[], short dimA[2], int B[], short dimB[2])
{
	int ret = 0, size, e;
	if (dimA[0] != dimB[0] || dimA[1] != dimB[1])
	{
		return -2; /* Can only compare like-dim'd matrices */
	}
	size = dimA[0] * dimA[1];
	for (e = 0; e < size; e += 1)
	{
		ret += A[e] != B[e];
	}
	return ret;
}

void matrix_print(char name[], int M[], short dimM[2])
{
	int r,c,s = dimM[0] * dimM[1];
	printf("Matrix %s =\n", name);
	for (r = 0; r < s; r += dimM[1])
	{
		printf(" ");
		for (c = 0; c < dimM[1]; c++)
		{
			printf(" %d ",M[r + c]);
		}
		printf("\n");
	}
	printf("\n");
}

