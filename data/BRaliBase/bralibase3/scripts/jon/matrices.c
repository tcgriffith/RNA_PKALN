#include <stdio.h>#include <stdlib.h>#include <string.h>#include "matrices.h"#include "precomp_definitions.h"/*--------------------------------------------------------------------------------------------------||	AllocMatrix||	Allocate and set up a 'rows' x 'cols' matrix of arbitrary type.  Storage is allocated for both|	a vector of row pointers and the matrix itself.  The matrix is allocated as a single block, so|	that its elements may be referenced either as a[i][j] or as (*a)[cols*i + j].*/int AllocMatrix (VoidPtr pA, size_t elSize, int rows, int cols){	int			i;	unsigned	rowBytes;	char		*p, **a;	if (*(VoidPtr *)pA != NULL)		DeallocMatrix(pA);	if ((a = (char **)calloc(rows, sizeof(VoidPtr))) == NULL)		return ERROR;	*(VoidPtr *)pA = a;	rowBytes = (unsigned)(cols * elSize);	if ((a[0] = (char *)calloc(rows, rowBytes)) == NULL)		{		free(a);		return ERROR;		}	for (i = 0, p = *a; i < rows; i++, p += rowBytes)		a[i] = p;			return NO_ERROR;}/*--------------------------------------------------------------------------------------------------||	DeallocMatrix||	Deallocate memory for matrix allocated by AllocMatrix.*/void DeallocMatrix (VoidPtr pA){	char		**a;		if (*(VoidPtr *)pA != NULL)		{		a = (char **)(*(VoidPtr *)pA);		free(a[0]);		free(a);		*(VoidPtr *)pA = NULL;		}}/*--------------------------------------------------------------------------------------------------||	ShortMatrix||	Set row pointers into a 'rows' x 'cols' matrix whose storage is at 'buffer'.  This storage|	can be defined as matrix[M][N] in the caller or can be any buffer of size at least rows x cols.*/short **ShortMatrix (int rows, int cols, short **a, short *buffer){	int			i;	short		*p;	for (i = 0, p = buffer; i < rows; i++, p += cols)		a[i] = p;			return a;}/*--------------------------------------------------------------------------------------------------||	DoubleMatrix||	Set row pointers into a 'rows' x 'cols' matrix whose storage is at 'buffer'.  This storage|	can be defined as matrix[M][N] in the caller or can be any buffer of size at least rows x cols.*/double **DoubleMatrix (int rows, int cols, double **a, double *buffer){	int			i;	double		*p;	for (i = 0, p = buffer; i < rows; i++, p += cols)		a[i] = p;			return a;}/*--------------------------------------------------------------------------------------------------||	SetIdentityMatrix||	Initialize a matrix to the identity matrix.*/void SetIdentityMatrix (double **a, int n){	int			i, j;	for (i = 0; i < n; i++)		{		for (j = 0; j < n; j++)			a[i][j] = 0.0;		a[i][i] = 1.0;		}}/*--------------------------------------------------------------------------------------------------||	CopyDoubleMatrix||	Copy matrix 'a' to 'b'.*/void CopyDoubleMatrix (double **a, double **b, int m, int n){	(void)memcpy(b[0], a[0], sizeof(double) * m * n);}void ListDoubleMatrix (char *title, double **a, int m, int n){	int		i, j;		printf("%s\n", title);	for (i = 0; i < m; i++)		{		for (j = 0; j < n; j++)			printf("%18.10g", a[i][j]);		printf("\n");		}}void ListFloatMatrix (char *title, float **a, int m, int n){	int		i, j;		printf("%s\n", title);	for (i = 0; i < n; i++)		{		for (j = 0; j < n; j++)			printf("%18.10g", a[i][j]);		printf("\n");		}}void ListDoubleVector (char *title, double *a, int n){	int		i;		printf("%s", title);	for (i = 0; i < n; i++)		printf("%18.10g", a[i]);	printf("\n");}double **psdmatrix (int dim)// allocate a complex square matrix with subscript // range m[0..dim-1][0..dim-1]{	int		i;	double	**m;		/* allocate pointers to rows */	m = (double **) malloc((size_t)((dim)*sizeof(double*)));	if (!m)		{		printf("allocation error in pscmatrix 1.\n");		exit(1);		}	// allocate rows and set pointers to them 	m[0] = (double *) malloc((size_t)((dim*dim)*sizeof(double)));	if (!m[0])		{		printf("allocation error in pscmatrix 2.\n");		exit(1);		}	for(i=1;i<dim;i++)		{		m[i]  =m[i-1] + dim;		}	// return pointer to array of pointers to rows	return m;	}void free_psdmatrix (double **m)// free a double matrix allocated by psdmatrix ()// this works - I don't know why its a char *, maybe it should be a void *?{	free((char *) (m[0]));	free((char *) (m));	}void copy_psdmatrix (double **from, double **to, int dim){	int		row, col;		for (row = 0; row < dim; row++)		{		for (col = 0; col < dim; col++)			{			to[row][col] = from[row][col];			}		}	}