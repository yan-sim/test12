#include "matrix.h"

#define SUCCESS 0
#define FAILURE 1


struct matrix_t
{
    size_t n_rows;
    size_t n_cols;
    float* data;
};

/** 
*   MatCreate
*   ------------
*   Creates a new matrix according to the user specifications.
*
*   Params
*   ------
*   rows - number of rows in the matrix.
*	cols - number of columns in the matrix.
*	data - an array of the elements that will be inserted into the matrix.
*   if data is NULL, a zero matrix will be created.
*
*   Return
*   ------
*   A pointer to the new matrix. NULL on failure.
* 
*   Asserts: n_rows > 0, n_cols > 0
*   Errors: malloc failure
*   Undefined Behaviour: data not allocated enough space or one of the dims are 0.
*/
matrix_t* MatCreate(size_t n_rows, size_t n_cols, const float* data)
{
  matrix_t* mat = (matrix_t*)malloc(sizeof(matrix_t));
  size_t total_elements = n_rows * n_cols;
  if (!mat)
  {
    return NULL;
  }
  assert(n_rows > 0 && n_cols > 0);
  mat->data = (float*)malloc(n_rows * n_cols* sizeof(float));
  if (!mat->data)
  {
    free(mat);
    return NULL;
  }
  if(data == NULL)
  {
    memset(mat->data, 0, total_elements * sizeof(float));
  }
  else
  {
   memcpy(mat->data, data, total_elements * sizeof(float));
  }
  mat->n_rows = n_rows;
  mat->n_cols = n_cols;
 
  return mat;
}



/** 
*   MatDestroy
*   ----------
*   Clears a matrix from memory.
*   
*   Params
*   ------
*   mat - a pointer to the matrix.
* 
*   Asserts:--
*   Errors:--
*   Undefined Behaviour:--
*/
void MatDestroy(matrix_t* mat)
{
  free(mat->data);
  free(mat);
}

/** 
*   MatSubmatrix
*   ------------
*   Creates a new submatrix by removing a row and col from a given matrix.
*
*   Params
*   ------
*   mat - a pointer to the matrix.
*   row - the row which won't be included in the submatrix.
*	cols - the column which won't be included in the submatrix.
*
*   Return
*   ------
*   A pointer to the submatrix. NULL on failure.
* 
*   Asserts: row/col out of dimensions
*   Errors: malloc failure
*   Undefined Behaviour: row/col out of dimensions
* 
*/
matrix_t* MatSubmatrix(const matrix_t* mat, size_t row, size_t col)
{
  matrix_t* new_mat = MatCreate(mat->n_rows - 1, mat->n_cols - 1, NULL);
  size_t i = 0, j = 0;
  float *src_p = mat->data, *dst_p = new_mat->data;
  if (new_mat == NULL || new_mat->data == NULL)
  {
    free(new_mat);
    return NULL;
  }
  for (i = 0; i < mat->n_rows; ++i)
  {
    if (i == row)
    {
      src_p += mat->n_cols;
      continue;
    }
    for (j = 0; j < n_cols; ++j)
    {
      if (j == col)
      {
        src_p++;
        continue;
      }
      *dst_p = *src_p;
      dst_p++;
      src_p++;     
    }
  }
  return new_mat;
}

/** 
*   MatI
*   ----
*   Creates an identity matrix of size n * n.
*
*   Params
*   ------
*   n - size of created identity matrix.
*
*   Return
*   ------
*   A pointer to the created identity matrix. NULL on failure.
* 
* 
*   Asserts: n > 0
*   Errors: malloc failure
*   Undefined Behaviour: n == 0
*/
matrix_t* MatI(size_t n)
{
  matrix_t* mat = MatCreate(n, n, NULL);
  float* dst_p = mat->data;
  size_t i = 0, j = o;
  if (n ==0)
  {
    return NULL;
  }
  if (!mat)
  {
    return NULL;
  }
  for (i = 0; i < n; ++i)
  {
    for (j = 0; j < n; ++j)
    {
      if (i == j)
      {
        *dst_p = 1.0f;
      }
      else
      {
        *dst_p = 0.0f;
      }
      dst_p++;
    }
  }
  return mat;  
}


/** 
*   MatAdd
*   ------
*   Adds two matrices.
*
*   Params
*   ------
*   mat1 - the 1st matrix of the calculation.
*   mat2 - the 2nd matrix of the calculation.
*
*   Return
*   ------
*   A pointer to the result matrix.
*   Returns NULL if the matrices dimensions don't match.
* 
* 
*   Asserts: mat1, mat2 not null
*   Errors: malloc failure, matrices dimensions should be the same.
*   Undefined Behaviour: null inputs.
* 
*/
matrix_t* MatAdd(const matrix_t* mat1, const matrix_t* mat2)
{
  matrix_t* new_mat;
  float* src1_p, *src2_p, *dst_p;
  size_t i = 0, j = 0;
  assert(mat1 != NULL && mat2 != NULL);
  if (mat1->n_rows != mat2->n_rows || mat1->n_cols != mat2->n_cols)
  {
    return NULL;
  }
  new_mat = MatCreate(mat1->n_rows, mat1->n_cols, NULL);
  if (result == NULL)
  {
    return NULL;
  }
  src1_p = mat1->data;
  src2_p = mat2->data;
  dst_p = new_mat->data;
  for (i = 0; i < mat1->n_rows; ++i)
  {
    for (j = 0; j < mat1->n_cols; ++j)
    {
      *dst_p = *src1_p + *src2_p;
      src1_p++;
      src2_p++;
      dst_p++;
    }
  }
  return new_mat;
}

/** 
*   MatMult
*   -------
*   Multiplies two matrices.
*
*   Params
*   ------
*   mat1 - the 1st matrix of the calculation.
*	mat2 - the 2nd matrix of the calculation.
*
*   Return
*   ------
*   A pointer to the result matrix.
*   Returns NULL on failure.
* 
*   Asserts: one or two of the matrices are NULL.
*   Errors: mat1 rows need to match mat2 columns.
*   Undefined Behaviour:  NULL as one of the matrices values.
*/
matrix_t* MatMult(const matrix_t* mat1, const matrix_t* mat2)
{
  matrix_t* new_mat;
  float* mat1_p, *mat2_p, *resukt_p, sum;
  size_t i = 0, j = 0, k = 0;
  assert(mat1 != NULL && mat2 != NULL);
  if (mat1->n_cols != mat2->n_rows)
  {
    return NULL;
  }
  new_mat = MatCreate(mat1->n_rows, mat2->n_cols, NULL);
  if (result == NULL)
  {
    return NULL;
  }
  mat1_p = mat1->data;
  mat2_p = mat2->data;
  result_p = new_mat->data;
  for (i = 0; i < mat1->n_rows; ++i)
  {
    for (j = 0; j < mat2->n_cols; ++j)
    {
      sum = 0.0f;

      for (k = 0; k < mat1->n_cols; ++k)
      {
         sum += *(mat1_p + i * mat1->n_cols + k) * *(mat2_p + k * mat2->n_cols + j);
      }
      *result_p = sum;
      result_p++;
    }
  }
    return new_mat;
}



/**
*
*   Params
*   ------
*   mat - a pointer to the matrix.
*	scalar - the scalar.
*
*   Return
*   ------
*   A pointer to the newly created matrix that holds the result.NULL on failure.
* 
* Asserts: matrix is NULL, scalar larger or equal to the percision rate.
* Errors: 
* Undefined Behaviour: scalar beyond percision rate.
*/
matrix_t* MatScalarMult(matrix_t* mat, float scalar)
{
  matrix_t* new_mat;
  float* src_p, *dst_p;
  size_t i = 0, j = 0;
  assert(mat != NULL);
  new_mat = MatCreate(mat->n_rows, mat->n_cols, NULL);
  if (new_mat == NULL)
  {
    return NULL;
  }
  src_pt = mat->data;
  dst_p = new_mat->data;
  for (i = 0; i < mat->n_rows; ++i)
  {
    for (j = 0; j < mat->n_cols; ++j)
    {
     *dst_p = *src_p * scalar;
      src_p++;
      dst_p++;
    }
  }
  return new_mat;
}






/** 
*   MatMult
*   -------
*   Compares mat1 to mat2.
*
*   Params
*   ------
*   mat1 - the 1st matrix of the comparison
*	mat1 - the 2nd matrix of the comparison
*
*   Return
*   ------
*   0 if the matrices are equal, nonzero otherwise.
* 
* 
*   Asserts: mat1 or mat2 is NULL
*   Errors: 
*   Undefined Behaviour: mat1 or mat2 is NULL
* 
*   check dimentions before checking the data
* 
*/
int MatCompare(const matrix_t* mat1, const matrix_t* mat2);

/**   
*  MatTranspose
*  ------------
*  flips a matrix over its diagonal
* 
*  Params
*  ------
*  mat - a pointer to the matrix.
* 
*  Return
*  ------
*  Transposed matrix. NULL on failure
*
* 
* 
*   Asserts: mat is null
*   Errors: malloc failure (assuming we create a new matrix)
*   Undefined Behaviour: 
*/
matrix_t* MatTranspose(const matrix_t* mat);

/** 
*   MatTrace
*   -------
*   Calculates the trace of the given matrix.
*
*   Params
*   ------
*   mat - a pointer to the matrix.
*
*   Return
*   ------
*   Trace of given matrix.
* 
*   Asserts:
*   Errors: 
*   Undefined Behaviour: (number or rows)!=(number of cols),
* 
*/
float MatTrace(const matrix_t* mat);

/** 
*   MatInvert
*   ---------
*   returns the inverse of a matrix.
*
*   Params
*   ------
*   mat - a pointer to the matrix.
*
*   Return
*   ------
*   pointer to inverted matrix. NULL on failure.
*   Undefined behavior if the matrix is not square.
*
*
* NUll pointer to mat
* matrix not square
* matrix not invertible
*
*
*   Asserts: pointer is not null, matrix is square
*   Errors: malloc failure when creating a new matrix
*   Undefined Behaviour: matrix is not square
*/
matrix_t* MatInvert(const matrix_t* mat);

/** 
*   MatDet
*   ------
*   calculates the determinant of given matrix.
*
*   Params
*   ------
*   mat - a pointer to the matrix.
*
*   Return
*   ------
*   float - value of matrices determinant.
*   Undefined behavior if the matrix is not square.
*
*
* NUll pointer to mat
* matrix not square
*/
float MatDet(const matrix_t* mat);

/** 
*   MatNorm
*   --------
*   calculates the frobinius norm (l2) of the matrix.
*
*   Params
*   ------
*   mat - a pointer to the matrix.
* 
*   Return
*   ------
*   The norm of the matrix.
*
*
*
* NUll pointer to mat
* 
*   Asserts:
*   Errors: 
*   Undefined Behaviour: mat is NULL
* 
*/
float MatNorm(const matrix_t* mat);

/** 
*   MatShape
*   --------
*   Get the shape of a matrix.
*
*   Params
*   ------
*   mat - a pointer to the matrix.
*   dims - a buffer where the result will be stored.
*   dims[0] will be the number of rows.
*   dims[1] will be the number of columns.
*
*
*
* 
* 
* NUll pointer to mat
*/
void MatShape(const matrix_t* mat, size_t dims[2]);

/** 
*   MatGetElem
*   ----------
*   Get a value from a matrix.
*
*   Params
*   ------
*   mat - a pointer to the matrix.
*   row - the value's row.
*   col - the value's column.
*
*   Return
*   ------
*   The value at the intersection of the row and column.
*
*
* 
* index outside of mat
* NUll pointer to mat
*/
float MatGetElem(const matrix_t* mat, size_t row, size_t col)
{
  float* mat_p = mat->data;
  assert(mat != NULL);
  if (row >= mat->n_rows || col >= mat->n_cols)
  {
    return -99999.9f;
  }
  return *(mat_p + row * mat->n_cols + col);
}




int test1(void)
{
  float* data = {1,2,3,1,2,3,1,2,3};
  matrix_t* mat = MatCreate(3, 3, data);
  if (mat == NULL)
  {
    return FAILURE;
  }
  MatDestroy(mat);
  return SUCCESS;
}


int main()
{
  if (SUCCESS = test1())
  {
    printf("test 1 passed!")
  }
  else
  {
  printf("test 1 failed! MatCreate/Destroy")
  }
  
  
  
  
}



























