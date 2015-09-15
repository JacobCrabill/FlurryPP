/*!
 * \file matrix.cpp
 * \brief Class for simplified matrix storage & manipulation
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 0.0.1
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2014 Jacob Crabill.
 *
 */
#include "../include/matrix.hpp"

#include <set>

template<typename T, uint N>
Array<T,N>::Array()
{
  data.resize(0);
  dims = {{0,0,0,0}};
}

template<typename T, uint N>
Array<T,N>::Array(uint inNDim0, uint inNDim1, uint inDim2, uint inDim3)
{
  data.resize(inNDim0*inNDim1*inDim2*inDim3);
  dims = {{inNDim0,inNDim1,inDim2,inDim3}};
}

template<typename T>
Array2D<T>::Array2D()
{
  this->data.resize(0);
  this->dims[0] = 0;
  this->dims[1] = 0;
}

template<typename T>
Array2D<T>::Array2D(uint inNDim0, uint inNDim1)
{
  this->data.resize(inNDim0*inNDim1);
  this->dims[0] = inNDim0;
  this->dims[1] = inNDim1;
}

template<typename T>
matrix<T>::matrix()
{
  this->data.resize(0);
  this->dims[0] = 0;
  this->dims[1] = 0;
}

template<typename T>
matrix<T>::matrix(uint inNDim0, uint inNDim1)
{
  this->data.resize(inNDim0*inNDim1);
  this->dims[0] = inNDim0;
  this->dims[1] = inNDim1;
}

template<typename T,uint N>
Array<T,N>::Array(const Array<T,N> &inMatrix)
{
  data = inMatrix.data;
  dims = inMatrix.dims;
}

template<typename T, uint N>
Array<T,N> Array<T,N>::operator=(const Array<T,N> &inMatrix)
{
  data = inMatrix.data;
  dims = inMatrix.dims;
  return *this;
}

template<typename T, uint N>
void Array<T,N>::setup(uint inDim0, uint inDim1, uint inDim2, uint inDim3)
{
  dims = {{inDim0,inDim1,inDim2,inDim3}};
  data.resize(inDim0*inDim1*inDim2*inDim3);
}


template<typename T>
void matrix<T>::addMatrix(matrix<T> &A, double a)
{
  if (A.dims[0] != this->dims[0] || A.dims[1] != this->dims[1])
    FatalErrorST("Incompatible matrix sizes for addMatrix.");

  for (uint i=0; i<this->dims[0]; i++)
    for (uint j=0; j<this->dims[1]; j++)
      this->data[i*this->dims[1]+j] += a*A(i,j);
}

template<>
matrix<double>& matrix<double>::operator+=(matrix<double> &A)
{
  if (A.dims[0] != this->dims[0] || A.dims[1] != this->dims[1])
    FatalErrorST("Incompatible matrix sizes for addMatrix.");

  for (uint i=0; i<this->dims[0]; i++)
    for (uint j=0; j<this->dims[1]; j++)
      this->data[i*this->dims[1]+j] += A(i,j);

  return *this;
}

template<>
matrix<double>& matrix<double>::operator-=(matrix<double> &A)
{
  if (A.dims[0] != this->dims[0] || A.dims[1] != this->dims[1])
    FatalErrorST("Incompatible matrix sizes for addMatrix.");

  for (uint i=0; i<this->dims[0]; i++)
    for (uint j=0; j<this->dims[1]; j++)
      this->data[i*this->dims[1]+j] -= A(i,j);

  return *this;
}

template<typename T, uint N>
T* Array<T,N>::operator[](int inRow)
{
  if (inRow < (int)this->dims[0] && inRow >= 0) {
    return &data[inRow*dims[1]];
  }
  else {
    cout << "inRow = " << inRow << ", nRows = " << this->dims[0] << endl;
    FatalErrorST("Operator[]: Attempted out-of-bounds access in matrix.");
  }
}

template<typename T, uint N>
T& Array<T,N>::operator()(int i, int j, int k, int l)
{
  if (i<(int)this->dims[0] && i>=0 && j<(int)this->dims[1] && j>=0 &&
      k<(int)this->dims[2] && k>=0 && l<(int)this->dims[3] && l>= 0)
  {
    return data[l+dims[3]*(k+dims[2]*(j+dims[1]*i))];
  }
  else {
    cout << "i=" << i << ", dim0=" << dims[0] << ", j=" << j << ", dim1=" << dims[1] << ", ";
    cout << "k=" << k << ", dim2=" << dims[2] << ", l=" << l << ", dim3=" << dims[3] << endl;
    FatalErrorST("Attempted out-of-bounds access in Array.");
  }
}

template<typename T>
T& Array2D<T>::operator()(int i, int j)
{
  if (i<(int)this->dims[0] && i>=0 && j<(int)this->dims[1] && j>=0) {
    return this->data[j+this->dims[1]*i];
  }
  else {
    cout << "i=" << i << ", dim0=" << this->dims[0] << ", j=" << j << ", dim1=" << this->dims[1] << endl;
    FatalErrorST("Attempted out-of-bounds access in matrix.");
  }
}

template<typename T>
void matrix<T>::initializeToZero(void)
{
  for (uint i=0; i<this->dims[0]; i++)
    for (uint j=0; j<this->dims[1]; j++)
      for (uint k=0; k<this->dims[2]; k++)
        for (uint l=0; l<this->dims[3]; l++)
          this->data[l+this->dims[3]*(k+this->dims[2]*(j+this->dims[1]*i))] = 0;
}

template<typename T>
void matrix<T>::initializeToValue(T val)
{
  for (uint i=0; i<this->dims[0]; i++)
    for (uint j=0; j<this->dims[1]; j++)
      for (uint k=0; k<this->dims[2]; k++)
        for (uint l=0; l<this->dims[3]; l++)
          this->data[l+this->dims[3]*(k+this->dims[2]*(j+this->dims[1]*i))] = val;
}

template <typename T>
void matrix<T>::timesMatrix(matrix<T> &A, matrix<T> &B)
{
  uint i, j, k, p;
  p = A.dims[1];

  if (A.dims[0] != this->dims[1]) FatalErrorST("Incompatible matrix sizes in matrix multiplication!");
  if (B.dims[0] != this->dims[0] || B.dims[1] != A.dims[1]) B.setup(this->dims[0], A.dims[1]);

  B.initializeToZero();

  for (i=0; i<this->dims[0]; i++) {
    for (j=0; j<this->dims[1]; j++) {
      for (k=0; k<p; k++) {
        B[i][k] += this->data[i*this->dims[1]+j]*A[j][k];
      }
    }
  }
}

/*matrix<int> operator*(matrix<int> &A, matrix<int> &B)
{
  FatalError("Not expecting to use operator* for integer matrices.");
}

matrix<double*> operator*(matrix<double*> &A, matrix<double*> &B)
{
  FatalError("operator* not supported for non-arithmatic data types.");
}

matrix<double> operator*(matrix<double> &A, matrix<double> &B)
{
  uint p = A.dims[1];

  if (A.dims[1] != B.dims[0]) FatalError("Incompatible matrix sizes in matrix multiplication!");

  matrix<double> out(A.dims[0],B.dims[1]);
  out.initializeToZero();

  for (uint i=0; i<A.dims[0]; i++) {
    for (uint j=0; j<B.dims[1]; j++) {
      for (uint k=0; k<p; k++) {
        out(i,j) += A(i,k)*B(k,j);
      }
    }
  }
}*/


template <typename T>
void matrix<T>::timesMatrixPlus(matrix<T> &A, matrix<T> &B)
{
  uint i, j, k, p;
  p = A.dims[1];

  if (A.dims[0] != this->dims[1]) FatalErrorST("Incompatible matrix sizes in matrix multiplication!");
  if (B.dims[0] != this->dims[0] || B.dims[1] != A.dims[1]) B.setup(this->dims[0], A.dims[1]);

  for (i=0; i<this->dims[0]; i++) {
    for (j=0; j<this->dims[1]; j++) {
      for (k=0; k<p; k++) {
        B[i][k] += this->data[i*this->dims[1]+j]*A[j][k];
      }
    }
  }
}

template <typename T>
void matrix<T>::timesVector(vector<T> &A, vector<T> &B)
{
  uint i, j;

  //if (A.size() != this->dims[1]) FatalError("Incompatible vector size");
  if (B.size() != this->dims[1]) B.resize(this->dims[1]);

  for (i=0; i<this->dims[0]; i++) {
    B[i] = 0;
    for (j=0; j<this->dims[1]; j++) {
      B[i] += this->data[i*this->dims[1]+j]*A[j];
    }
  }
}

template<typename T>
void Array2D<T>::appendRows(Array2D<T> &mat)
{
  if (this->dims[1]!= 0 && mat.getDim1()!=this->dims[1])
    FatalErrorST("Attempting to append rows of wrong size to matrix.");

  this->data.insert(this->data.end(),mat.getData(),mat.getData()+mat.getSize());

  if (this->dims[1]==0) this->dims[1]=mat.getDim1();
  this->dims[0]+=mat.getDim0();
}

template<typename T>
void Array2D<T>::insertRow(const vector<T> &vec, int rowNum)
{
  if (this->dims[1]!= 0 && vec.size()!=this->dims[1])
    FatalErrorST("Attempting to assign row of wrong size to matrix.");

  if (rowNum==INSERT_AT_END || rowNum==(int)this->dims[0]) {
    // Default action - add to end
    this->data.insert(this->data.end(),vec.begin(),vec.end());
  }else{
    // Insert at specified location
    this->data.insert(this->data.begin()+rowNum*this->dims[1],vec.begin(),vec.end());
  }

  if (this->dims[1]==0) this->dims[1]=vec.size(); // This may not be needed (i.e. may never have this->dims[1]==0). need to verify how I set up this->dims[0], this->dims[1]...
  this->dims[0]++;
}

template<typename T>
void Array2D<T>::insertRow(T *vec, int rowNum, int length)
{
  if (this->dims[1]!=0 && length!=(int)this->dims[1])
    FatalErrorST("Attempting to assign row of wrong size to matrix.");

  if (rowNum==INSERT_AT_END || rowNum==(int)this->dims[0]) {
    // Default action - add to end
    this->data.insert(this->data.end(),vec,vec+length);
  }else{
    // Insert at specified location
    this->data.insert(this->data.begin()+rowNum*this->dims[1],vec,vec+length);
  }

  if (this->dims[0]==0)
    this->dims = {{0,(uint)length,1,1}};

  this->dims[0]++;
}


template<typename T>
void Array2D<T>::insertRowUnsized(const vector<T> &vec)
{
  // Add row to end, and resize matrix (add columns) if needed
  if (vec.size() > this->dims[1]) addCols(vec.size()-this->dims[1]);

  this->data.insert(this->data.end(),vec.begin(),vec.end());

  // If row too short, finish filling with 0's
  if (vec.size() < this->dims[1]) this->data.insert(this->data.end(),this->dims[1]-vec.size(),(T)0);

  this->dims[0]++;
}

template<typename T>
void Array2D<T>::insertRowUnsized(T* vec, uint length)
{
  // Add row to end, and resize matrix (add columns) if needed
  if (length > this->dims[1]) addCols(length-this->dims[1]);

  this->data.insert(this->data.end(),vec,vec+length);

  // If row too short, finish filling with 0's
  if (length < this->dims[1]) this->data.insert(this->data.end(),this->dims[1]-length,(T)0);

  this->dims[0]++;
}

template<typename T>
void Array2D<T>::addCol(void)
{
  typename vector<T>::iterator it;
  for (uint row=0; row<this->dims[0]; row++) {
    it = this->data.begin() + (row+1)*(this->dims[1]+1) - 1;
    this->data.insert(it,1,(T)0);
  }
  this->dims[1]++;
}

template<typename T>
void Array2D<T>::addCols(int nCols)
{
  typename vector<T>::iterator it;
  for (uint row=0; row<this->dims[0]; row++) {
    it = this->data.begin() + (row+1)*(this->dims[1]+nCols) - nCols;
    this->data.insert(it,nCols,(T)0);
  }
  this->dims[1] += nCols;
}

template<typename T>
void Array2D<T>::removeCols(int nCols)
{
  if (nCols == 0) return;

  typename vector<T>::iterator it;
  for (uint row=this->dims[0]; row>0; row--) {
    it = this->data.begin() + row*this->dims[1];
    this->data.erase(it-nCols,it);
  }
  this->dims[1] -= nCols;
}

template<typename T>
vector<T> Array2D<T>::getRow(uint row)
{
  if (row > this->dims[0]) FatalErrorST("Attempting to grab row beyond end of matrix.");

  vector<T> out;
  out.assign(&(this->data[row*this->dims[1]]),&(this->data[row*this->dims[1]])+this->dims[1]);
  return out;
}

template<typename T>
Array2D<T> Array2D<T>::getRows(vector<int> ind)
{
  matrix<T> out;
  for (auto& i:ind) out.insertRow(&(this->data[i*this->dims[1]]),-1,this->dims[1]);
  return out;
}

template<typename T>
vector<T> Array2D<T>::getCol(int col)
{
  vector<T> out;
  for (uint i=0; i<this->dims[0]; i++) out.push_back(this->data[i*this->dims[1]+col]);
  return  out;
}

template<typename T>
void matrix<T>::print()
{
  cout.precision(8);
  for (uint i=0; i<this->dims[0]; i++) {
    for (uint j=0; j<this->dims[1]; j++) {
      cout << setw(10) << left << this->data[i*this->dims[1]+j] << " ";
    }
    cout << endl;
  }
}

template<typename T>
bool matrix<T>::checkNan(void)
{
  for (auto& i:this->data)
    if (std::isnan(i)) return true;

  return false;
}

template<typename T>
void matrix<T>::unique(matrix<T>& out, vector<int> &iRow)
{
  out.setup(0,0);

  // Setup vector for rows of final matrix that map to each row of initial matrix
  iRow.resize(this->dims[0]);
  iRow.assign(this->dims[0],-1);

  /* --- For each row in the matrix, compare to all
     previous rows to get first unique occurence --- */
  typename vector<T>::iterator itI, itJ;
  for (uint i=0; i<this->dims[0]; i++) {
    itI = this->data.begin() + i*this->dims[1];
    for (uint j=0; j<i; j++) {
      itJ = this->data.begin() + j*this->dims[1];
      if (equal(itI,itI+this->dims[1],itJ)) {
        iRow[i] = iRow[j];
        break;
      }
    }

    // If no prior occurance found, put in 'out' matrix
    if (iRow[i]==-1) {
      out.insertRow(&(this->data[i*this->dims[1]]),-1,this->dims[1]);
      iRow[i] = out.getDim0() - 1;
    }
  }
}

template<typename T, uint N>
T *Array<T,N>::getData(void)
{
  return data.data();
}

// Method to invert matrix - returns the inverse
template<typename T>
matrix<T> matrix<T>::invertMatrix(void)
{
  if (this->dims[0] != this->dims[1])
    FatalErrorST("Can only obtain inverse of a square matrix.");

  // Gaussian elimination with full pivoting
  // not to be used where speed is paramount

  int pivot_i = 0;
  int pivot_j = 0;
  int itemp_0;
  double mag;
  double max;
  double dtemp_0;
  double first;
  matrix <T> atemp_0(this->dims[0],1);
  matrix <T> identity(this->dims[0],this->dims[0]);
  matrix <T> input(this->dims[0],this->dims[0]);
  matrix <T> inverse(this->dims[0],this->dims[0]);
  matrix <T> inverse_out(this->dims[0],this->dims[0]);
  matrix <int> swap_0(this->dims[0],1);
  matrix <int> swap_1(this->dims[0],1);

  input = *this;

  // setup swap arrays
  for(uint i=0;i<this->dims[0];i++) {
    swap_0(i)=i;
    swap_1(i)=i;
  }

  // setup identity array
  for(uint i=0; i<this->dims[0]; i++) {
    for(uint j=0; j<this->dims[0]; j++) {
      identity(i,j)=0.0;
    }
    identity(i,i)=1.0;
  }

  // make triangular
  for(uint k=0; k<this->dims[0]-1; k++) {
    max=0;

    // find pivot
    for(uint i=k; i<this->dims[0]; i++) {
      for(uint j=k; j<this->dims[0]; j++) {
        mag=input(i,j)*input(i,j);
        if(mag>max) {
          pivot_i=i;
          pivot_j=j;
          max=mag;
        }
      }
    }

    // swap the swap arrays
    itemp_0=swap_0(k);
    swap_0(k)=swap_0(pivot_i);
    swap_0(pivot_i)=itemp_0;
    itemp_0=swap_1(k);
    swap_1(k)=swap_1(pivot_j);
    swap_1(pivot_j)=itemp_0;

    // swap the columns
    for(uint i=0; i<this->dims[0]; i++) {
      atemp_0(i)=input(i,pivot_j);
      input(i,pivot_j)=input(i,k);
      input(i,k)=atemp_0(i);
    }

    // swap the rows
    for(uint j=0; j<this->dims[0]; j++) {
      atemp_0(j)=input(pivot_i,j);
      input(pivot_i,j)=input(k,j);
      input(k,j)=atemp_0(j);
      atemp_0(j)=identity(pivot_i,j);
      identity(pivot_i,j)=identity(k,j);
      identity(k,j)=atemp_0(j);
    }

    // subtraction
    for(uint i=k+1; i<this->dims[0]; i++) {
      first=input(i,k);
      for(uint j=0; j<this->dims[0]; j++) {
        if(j>=k) {
          input(i,j)=input(i,j)-((first/input(k,k))*input(k,j));
        }
        identity(i,j)=identity(i,j)-((first/input(k,k))*identity(k,j));
      }
    }

    //exact zero
    for(uint j=0; j<k+1; j++) {
      for(uint i=j+1; i<this->dims[0]; i++) {
        input(i,j)=0.0;
      }
    }
  }

  // back substitute
  for(int i=(int)this->dims[0]-1; i>=0; i--) {
    for(uint j=0; j<this->dims[0]; j++) {
      dtemp_0=0.0;
      for(uint k=i+1; k<this->dims[0]; k++) {
        dtemp_0=dtemp_0+(input(i,k)*inverse(k,j));
      }
      inverse(i,j)=(identity(i,j)-dtemp_0)/input(i,i);
    }
  }

  // swap solution rows
  for(uint i=0; i<this->dims[0]; i++) {
    for(uint j=0; j<this->dims[0]; j++) {
      inverse_out(swap_1(i),j)=inverse(i,j);
    }
  }

  return inverse_out;
}

template<typename T>
double matrix<T>::det()
{
  if (this->dims[0]!=this->dims[1])
    FatalErrorST("Determinant only meaningful for square matrices.");

  if (this->dims[0] == 2) {
    // Base case
    return this->data[0]*this->data[3] - this->data[1]*this->data[2];
  }
  else {
    // Use minor-matrix recursion
    double Det = 0;
    int sign = -1;
    matrix<T> Minor(this->dims[0]-1,this->dims[1]-1);
    for (int row=0; row<this->dims[0]; row++) {
      sign *= -1;
      // Setup the minor matrix (expanding along first column)
      int i0 = 0;
      for (int i=0; i<this->dims[0]; i++) {
        if (i==row) continue;
        for (int j=1; j<this->dims[1]; j++) {
          Minor(i0,j-1) = this->operator()(i,j);
        }
        i0++;
      }
      // Add in the minor's determinant
      Det += sign*Minor.det()*this->operator()(row,0);
    }

    return Det;
  }
}

template<typename T>
matrix<T> matrix<T>::adjoint(void)
{
  if (this->dims[0]!=this->dims[1])
    FatalErrorST("Adjoint only meaningful for square matrices.");

  matrix<T> adj(this->dims[0],this->dims[1]);

  int signRow = -1;
  matrix<T> Minor(this->dims[0]-1,this->dims[1]-1);
  for (int row=0; row<this->dims[0]; row++) {
    signRow *= -1;
    int sign = -1*signRow;
    for (int col=0; col<this->dims[1]; col++) {
      sign *= -1;
      // Setup the minor matrix (expanding along row, col)
      int i0 = 0;
      for (int i=0; i<this->dims[0]; i++) {
        if (i==row) continue;
        int j0 = 0;
        for (int j=0; j<this->dims[1]; j++) {
          if (j==col) continue;
          Minor(i0,j0) = this->operator()(i,j);
          j0++;
        }
        i0++;
      }
      // Recall: adjoint is TRANSPOSE of cofactor matrix
      adj(col,row) = sign*Minor.det();
    }
  }

  return adj;
}


// Method to resize a std Vector to a matrix
template<typename T>
void matrix<T>::vecToMatrixResize(vector<T> &A)
{
  if(A.size() != this->dims[0]*this->dims[1]) {
    FatalErrorST("Cannot fill vector into matrix - must have vector size == this->dims[0]*this->dims[1].");
  }
  else {
    for(uint i=0; i<A.size(); i++)
      this->data[i] = A[i];
  }
}


// Fix for compiler to know which template types will be needed later (and therefore must now be compiled):
template class Array<int,1>;
template class Array<int,2>;
template class Array<int,3>;
template class Array<int,4>;

template class Array<double,1>;
template class Array<double,2>;
template class Array<double,3>;
template class Array<double,4>;

template class Array<double*,1>;
template class Array<double*,2>;
template class Array<double*,3>;
template class Array<double*,4>;

template class Array<point,1>;
template class Array<point,2>;
template class Array<point,3>;
template class Array<point,4>;

template class Array<set<int>,1>;
template class Array<set<int>,2>;
template class Array<set<int>,3>;
template class Array<set<int>,4>;

template class Array2D<int>;
template class Array2D<double>;
template class Array2D<point>;

template class matrix<int>;
template class matrix<double>;

template class Array<matrix<double>,2>;
