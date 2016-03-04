/*!
 * \file matrix.cpp
 * \brief Class for simplified matrix storage & manipulation
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 1.0.0
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2015 Jacob Crabill
 *
 * Flurry++ is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * Flurry++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Flurry++; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA..
 *
 */
#include "../include/matrix.hpp"

#include <set>

#ifndef _NO_MPI
#include "mpi.h"
#endif

template<typename T, uint N>
Array<T,N>::Array()
{
  data.resize(0);
  dims[0] = 0;
  dims[1] = 0;
  dims[2] = 0;
  dims[3] = 0;
}

template<typename T, uint N>
Array<T,N>::Array(uint inNDim0, uint inNDim1, uint inDim2, uint inDim3)
{
  data.resize(inNDim0*inNDim1*inDim2*inDim3);
  dims[0] = inNDim0;
  dims[1] = inNDim1;
  dims[2] = inDim2;
  dims[3] = inDim3;
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

  if (this->data.size()>0) {
    this->dims[2] = 1;
    this->dims[3] = 1;
  }
}

template<typename T>
matrix<T>::matrix(const Array2D<T> &inMatrix)
{
  this->data = inMatrix.data;
  for (uint i = 0; i < 4; i++)
    this->dims[i] = inMatrix.dims[i];
}

template<typename T,uint N>
Array<T,N>::Array(const Array<T,N> &inMatrix)
{
  data = inMatrix.data;
  for (uint i = 0; i < 4; i++)
    dims[i] = inMatrix.dims[i];
}

template<typename T, uint N>
Array<T,N> Array<T,N>::operator=(const Array<T,N> &inMatrix)
{
  data = inMatrix.data;
  for (uint i = 0; i < 4; i++)
    dims[i] = inMatrix.dims[i];
  return *this;
}

template<typename T, uint N>
void Array<T,N>::setup(uint inDim0, uint inDim1, uint inDim2, uint inDim3)
{
  data.resize(inDim0*inDim1*inDim2*inDim3);
  dims[0] = inDim0;
  dims[1] = inDim1;
  dims[2] = inDim2;
  dims[3] = inDim3;
}

template <typename T, uint N>
void Array<T,N>::add_dim_0(uint ind, const T &val)
{
  /* Insert new 'book' of memory */
  uint bookSize = dims[3]*dims[2]*dims[1];
  auto it = data.begin() + bookSize*ind;
  data.insert(it, bookSize, val);
  dims[0]++;
}

template <typename T, uint N>
void Array<T,N>::add_dim_1(uint ind, const T &val)
{
  /* Insert new 'page' of memory */
  uint stride = dims[3]*dims[2]*dims[1];
  uint pageSize = dims[2]*dims[3];
  uint offset = ind*pageSize;
  for (int i = dims[0]-1; i >= 0; i--) {
    auto it = data.begin() + i*stride + offset;
    data.insert(it, pageSize, val);
  }

  dims[1]++;
}

template <typename T, uint N>
void Array<T,N>::add_dim_2(uint ind, const T &val)
{
  /* Insert new 'column' of memory */
  uint stride0 = dims[3]*dims[2]*dims[1];
  uint stride1 = dims[3]*dims[2];
  uint colSize = dims[3];
  uint offset = ind*colSize;
  for (int i = dims[0]-1; i >= 0; i--) {
    for (int j = dims[1]-1; j >= 0; j--) {
      auto it = data.begin() + i*stride0 + j*stride1 + offset;
      data.insert(it, colSize, val);
    }
  }

  dims[2]++;
}

template <typename T, uint N>
void Array<T,N>::add_dim_3(uint ind, const T &val)
{
  /* Insert new 'row' of memory */
  uint stride0 = dims[3]*dims[2]*dims[1];
  uint stride1 = dims[3]*dims[2];
  uint stride2 = dims[3];
  uint rowSize = 1;
  uint offset = ind*rowSize;
  for (int i = dims[0]-1; i >= 0; i--) {
    for (int j = dims[1]-1; j >= 0; j--) {
      for (int k = dims[2]-1; k >= 0; k--) {
        auto it = data.begin() + i*stride0 + j*stride1 + k*stride2 + offset;
        data.insert(it, rowSize, val);
      }
    }
  }

  dims[3]++;
}

template <typename T, uint N>
void Array<T,N>::remove_dim_0(uint ind)
{
  /* Insert new 'book' of memory */
  uint bookSize = dims[3]*dims[2]*dims[1];
  auto it = data.begin() + bookSize*ind;
  data.erase(it, it+bookSize);
  dims[0]--;
}

template <typename T, uint N>
void Array<T,N>::remove_dim_1(uint ind)
{
  /* Insert new 'page' of memory */
  uint stride = dims[3]*dims[2]*dims[1];
  uint pageSize = dims[2]*dims[3];
  uint offset = ind*pageSize;
  for (int i = dims[0]-1; i >= 0; i--) {
    auto it = data.begin() + i*stride + offset;
    data.erase(it, it+pageSize);
  }

  dims[1]--;
}

template <typename T, uint N>
void Array<T,N>::remove_dim_2(uint ind)
{
  /* Insert new 'column' of memory */
  uint stride0 = dims[3]*dims[2]*dims[1];
  uint stride1 = dims[3]*dims[2];
  uint colSize = dims[3];
  uint offset = ind*colSize;
  for (int i = dims[0]-1; i >= 0; i--) {
    for (int j = dims[1]-1; j >= 0; j--) {
      auto it = data.begin() + i*stride0 + j*stride1 + offset;
      data.erase(it, it+colSize);
    }
  }

  dims[2]--;
}

template <typename T, uint N>
void Array<T,N>::remove_dim_3(uint ind)
{
  /* Insert new 'row' of memory */
  uint stride0 = dims[3]*dims[2]*dims[1];
  uint stride1 = dims[3]*dims[2];
  uint stride2 = dims[3];
  uint rowSize = 1;
  uint offset = ind*rowSize;
  for (int i = dims[0]-1; i >= 0; i--) {
    for (int j = dims[1]-1; j >= 0; j--) {
      for (int k = dims[2]-1; k >= 0; k--) {
        auto it = data.begin() + i*stride0 + j*stride1 + k*stride2 + offset;
        data.erase(it, it+rowSize);
      }
    }
  }

  dims[3]--;
}

template<typename T>
void matrix<T>::addMatrix(matrix<T> &A, double a)
{
#ifdef _DEBUG
  if (A.dims[0] != this->dims[0] || A.dims[1] != this->dims[1])
    FatalErrorST("Incompatible matrix sizes for addMatrix.");
#endif
  for (uint i=0; i<this->dims[0]; i++)
    for (uint j=0; j<this->dims[1]; j++)
      this->data[i*this->dims[1]+j] += a*A(i,j);
}

template<>
matrix<double>& matrix<double>::operator+=(matrix<double> &A)
{
#ifdef _DEBUG
  if (A.dims[0] != this->dims[0] || A.dims[1] != this->dims[1])
    FatalErrorST("Incompatible matrix sizes for addMatrix.");
#endif
  for (uint i=0; i<this->dims[0]; i++)
    for (uint j=0; j<this->dims[1]; j++)
      this->data[i*this->dims[1]+j] += A(i,j);

  return *this;
}

template<>
matrix<double>& matrix<double>::operator-=(matrix<double> &A)
{
#ifdef _DEBUG
  if (A.dims[0] != this->dims[0] || A.dims[1] != this->dims[1])
    FatalErrorST("Incompatible matrix sizes for addMatrix.");
#endif
  for (uint i=0; i<this->dims[0]; i++)
    for (uint j=0; j<this->dims[1]; j++)
      this->data[i*this->dims[1]+j] -= A(i,j);

  return *this;
}

template<>
matrix<double>& matrix<double>::operator*=(double a)
{
  for (auto &val:this->data) val *= a;
  return *this;
}

template<>
matrix<double>& matrix<double>::operator/=(double a)
{
  for (auto &val:this->data) val /= a;
  return *this;
}

template<typename T, uint N>
T* Array<T,N>::operator[](int inRow)
{
#ifdef _DEBUG
  if (inRow >= (int)this->dims[0] || inRow < 0)
  {
    cout << "inRow = " << inRow << ", nRows = " << this->dims[0] << endl;
    FatalErrorST("Operator[]: Attempted out-of-bounds access in matrix.");
  }
#endif
  return &data[inRow*this->dims[1]];
}

template<typename T, uint N>
T& Array<T,N>::operator()(int i, int j, int k, int l)
{
#ifdef _DEBUG
  if (i>=(int)this->dims[0] || i<0 || j>=(int)this->dims[1] || j<0 ||
      k>=(int)this->dims[2] || k<0 || l>=(int)this->dims[3] || l<0)
  {
    cout << "i=" << i << ", dim0=" << dims[0] << ", j=" << j << ", dim1=" << dims[1] << ", ";
    cout << "k=" << k << ", dim2=" << dims[2] << ", l=" << l << ", dim3=" << dims[3] << endl;
    FatalErrorST("Attempted out-of-bounds access in Array.");
  }
#endif
  return data[l+dims[3]*(k+dims[2]*(j+dims[1]*i))];
}

template<typename T, uint N>
T Array<T,N>::operator()(int i, int j, int k, int l) const
{
#ifdef _DEBUG
  if (i>=(int)this->dims[0] || i<0 || j>=(int)this->dims[1] || j<0 ||
      k>=(int)this->dims[2] || k<0 || l>=(int)this->dims[3] || l<0)
  {
    cout << "i=" << i << ", dim0=" << dims[0] << ", j=" << j << ", dim1=" << dims[1] << ", ";
    cout << "k=" << k << ", dim2=" << dims[2] << ", l=" << l << ", dim3=" << dims[3] << endl;
    FatalErrorST("Attempted out-of-bounds access in Array.");
  }
#endif
  return data[l+dims[3]*(k+dims[2]*(j+dims[1]*i))];
}

template<typename T>
T& Array2D<T>::operator()(int i, int j)
{
#ifdef _DEBUG
  if (i>=(int)this->dims[0] || i<0 || j>=(int)this->dims[1] || j<0)
  {
    cout << "i=" << i << ", dim0=" << this->dims[0] << ", j=" << j << ", dim1=" << this->dims[1] << endl;
    FatalErrorST("Attempted out-of-bounds access in matrix.");
  }
#endif
  return this->data[j+this->dims[1]*i];
}

template<typename T>
T Array2D<T>::operator()(int i, int j) const
{
#ifdef _DEBUG
  if (i>=(int)this->dims[0] || i<0 || j>=(int)this->dims[1] || j<0)
  {
    cout << "i=" << i << ", dim0=" << this->dims[0] << ", j=" << j << ", dim1=" << this->dims[1] << endl;
    FatalErrorST("Attempted out-of-bounds access in matrix.");
  }
#endif
  return this->data[j+this->dims[1]*i];
}

template<typename T>
void matrix<T>::initializeToZero(void)
{
  for (auto &val:this->data) val = 0;
}

template<typename T, uint N>
void Array<T,N>::initializeToValue(const T &_val)
{
  for (auto &val:this->data) val = _val;
}

template <typename T>
void matrix<T>::timesMatrix(matrix<T> &A, matrix<T> &B)
{
  uint i, j, k, p;
  p = A.dims[1];
#ifdef _DEBUG
  if (A.dims[0] != this->dims[1]) FatalErrorST("Incompatible matrix sizes in matrix multiplication!");
#endif
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

template <typename T>
void matrix<T>::timesMatrixPlus(matrix<T> &A, matrix<T> &B)
{
  uint p = A.dims[1];
#ifdef _DEBUG
  if (A.dims[0] != this->dims[1]) FatalErrorST("Incompatible matrix sizes in matrix multiplication!");
#endif
  if (B.dims[0] != this->dims[0] || B.dims[1] != A.dims[1]) B.setup(this->dims[0], A.dims[1]);

  for (uint i=0; i<this->dims[0]; i++) {
    for (uint j=0; j<this->dims[1]; j++) {
      for (uint k=0; k<p; k++) {
        B[i][k] += this->data[i*this->dims[1]+j]*A[j][k];
      }
    }
  }
}

template <typename T>
void matrix<T>::timesVector(vector<T> &A, vector<T> &B)
{
#ifdef _DEBUG
  if (A.size() != this->dims[1]) FatalError("Incompatible vector size");
#endif
  if (B.size() != this->dims[1]) B.resize(this->dims[1]);

  for (uint i=0; i<this->dims[0]; i++) {
    B[i] = 0;
    for (uint j=0; j<this->dims[1]; j++) {
      B[i] += this->data[i*this->dims[1]+j]*A[j];
    }
  }
}

template<typename T>
void Array2D<T>::appendRows(Array2D<T> &mat)
{
#ifdef _DEBUG
  if (this->dims[1]!= 0 && mat.getDim1()!=this->dims[1])
    FatalErrorST("Attempting to append rows of wrong size to matrix.");
#endif
  this->data.insert(this->data.end(),mat.getData(),mat.getData()+mat.getSize());

  if (this->dims[1]==0) this->dims[1]=mat.getDim1();
  this->dims[0]+=mat.getDim0();
}

template<typename T>
void Array2D<T>::insertRow(const vector<T> &vec, int rowNum)
{
#ifdef _DEBUG
  if (this->dims[1]!= 0 && vec.size()!=this->dims[1])
    FatalErrorST("Attempting to assign row of wrong size to matrix.");
#endif
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
#ifdef _DEBUG
  if (this->dims[1]!=0 && length!=(int)this->dims[1])
    FatalErrorST("Attempting to assign row of wrong size to matrix.");
#endif
  if (rowNum==INSERT_AT_END || rowNum==(int)this->dims[0]) {
    // Default action - add to end
    this->data.insert(this->data.end(),vec,vec+length);
  }else{
    // Insert at specified location
    this->data.insert(this->data.begin()+rowNum*this->dims[1],vec,vec+length);
  }

  if (this->dims[0]==0) {
    this->dims[1] = (uint)length;
    this->dims[2] = 1;
    this->dims[3] = 1;
  }

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

  if (this->dims[0]==0) {
    this->dims[1] = vec.size();
    this->dims[2] = 1;
    this->dims[3] = 1;
  }

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

  if (this->dims[0]==0) {
    this->dims[1] = (uint)length;
    this->dims[2] = 1;
    this->dims[3] = 1;
  }

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
#ifdef _DEBUG
  if (row > this->dims[0]) FatalErrorST("Attempting to grab row beyond end of matrix.");
#endif
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
Array2D<T> Array2D<T>::transpose(void)
{
  Array2D<T> out(this->dims[1],this->dims[0]);
  for (int i=0; i<this->dims[0]; i++) {
    for (int j=0; j<this->dims[1]; j++) {
      out(j,i) = this->operator()(i,j);
    }
  }
  return out;
}

template<typename T>
Array2D<T> Array2D<T>::slice(array<int,2> rows, array<int,2> cols)
{
  Array2D<T> out(rows[1]-rows[0]+1,cols[1]-cols[0]+1);
  int ii=0;
  for (int i=rows[0]; i<=rows[1]; i++) {
    int jj = 0;
    for (int j=cols[0]; j<=cols[1]; j++) {
      out(ii,jj) = this->operator()(i,j);
      jj++;
    }
    ii++;
  }
  return out;
}

template<typename T>
void matrix<T>::print(int prec) const
{
  if (this->data.size()>0) {
    cout << endl;
    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(prec);
    for (uint i=0; i<this->dims[0]; i++) {
      for (uint j=0; j<this->dims[1]; j++) {
        cout << setw(prec+4) << left << this->data[i*this->dims[1]+j] << " ";
      }
      cout << endl;
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

template<typename T>
matrix<T> matrix<T>::invertMatrix(void)
{
#ifdef _DEBUG
  if (this->dims[0] != this->dims[1])
    FatalErrorST("Can only obtain inverse of a square matrix.");
#endif
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
vector<T> matrix<T>::solve(vector<T> RHS)
{
#ifdef _DEBUG
  if (this->dims[0] != this->dims[1])
    FatalErrorST("Can only use Gaussian elimination on a square matrix.");
#endif
  // Gaussian elimination with full pivoting
  // not to be used where speed is paramount

  vector<T> vec(this->dims[0]);
  vector<T> solution(this->dims[0]);
  vector<int> swap_0(this->dims[0],1);
  vector<int> swap_1(this->dims[0],1);
  vector<T> tmpRow(this->dims[0]);

  matrix<T> LHS(*this);

  // setup swap arrays
  for(uint i=0; i<this->dims[0]; i++) {
    swap_0[i]=i;
    swap_1[i]=i;
  }

  // make triangular
  for(uint k=0; k<this->dims[0]-1; k++) {
    T max = 0;
    T mag;

    // find pivot
    int pivot_i = 0;
    int pivot_j = 0;
    for(uint i=k; i<this->dims[0]; i++) {
      for(uint j=k; j<this->dims[0]; j++) {
        mag = LHS(i,j)*LHS(i,j);
        if(mag>max) {
          pivot_i = i;
          pivot_j = j;
          max = mag;
        }
      }
    }

    // swap the swap arrays
    int itemp_0 = swap_0[k];
    swap_0[k] = swap_0[pivot_i];
    swap_0[pivot_i] = itemp_0;
    itemp_0 = swap_1[k];
    swap_1[k] = swap_1[pivot_j];
    swap_1[pivot_j] = itemp_0;

    // swap the columns
    for(uint i=0; i<this->dims[0]; i++) {
      tmpRow[i] = LHS(i,pivot_j);
      LHS(i,pivot_j) = LHS(i,k);
      LHS(i,k) = tmpRow[i];
    }

    // swap the rows
    for(uint j=0; j<this->dims[0]; j++) {
      tmpRow[j] = LHS(pivot_i,j);
      LHS(pivot_i,j) = LHS(k,j);
      LHS(k,j) = tmpRow[j];
    }
    T tmp = RHS[pivot_i];
    RHS[pivot_i] = RHS[k];
    RHS[k] = tmp;

    // subtraction
    for(uint i=k+1; i<this->dims[0]; i++) {
      T first = LHS(i,k);
      RHS[i] = RHS[i] - first/LHS(k,k)*RHS[k];
      for(uint j=k; j<this->dims[0]; j++)
        LHS(i,j) = LHS(i,j) - (first/LHS(k,k))*LHS(k,j);
    }

    // exact zero
    for(uint j=0; j<k+1; j++) {
      for(uint i=j+1; i<this->dims[0]; i++) {
        LHS(i,j)=0.0;
      }
    }
  }

  // back substitute
  for(int i=(int)this->dims[0]-1; i>=0; i--) {
    T dtemp_0 = 0.0;
    for(uint k = i+1; k<this->dims[0]; k++) {
      dtemp_0 = dtemp_0 + (LHS(i,k)*vec[k]);
    }
    vec[i] = (RHS[i]-dtemp_0)/LHS(i,i);
  }

  // swap solution rows
  for(uint i=0; i<this->dims[0]; i++)
    solution[swap_1[i]] = vec[i];

  return solution;
}

template<typename T>
double matrix<T>::det()
{
#ifdef _DEBUG
  if (this->dims[0]!=this->dims[1])
    FatalErrorST("Determinant only meaningful for square matrices.");
#endif
  if (this->dims[0] == 1) {
    return this->data[0];
  }
  else if (this->dims[0] == 2) {
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
#ifdef _DEBUG
  if (this->dims[0]!=this->dims[1])
    FatalErrorST("Adjoint only meaningful for square matrices.");
#endif
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

template<typename T>
T matrix<T>::frobNorm(void)
{
  T norm = 0;
  for (auto &val:this->data) norm += val*val;
  return norm;
}

// Method to resize a std Vector to a matrix
template<typename T>
void matrix<T>::vecToMatrixResize(vector<T> &A)
{
#ifdef _DEBUG
  if(A.size() != this->dims[0]*this->dims[1])
    FatalErrorST("Cannot fill vector into matrix - must have vector size == this->dims[0]*this->dims[1].");
#endif
  for(uint i=0; i<A.size(); i++)
    this->data[i] = A[i];
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
template class Array2D<Array<double,3>>;
template class Array<Array<double,3>,2>;

template class matrix<int>;
template class matrix<double>;

template class Array<matrix<double>,2>;
