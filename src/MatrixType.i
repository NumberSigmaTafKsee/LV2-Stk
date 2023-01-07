////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 1996-2022 The Octave Project Developers
//
// See the file COPYRIGHT.md in the top-level directory of this
// distribution or <https://octave.org/copyright/>.
//
// This file is part of Octave.
//
// Octave is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Octave is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Octave; see the file COPYING.  If not, see
// <https://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////

%{
#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>
#include <octave/interpreter.h>
%}

class MatrixType
{
public:
  
  enum matrix_type
  {
    Unknown = 0,
    Full,
    Diagonal,
    Permuted_Diagonal,
    Upper,
    Lower,
    Permuted_Upper,
    Permuted_Lower,
    Banded,
    Hermitian,
    Banded_Hermitian,
    Tridiagonal,
    Tridiagonal_Hermitian,
    Rectangular
  };

  MatrixType (void);

  MatrixType (const MatrixType& a);

  MatrixType (const Matrix& a);

  MatrixType (const ComplexMatrix& a);

  MatrixType (const FloatMatrix& a);

  MatrixType (const FloatComplexMatrix& a);

  template <typename T>
  OCTAVE_API
  MatrixType (const MSparse<T>& a);

  MatrixType (const matrix_type t, bool _full = false);

  MatrixType (const matrix_type t, const octave_idx_type np,
                         const octave_idx_type *p, bool _full = false);

  MatrixType (const matrix_type t, const octave_idx_type ku,
                         const octave_idx_type kl, bool _full = false);

  ~MatrixType (void);

  MatrixType& operator = (const MatrixType& a);

  int type (bool quiet = true);

  int type (const Matrix& a);

  int type (const ComplexMatrix& a);

  int type (const FloatMatrix& a);

  int type (const FloatComplexMatrix& a);

  int type (const SparseMatrix& a);

  int type (const SparseComplexMatrix& a);

  double band_density (void) const;

  int nupper (void) const;

  int nlower (void) const;

  bool is_dense (void) const;

  bool isdiag (void) const;
  bool istriu (void) const;
  bool istril (void) const;
  bool isbanded (void) const;
  bool is_tridiagonal (void) const;
  bool ishermitian (void) const;

  bool is_rectangular (void) const;

  bool is_known (void) const;

  bool is_unknown (void) const;

  void info (void) const;

  octave_idx_type * triangular_perm (void) const;

  void invalidate_type (void);

  void mark_as_diagonal (void);

  void mark_as_permuted_diagonal (void);

  void mark_as_upper_triangular (void);

  void mark_as_lower_triangular (void);

  void mark_as_tridiagonal (void);

  void mark_as_banded (const octave_idx_type ku, const octave_idx_type kl);
  void mark_as_full (void);
  void mark_as_rectangular (void);
  void mark_as_dense (void);
  void mark_as_not_dense (void);
  void mark_as_symmetric (void);
  void mark_as_unsymmetric (void);
  void mark_as_permuted (const octave_idx_type np, const octave_idx_type *p);
  void mark_as_unpermuted (void);
  MatrixType transpose (void) const;
};
