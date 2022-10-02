# Copyright (c) 2022 Stefan Zobel
# Distributed under the Apache v2 License (license terms are at https://www.apache.org/licenses/LICENSE-2.0).

import ../nimmkl,common
import std/complex

proc ONE[T](): T =
  when T is float64:
    result = 1.0'f64
  when T is float32:
    result = 1.0'f32
  when T is Complex64:
    result = complex64(1.0'f64)
  when T is Complex32:
    result = complex32(1.0'f32)


# generic ?getrf

proc getrf*[T: SomeElementType](matrixLayout: TLayout; m: int; n: int; a: seq[T];
            lda: int; indices: seq[int32]): int =
  when T is float64:
    result = dgetrf(matrixLayout, m, n, a, lda, indices)
  when T is float32:
    result = sgetrf(matrixLayout, m, n, a, lda, indices)
  when T is Complex64:
    result = zgetrf(matrixLayout, m, n, a, lda, indices)
  when T is Complex32:
    result = cgetrf(matrixLayout, m, n, a, lda, indices)


proc copyIntoL[T](L: var seq[T], rowsL: int, colsL: int, A: seq[T], rowsA: int) =
  for col in 0..<colsL:
    for row in 0..<rowsL:
      if row == col:
        L[col * rowsL + row] = ONE[T]()
      elif row > col:
        L[col * rowsL + row] = A[col * rowsA + row]
      else:
        discard

proc copyIntoU[T](U: var seq[T], rowsU: int, colsU: int, A: seq[T], rowsA: int) =
  for col in 0..<colsU:
    for row in 0..<rowsU:
      if row <= col:
        U[col * rowsU + row] = A[col * rowsA + row]

proc getMaxRowIndex(pivot: seq[int32]): int =
  var maxRowIdx = -1
  for i in 0..<pivot.len:
    let row = pivot[i]
    if row > maxRowIdx:
      maxRowIdx = row
  result = maxRowIdx

proc initPermVector(pivot: seq[int32]): seq[int32] =
  var perm = newSeq[int32](max(pivot.len, getMaxRowIndex(pivot)))
  for i in 0..<perm.len:
    perm[i] = i.int32 + 1
  result = perm

proc swapRows(pivot: seq[int32], perm: var seq[int32]): bool =
  var swapped = false
  for j in 0..<pivot.len:
    let rowIndex = j + 1
    let newRow = pivot[j]
    if newRow != rowIndex:
      let oldRow = perm[rowIndex - 1]
      perm[rowIndex - 1] = perm[newRow - 1]
      perm[newRow - 1] = oldRow
      swapped = true

  if swapped:
    for i in 0..<perm.len:
      if perm[i] != i + 1:
        return true

  return false

proc checkPivotVector(pivot: seq[int32]): seq[int32] =
  if pivot.len > 0:
    var rowSwapDetected = false
    for i in 0..<pivot.len:
      if pivot[i] != i + 1:
        rowSwapDetected = true
        break
    if rowSwapDetected:
      var perm = initPermVector(pivot)
      let positionsChanged = swapRows(pivot, perm)
      if positionsChanged:
        return perm
  return @[]

proc buildMatrix[T](perm: seq[int32], dim: int): seq[T] =
  let n = perm.len
  var permMatrix = newSeq[T](dim * dim)
  for i in 0..<dim:
    if i < n:
      permMatrix[i * dim + (perm[i] - 1)] = ONE[T]()
    else:
      permMatrix[i * dim + i] = ONE[T]()
  result = permMatrix

proc genEye[T](dim: int): seq[T] =
  var identity = newSeq[T](dim * dim)
  for i in 0..<dim:
    identity[i * dim + i] = ONE[T]()
  result = identity

proc genPermutationMatrix[T](partialPivot: seq[int32], dim: int): seq[T] =
  let permutation = checkPivotVector(partialPivot)
  if permutation.len > 0:
    return buildMatrix[T](permutation, dim)
  return genEye[T](dim)


# convenience LU decomposition

proc lu*[T: SomeElementType](m: int; n: int; a: seq[T]): (int, seq[T], seq[T], seq[T], seq[int32]) =
  let lda = max(1, m)
  let pivot : seq[int32] = newSeq[int32](min(m, n))
  var P : seq[T]
  var L : seq[T]
  var U : seq[T]
  var colsL : int
  var rowsU : int
  if m >= n:
    L = newSeq[T](m * n)
    U = newSeq[T](n * n)
    colsL = n
    rowsU = n
  else:
    L = newSeq[T](m * m)
    U = newSeq[T](m * n)
    colsL = m
    rowsU = m
  let copy = a

  let r = getrf[T](ColMajorOrder, m, n, copy, lda, pivot)
  if r >= 0: # we also compute if U is exactly singular
    copyIntoL(L, m, colsL, copy, m)
    copyIntoU(U, rowsU, n, copy, m)
    P = genPermutationMatrix[T](pivot, m)
  (r, P, L, U, pivot)

