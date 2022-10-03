# Copyright (c) 2022 Stefan Zobel
# Distributed under the Apache v2 License (license terms are at https://www.apache.org/licenses/LICENSE-2.0).

import ../nimmkl,common
import std/complex


# generic ?gels

proc gels*[T: SomeElementType](matrixLayout: TLayout; trans: TTrans; m: int;
           n: int; nrhs: int; a: seq[T]; lda: int; b: seq[T]; ldb: int): int =

  when T is float64:
    result = dgels(matrixLayout, trans, m, n, nrhs, a, lda, b, ldb)
  when T is float32:
    result = sgels(matrixLayout, trans, m, n, nrhs, a, lda, b, ldb)
  when T is Complex64:
    result = zgels(matrixLayout, trans, m, n, nrhs, a, lda, b, ldb)
  when T is Complex32:
    result = cgels(matrixLayout, trans, m, n, nrhs, a, lda, b, ldb)


# convenience qrsolve() implementation

proc qrsolve*[T: SomeElementType](m: int, n: int, rhsCount: int, a: seq[T],
              x: var seq[T], b: seq[T]): (int, seq[T]) =
  let lda = max(1, m)
  let ldb = max(1, max(m, n))
  var tmp : seq[T] = newSeq[T](ldb * rhsCount)

  for j in 0..<rhsCount:
    for i in 0..<m:
      tmp[j * ldb + i] = b[j * m + i]

  # use `var` instead of `let` to avoid `deepCopy(a)`
  var copy = a

  let r = gels[T](ColMajorOrder, NoTrans, m, n, rhsCount, copy, lda, tmp, ldb)

  for j in 0..<rhsCount:
    for i in 0..<n:
      x[j * n + i] = tmp[j * ldb + i]

  tmp.setLen(0)
  (r, x)


