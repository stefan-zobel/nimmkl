# Copyright (c) 2022 Stefan Zobel
# Distributed under the Apache v2 License (license terms are at https://www.apache.org/licenses/LICENSE-2.0).

import ../nimmkl,common
import std/complex


# generic ?gesv

proc gesv*[T: SomeElementType](matrixLayout: TLayout; n: int; nrhs: int; a: seq[T];
           lda: int; indices: seq[int32]; b: seq[T]; ldb: int): int =

  when T is float64:
    result = dgesv(matrixLayout, n, nrhs, a, lda, indices, b, ldb)
  when T is float32:
    result = sgesv(matrixLayout, n, nrhs, a, lda, indices, b, ldb)
  when T is Complex64:
    result = zgesv(matrixLayout, n, nrhs, a, lda, indices, b, ldb)
  when T is Complex32:
    result = cgesv(matrixLayout, n, nrhs, a, lda, indices, b, ldb)


proc SIZE[T](): Natural =
  when T is float64:
    result = 8
  when T is float32:
    result = 4
  when T is Complex64:
    result = 16
  when T is Complex32:
    result = 8


# convenience lusolve() implementation

proc lusolve*[T: SomeElementType](n: int; rhsCount: int, a: seq[T], x: seq[T],
              b: seq[T]): (int, seq[T]) =
  let lda = max(1, n)
  let ldb = lda
  let xSize = n * rhsCount * SIZE[T]()
  let ipiv: seq[int32] = newSeq[int32](n)
  moveMem(unsafeAddr(x[0]), unsafeAddr(b[0]), xSize)
  # use `var` instead of `let` to avoid `deepCopy(a)`
  var copy = a
  let r = gesv[T](ColMajorOrder, n, rhsCount, copy, lda, ipiv, x, ldb)
  (r, x)

