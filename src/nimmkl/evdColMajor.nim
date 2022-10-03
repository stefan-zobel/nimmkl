# Copyright (c) 2022 Stefan Zobel
# Distributed under the Apache v2 License (license terms are at https://www.apache.org/licenses/LICENSE-2.0).

import ../nimmkl,common
import std/complex

template ty(T: typedesc[float64]): typedesc[float64] = typedesc[float64]
template ty(T: typedesc[float32]): typedesc[float32] = typedesc[float32]
template ty(T: typedesc[Complex64]): typedesc[float64] = typedesc[float64]
template ty(T: typedesc[Complex32]): typedesc[float32] = typedesc[float32]

# generic ?geev

proc geev*[T: SomeElementType, U: SomeFloat](matrixLayout: TLayout; jobvl: TEigJob;
           jobvr: TEigJob; n: int; a: seq[T]; lda: int; w: var seq[Complex[U]];
           vl: seq[T]; ldvl: int; vr: seq[T]; ldvr: int): int =

  when T is float64:
    let wr = newSeq[float64](w.len)
    let wi = newSeq[float64](w.len)
    result = dgeev(matrixLayout, jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr)
    for i in 0..<wr.len:
      w[i] = complex64(wr[i], wi[i])
  when T is float32:
    let wr = newSeq[float32](w.len)
    let wi = newSeq[float32](w.len)
    result = sgeev(matrixLayout, jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr)
    for i in 0..<wr.len:
      w[i] = complex32(wr[i], wi[i])
  when T is Complex64:
    result = zgeev(matrixLayout, jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr)
  when T is Complex32:
    result = cgeev(matrixLayout, jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr)


proc newVlVrW[T](jobvr: TEigJob; n: int) : auto =
  let ld = max(1, n)
  let vl : seq[T] = @[]
  let w = newSeq[Complex[ty(T)]](n)
  var vr : seq[T]
  if jobvr == EigJobAll: vr = newSeq[T](n * n) else: vr = @[]
  (ld, vl, vr, w)

proc evdFullAuto[T: SomeElementType](n: int; a: seq[T]): auto =
  var (ld, vl, vr, w) = newVlVrW[T](EigJobAll, n)
  # use `var` instead of `let` to avoid `deepCopy(a)`
  var copy = a
  let r = geev[T, ty(T)](ColMajorOrder, ValuesOnly, EigJobAll, n, copy, ld, w, vl,
    ld, vr, ld)
  result = (r, vr, w)

proc evdAuto[T: SomeElementType](n: int; a: seq[T]): auto =
  var (ld, vl, vr, w) = newVlVrW[T](ValuesOnly, n)
  # use `var` instead of `let` to avoid `deepCopy(a)`
  var copy = a
  let r = geev[T, ty(T)](ColMajorOrder, ValuesOnly, ValuesOnly, n, copy, ld, w, vl,
    ld, vr, ld)
  result = (r, w)


# convenience eigenvalue decomposition

proc evd*[T: float64|Complex64](ty: typedesc[float64|Complex64], n: int; a: seq[T]): (int, seq[Complex64]) =
  evdAuto[T](n, a)

proc evd*[T: float32|Complex32](ty: typedesc[float32|Complex32], n: int; a: seq[T]): (int, seq[Complex32]) =
  evdAuto[T](n, a)

proc evdFull*[T: float64|Complex64](ty: typedesc[float64|Complex64], n: int; a: seq[T]): (int, seq[T], seq[Complex64]) =
  evdFullAuto[T](n, a)

proc evdFull*[T: float32|Complex32](ty: typedesc[float32|Complex32], n: int; a: seq[T]): (int, seq[T], seq[Complex32]) =
  evdFullAuto[T](n, a)

