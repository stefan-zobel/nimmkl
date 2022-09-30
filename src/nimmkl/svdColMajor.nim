import ../nimmkl,common
import std/complex


template ty(T: typedesc[float64]): typedesc[float64] = typedesc[float64]
template ty(T: typedesc[float32]): typedesc[float32] = typedesc[float32]
template ty(T: typedesc[Complex64]): typedesc[float64] = typedesc[float64]
template ty(T: typedesc[Complex32]): typedesc[float32] = typedesc[float32]

proc newUSVt[T, U](m: int, n: int, jobz: TSvdJob) : (seq[T], seq[U], seq[T]) =
  var u : seq[T]
  var s : seq[U]
  var vt: seq[T]

  case jobz
  of SvdJobAll:
    u = newSeq[T](m * m)
    s = newSeq[U](min(m, n))
    vt = newSeq[T](n * n)
  of None:
    s = newSeq[U](min(m, n))
  of Part:
    u = newSeq[T](m * min(m, n))
    s = newSeq[U](min(m, n))
    vt = newSeq[T](min(m, n) * n)
  else:
    discard

  (u, s, vt)


proc gesdd*[T: SomeElementType, U: SomeFloat](matrixLayout: TLayout;
  jobz: TSvdJob; m: int; n: int; a: seq[T]; lda: int; s: seq[U]; u: seq[T];
  ldu: int; vt: seq[T]; ldvt: int): int =

  when T is float64:
    result = dgesdd(matrixLayout, jobz, m, n, a, lda, s, u, ldu, vt, ldvt)
  when T is float32:
    result = sgesdd(matrixLayout, jobz, m, n, a, lda, s, u, ldu, vt, ldvt)
  when T is Complex64:
    result = zgesdd(matrixLayout, jobz, m, n, a, lda, s, u, ldu, vt, ldvt)
  when T is Complex32:
    result = cgesdd(matrixLayout, jobz, m, n, a, lda, s, u, ldu, vt, ldvt)


proc svdAllAuto[T: SomeElementType](m: int; n: int; a: seq[T]): auto =
  let (u, s, vt) = newUSVt[T, ty(T)](m, n, TSvdJob.SvdJobAll)
  let r = gesdd[T, ty(T)](TLayout.ColMajorOrder, TSvdJob.SvdJobAll, m, n, a, max(1, m),
    s, u, max(1, m), vt, max(1, n))
  (r, u, s, vt)

proc svdNoneAuto[T: SomeElementType](m: int; n: int; a: seq[T]): auto =
  let (u, s, vt) = newUSVt[T, ty(T)](m, n, TSvdJob.None)
  let r = gesdd[T, ty(T)](TLayout.ColMajorOrder, TSvdJob.None, m, n, a, max(1, m), s, u,
    1, vt, 1)
  (r, u, s, vt)

proc svdEconAuto[T: SomeElementType](m: int; n: int; a: seq[T]): auto =
  let (u, s, vt) = newUSVt[T, ty(T)](m, n, TSvdJob.Part)
  let r = gesdd[T, ty(T)](TLayout.ColMajorOrder, TSvdJob.Part, m, n, a, max(1, m), s, u,
    max(1, m), vt, max(1, min(m, n)))
  (r, u, s, vt)


proc svdAll*[T: float64|Complex64](m: int; n: int; a: seq[T]): (int, seq[T], seq[float64], seq[T]) =
  svdAllAuto[T](m, n, a)

proc svdAll*[T: float32|Complex32](m: int; n: int; a: seq[T]): (int, seq[T], seq[float32], seq[T]) =
  svdAllAuto[T](m, n, a)

proc svdNone*[T: float64|Complex64](m: int; n: int; a: seq[T]): (int, seq[T], seq[float64], seq[T]) =
  svdNoneAuto[T](m, n, a)

proc svdNone*[T: float32|Complex32](m: int; n: int; a: seq[T]): (int, seq[T], seq[float32], seq[T]) =
  svdNoneAuto[T](m, n, a)

proc svdEcon*[T: float64|Complex64](m: int; n: int; a: seq[T]): (int, seq[T], seq[float64], seq[T]) =
  svdEconAuto[T](m, n, a)

proc svdEcon*[T: float32|Complex32](m: int; n: int; a: seq[T]): (int, seq[T], seq[float32], seq[T]) =
  svdEconAuto[T](m, n, a)

