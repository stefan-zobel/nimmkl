import ../nimmkl,common
import std/complex


proc geqrf*[T: SomeElementType](matrixLayout: TLayout; m: int; n: int; a: seq[T];
            lda: int; tau: seq[T]): int =

  when T is float64:
    result = dgeqrf(matrixLayout, m, n, a, lda, tau)
  when T is float32:
    result = sgeqrf(matrixLayout, m, n, a, lda, tau)
  when T is Complex64:
    result = zgeqrf(matrixLayout, m, n, a, lda, tau)
  when T is Complex32:
    result = cgeqrf(matrixLayout, m, n, a, lda, tau)


proc orgqr*[T: SomeElementType](matrixLayout: TLayout; m: int; n: int; k: int;
            a: seq[T]; lda: int; tau: seq[T]): int =

  when T is float64:
    result = dorgqr(matrixLayout, m, n, k, a, lda, tau)
  when T is float32:
    result = sorgqr(matrixLayout, m, n, k, a, lda, tau)
  when T is Complex64:
    result = zungqr(matrixLayout, m, n, k, a, lda, tau)
  when T is Complex32:
    result = cungqr(matrixLayout, m, n, k, a, lda, tau)


proc setUpperTrapezoidal[T: SomeElementType](rowsDest: int, colsDest: int,
                         dest: var seq[T], src: seq[T]) =
  for col in 0..<colsDest:
    for row in 0..<rowsDest:
      if row <= col:
        let idx = col * rowsDest + row
        dest[idx] = src[idx]


proc qr*[T: SomeElementType](m: int; n: int; a: seq[T]) : (int, seq[T], seq[T]) =
  let k = min(m, n)
  let lda = max(1, m)
  let tau = newSeq[T](k)
  var R = newSeq[T](n * n)
  let Q = a
  var rc = geqrf[T](ColMajorOrder, m, n, Q, lda, tau)

  if rc == 0:
    setUpperTrapezoidal[T](n, n, R, Q)
    rc = orgqr[T](ColMajorOrder, m, n, k, Q, lda, tau)

  (rc, Q, R)

