import ../nimmkl,common,mklTypes
import std/complex

# generic ?gemm

proc gemm*[T: SomeElementType](layout: Cblas_Layout; transA: Cblas_Transpose;
           transB: Cblas_Transpose; m: int; n: int; k: int; alpha: T; a: seq[T];
           lda: int; b: seq[T]; ldb: int; beta: T; c: seq[T]; ldc: int) =

  when T is float64:
    dgemm(layout, transA, transB, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
  when T is float32:
    sgemm(layout, transA, transB, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
  when T is Complex64:
    zgemm3m(layout, transA, transB, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
  when T is Complex32:
    cgemm3m(layout, transA, transB, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)


# convenience multiplication A * B

proc mul*[T: SomeElementType](m: int; n: int; k: int; a: seq[T]; b: seq[T]): seq[T] =
  when T is float64:
    let alpha = 1'f64
    let beta = alpha
  when T is float32:
    let alpha = 1'f32
    let beta = alpha
  when T is Complex64:
    let alpha = complex64(1'f64)
    let beta = alpha
  when T is Complex32:
    let alpha = complex32(1'f32)
    let beta = alpha

  result = newseq[T](m * n)
  gemm[T](CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, a, max(1, n),
    b, max(1, k), beta, result, max(1, m))

