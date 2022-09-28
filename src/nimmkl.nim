# Copyright (c) 2022 Stefan Zobel
# Distributed under the Apache v2 License (license terms are at https://www.apache.org/licenses/LICENSE-2.0).

import nimmkl/[mklTypes, mklTrans, mklCblas, mklLapacke, mklService]
import std/complex

type
  TLayout* = enum
    RowMajorOrder = LAPACK_ROW_MAJOR,
    ColMajorOrder = LAPACK_COL_MAJOR

type
  TDiag* = enum
    NonUnit = 'N',
    Unit = 'U'

  TEigJob* = enum
    ValuesOnly = 'N',
    EigJobAll = 'V'

  TNorm* = enum
    One = '1',
    Inf = 'I'

  TRange* = enum
    RangeAll = 'A',
    Indices = 'I',
    Interval = 'V'

  TSide* = enum
    Left = 'L',
    Right = 'R'

  TSvdJob* = enum
    SvdJobAll = 'A',
    None = 'N',
    Overwrite = 'O',
    Part = 'S'

  TTrans* = enum
    ConjTrans = 'C',
    NoTrans = 'N',
    Conj = 'R',
    Trans = 'T'

  TUpLo* = enum
    Lower = 'L',
    Upper = 'U'

  TOrder* = enum
    ColMajor = 'C',
    RowMajor = 'R'


proc toPtr(a: seq[int]): ptr cint =
  cast[ptr cint](unsafeAddr(a[0]))

proc toPtr(a: seq[float64]): ptr cdouble =
  unsafeAddr(a[0])

proc toPtr(a: seq[float32]): ptr cfloat =
  unsafeAddr(a[0])

proc toPtr(a: seq[Complex64]): ptr MKL_Complex16 =
  cast[ptr MKL_Complex16](unsafeAddr(a[0]))

proc toPtr(a: seq[Complex32]): ptr MKL_Complex8 =
  cast[ptr MKL_Complex8](unsafeAddr(a[0]))


proc dgesdd*(matrixLayout: TLayout; jobz: TSvdJob; m: int; n: int; a: seq[float64];
            lda: int; s: seq[float64]; u: seq[float64]; ldu: int; vt: seq[float64];
            ldvt: int): int =
  dgesdd(matrixLayout.cint, jobz.cchar, m.cint, n.cint, toPtr(a), lda.cint, toPtr(s),
    toPtr(u), ldu.cint, toPtr(vt), ldvt.cint)

proc sgesdd*(matrixLayout: TLayout; jobz: TSvdJob; m: int; n: int; a: seq[float32];
            lda: int; s: seq[float32]; u: seq[float32]; ldu: int; vt: seq[float32];
            ldvt: int): int =
  sgesdd(matrixLayout.cint, jobz.cchar, m.cint, n.cint, toPtr(a), lda.cint, toPtr(s),
    toPtr(u), ldu.cint, toPtr(vt), ldvt.cint)

proc dgeev*(matrixLayout: TLayout; jobvl: TEigJob; jobvr: TEigJob; n: int; a: seq[float64];
           lda: int; wr: seq[float64]; wi: seq[float64]; vl: seq[float64]; ldvl: int;
           vr: seq[float64]; ldvr: int): int =
  dgeev(matrixLayout.cint, jobvl.cchar, jobvr.cchar, n.cint, toPtr(a), lda.cint, toPtr(wr),
    toPtr(wi), toPtr(vl), ldvl.cint, toPtr(vr), ldvr.cint)

proc sgeev*(matrixLayout: TLayout; jobvl: TEigJob; jobvr: TEigJob; n: int; a: seq[float32];
           lda: int; wr: seq[float32]; wi: seq[float32]; vl: seq[float32]; ldvl: int;
           vr: seq[float32]; ldvr: int): int =
  sgeev(matrixLayout.cint, jobvl.cchar, jobvr.cchar, n.cint, toPtr(a), lda.cint, toPtr(wr),
    toPtr(wi), toPtr(vl), ldvl.cint, toPtr(vr), ldvr.cint)

proc dgetrf*(matrixLayout: TLayout; m: int; n: int; a: seq[float64]; lda: int;
            indices: seq[int]): int =
  dgetrf(matrixLayout.cint, m.cint, n.cint, toPtr(a), lda.cint, toPtr(indices))

proc sgetrf*(matrixLayout: TLayout; m: int; n: int; a: seq[float32]; lda: int;
            indices: seq[int]): int =
  sgetrf(matrixLayout.cint, m.cint, n.cint, toPtr(a), lda.cint, toPtr(indices))

proc dgeqrf*(matrixLayout: TLayout; m: int; n: int; a: seq[float64]; lda: int;
            tau: seq[float64]): int =
  dgeqrf(matrixLayout.cint, m.cint, n.cint, toPtr(a), lda.cint, toPtr(tau))

proc sgeqrf*(matrixLayout: TLayout; m: int; n: int; a: seq[float32]; lda: int;
            tau: seq[float32]): int =
  sgeqrf(matrixLayout.cint, m.cint, n.cint, toPtr(a), lda.cint, toPtr(tau))

proc dorgqr*(matrixLayout: TLayout; m: int; n: int; k: int; a: seq[float64];
            lda: int; tau: seq[float64]): int =
  dorgqr(matrixLayout.cint, m.cint, n.cint, k.cint, toPtr(a), lda.cint, toPtr(tau))

proc sorgqr*(matrixLayout: TLayout; m: int; n: int; k: int; a: seq[float32];
            lda: int; tau: seq[float32]): int =
  sorgqr(matrixLayout.cint, m.cint, n.cint, k.cint, toPtr(a), lda.cint, toPtr(tau))

proc dgesv*(matrixLayout: TLayout; n: int; nrhs: int; a: seq[float64]; lda: int;
           indices: seq[int]; b: seq[float64]; ldb: int): int =
  dgesv(matrixLayout.cint, n.cint, nrhs.cint, toPtr(a), lda.cint, toPtr(indices),
    toPtr(b), ldb.cint)

proc sgesv*(matrixLayout: TLayout; n: int; nrhs: int; a: seq[float32]; lda: int;
           indices: seq[int]; b: seq[float32]; ldb: int): int =
  sgesv(matrixLayout.cint, n.cint, nrhs.cint, toPtr(a), lda.cint, toPtr(indices),
    toPtr(b), ldb.cint)

proc dgels*(matrixLayout: TLayout; trans: TTrans; m: int; n: int; nrhs: int; a: seq[float64];
           lda: int; b: seq[float64]; ldb: int): int =
  dgels(matrixLayout.cint, trans.cchar, m.cint, n.cint, nrhs.cint, toPtr(a),
    lda.cint, toPtr(b), ldb.cint)

proc sgels*(matrixLayout: TLayout; trans: TTrans; m: int; n: int; nrhs: int; a: seq[float32];
           lda: int; b: seq[float32]; ldb: int): int =
  sgels(matrixLayout.cint, trans.cchar, m.cint, n.cint, nrhs.cint, toPtr(a),
    lda.cint, toPtr(b), ldb.cint)

proc dgemm*(layout: Cblas_Layout; transA: Cblas_Transpose; transB: Cblas_Transpose; m: int;
           n: int; k: int; alpha: float64; a: seq[float64]; lda: int; b: seq[float64];
           ldb: int; beta: float64; c: seq[float64]; ldc: int) =
  dgemm(layout, transA, transB, m.cint, n.cint, k.cint, alpha, toPtr(a),
    lda.cint, toPtr(b), ldb.cint, beta, toPtr(c), ldc.cint)

proc sgemm*(layout: Cblas_Layout; transA: Cblas_Transpose; transB: Cblas_Transpose; m: int;
           n: int; k: int; alpha: float32; a: seq[float32]; lda: int; b: seq[float32];
           ldb: int; beta: float32; c: seq[float32]; ldc: int) =
  sgemm(layout, transA, transB, m.cint, n.cint, k.cint, alpha, toPtr(a),
    lda.cint, toPtr(b), ldb.cint, beta, toPtr(c), ldc.cint)

proc zgemm3m*(layout: Cblas_Layout; transA: Cblas_Transpose; transB: Cblas_Transpose;
             m: int; n: int; k: int; alpha: Complex64; a: seq[Complex64]; lda: int;
             b: seq[Complex64]; ldb: int; beta: Complex64; c: seq[Complex64]; ldc: int) =
  var alph = MKL_Complex16(real: alpha.re, imag: alpha.im)
  var bet = MKL_Complex16(real: beta.re, imag: beta.im)
  zgemm3m(layout, transA, transB, m.cint, n.cint, k.cint, addr alph, toPtr(a),
    lda.cint, toPtr(b), ldb.cint, addr bet, toPtr(c), ldc.cint)

proc cgemm3m*(layout: Cblas_Layout; transA: Cblas_Transpose; transB: Cblas_Transpose;
             m: int; n: int; k: int; alpha: Complex32; a: seq[Complex32]; lda: int;
             b: seq[Complex32]; ldb: int; beta: Complex32; c: seq[Complex32]; ldc: int) =
  var alph = MKL_Complex8(real: alpha.re, imag: alpha.im)
  var bet = MKL_Complex8(real: beta.re, imag: beta.im)
  cgemm3m(layout, transA, transB, m.cint, n.cint, k.cint, addr alph, toPtr(a),
    lda.cint, toPtr(b), ldb.cint, addr bet, toPtr(c), ldc.cint)

proc zungqr*(matrixLayout: TLayout; m: int; n: int; k: int; a: seq[Complex64];
            lda: int; tau: seq[Complex64]): int =
  zungqr(matrixLayout.cint, m.cint, n.cint, k.cint, toPtr(a), lda.cint, toPtr(tau))

proc cungqr*(matrixLayout: TLayout; m: int; n: int; k: int; a: seq[Complex32];
            lda: int; tau: seq[Complex32]): int =
  cungqr(matrixLayout.cint, m.cint, n.cint, k.cint, toPtr(a), lda.cint, toPtr(tau))

proc zgeqrf*(matrixLayout: TLayout; m: int; n: int; a: seq[Complex64]; lda: int;
            tau: seq[Complex64]): int =
  zgeqrf(matrixLayout.cint, m.cint, n.cint, toPtr(a), lda.cint, toPtr(tau))

proc cgeqrf*(matrixLayout: TLayout; m: int; n: int; a: seq[Complex32]; lda: int;
            tau: seq[Complex32]): int =
  cgeqrf(matrixLayout.cint, m.cint, n.cint, toPtr(a), lda.cint, toPtr(tau))

proc zgetrf*(matrixLayout: TLayout; m: int; n: int; a: seq[Complex64]; lda: int;
            indices: seq[int]): int =
  zgetrf(matrixLayout.cint, m.cint, n.cint, toPtr(a), lda.cint, toPtr(indices))

proc cgetrf*(matrixLayout: TLayout; m: int; n: int; a: seq[Complex32]; lda: int;
            indices: seq[int]): int =
  cgetrf(matrixLayout.cint, m.cint, n.cint, toPtr(a), lda.cint, toPtr(indices))

proc zgeev*(matrixLayout: TLayout; jobvl: TEigJob; jobvr: TEigJob; n: int; a: seq[Complex64];
           lda: int; w: seq[Complex64]; vl: seq[Complex64]; ldvl: int; vr: seq[Complex64];
           ldvr: int): int =
  zgeev(matrixLayout.cint, jobvl.cchar, jobvr.cchar, n.cint, toPtr(a), lda.cint,
    toPtr(w), toPtr(vl), ldvl.cint, toPtr(vr), ldvr.cint)

proc cgeev*(matrixLayout: TLayout; jobvl: TEigJob; jobvr: TEigJob; n: int; a: seq[Complex32];
           lda: int; w: seq[Complex32]; vl: seq[Complex32]; ldvl: int; vr: seq[Complex32];
           ldvr: int): int =
  cgeev(matrixLayout.cint, jobvl.cchar, jobvr.cchar, n.cint, toPtr(a), lda.cint,
    toPtr(w), toPtr(vl), ldvl.cint, toPtr(vr), ldvr.cint)

proc zgels*(matrixLayout: TLayout; trans: TTrans; m: int; n: int; nrhs: int;
           a: seq[Complex64]; lda: int; b: seq[Complex64]; ldb: int): int =
  zgels(matrixLayout.cint, trans.cchar, m.cint, n.cint, nrhs.cint, toPtr(a), lda.cint,
    toPtr(b), ldb.cint)

proc cgels*(matrixLayout: TLayout; trans: TTrans; m: int; n: int; nrhs: int;
           a: seq[Complex32]; lda: int; b: seq[Complex32]; ldb: int): int =
  cgels(matrixLayout.cint, trans.cchar, m.cint, n.cint, nrhs.cint, toPtr(a), lda.cint,
    toPtr(b), ldb.cint)

proc zgesv*(matrixLayout: TLayout; n: int; nrhs: int; a: seq[Complex64]; lda: int;
           indices: seq[int]; b: seq[Complex64]; ldb: int): int =
  zgesv(matrixLayout.cint, n.cint, nrhs.cint, toPtr(a), lda.cint, toPtr(indices), toPtr(b), ldb.cint)

proc cgesv*(matrixLayout: TLayout; n: int; nrhs: int; a: seq[Complex32]; lda: int;
           indices: seq[int]; b: seq[Complex32]; ldb: int): int =
  cgesv(matrixLayout.cint, n.cint, nrhs.cint, toPtr(a), lda.cint, toPtr(indices), toPtr(b), ldb.cint)

proc zgesdd*(matrixLayout: TLayout; jobz: TSvdJob; m: int; n: int; a: seq[Complex64];
            lda: int; s: seq[float64]; u: seq[Complex64]; ldu: int; vt: seq[Complex64];
            ldvt: int): int =
  zgesdd(matrixLayout.cint, jobz.cchar, m.cint, n.cint, toPtr(a), lda.cint, toPtr(s),
    toPtr(u), ldu.cint, toPtr(vt), ldvt.cint)

proc cgesdd*(matrixLayout: TLayout; jobz: TSvdJob; m: int; n: int; a: seq[Complex32];
            lda: int; s: seq[float32]; u: seq[Complex32]; ldu: int; vt: seq[Complex32];
            ldvt: int): int =
  cgesdd(matrixLayout.cint, jobz.cchar, m.cint, n.cint, toPtr(a), lda.cint, toPtr(s),
    toPtr(u), ldu.cint, toPtr(vt), ldvt.cint)


