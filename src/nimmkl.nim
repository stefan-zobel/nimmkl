# Copyright (c) 2022 Stefan Zobel
# Distributed under the Apache v2 License (license terms are at https://www.apache.org/licenses/LICENSE-2.0).

import nimmkl/[mklTypes, mklTrans, mklCblas, mklLapacke, mklService]

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


