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
