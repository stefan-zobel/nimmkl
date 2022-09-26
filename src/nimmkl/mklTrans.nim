import mklTypes

##  In-place transposition routines

proc simatcopy*(ordering: char; trans: char; rows: csize_t; cols: csize_t;
               alpha: cfloat; ab: ptr cfloat; lda: csize_t; ldb: csize_t) {.
    importc: "MKL_Simatcopy", dynlib: libName.}

proc dimatcopy*(ordering: char; trans: char; rows: csize_t; cols: csize_t;
               alpha: cdouble; ab: ptr cdouble; lda: csize_t; ldb: csize_t) {.
    importc: "MKL_Dimatcopy", dynlib: libName.}

proc cimatcopy*(ordering: char; trans: char; rows: csize_t; cols: csize_t;
               alpha: MKL_Complex8; ab: ptr MKL_Complex8; lda: csize_t; ldb: csize_t) {.
    importc: "MKL_Cimatcopy", dynlib: libName.}

proc zimatcopy*(ordering: char; trans: char; rows: csize_t; cols: csize_t;
               alpha: MKL_Complex16; ab: ptr MKL_Complex16; lda: csize_t; ldb: csize_t) {.
    importc: "MKL_Zimatcopy", dynlib: libName.}

##  Out-of-place transposition routines

proc somatcopy*(ordering: char; trans: char; rows: csize_t; cols: csize_t;
               alpha: cfloat; a: ptr cfloat; lda: csize_t; b: ptr cfloat; ldb: csize_t) {.
    importc: "MKL_Somatcopy", dynlib: libName.}

proc domatcopy*(ordering: char; trans: char; rows: csize_t; cols: csize_t;
               alpha: cdouble; a: ptr cdouble; lda: csize_t; b: ptr cdouble;
               ldb: csize_t) {.importc: "MKL_Domatcopy", dynlib: libName.}

proc comatcopy*(ordering: char; trans: char; rows: csize_t; cols: csize_t;
               alpha: MKL_Complex8; a: ptr MKL_Complex8; lda: csize_t; b: ptr MKL_Complex8;
               ldb: csize_t) {.importc: "MKL_Comatcopy", dynlib: libName.}

proc zomatcopy*(ordering: char; trans: char; rows: csize_t; cols: csize_t;
               alpha: MKL_Complex16; a: ptr MKL_Complex16; lda: csize_t; b: ptr MKL_Complex16;
               ldb: csize_t) {.importc: "MKL_Zomatcopy", dynlib: libName.}

##  Out-of-place transposition routines (all-strided case)

proc somatcopy2*(ordering: char; trans: char; rows: csize_t; cols: csize_t;
                alpha: cfloat; a: ptr cfloat; lda: csize_t; stridea: csize_t;
                b: ptr cfloat; ldb: csize_t; strideb: csize_t) {.
    importc: "MKL_Somatcopy2", dynlib: libName.}

proc domatcopy2*(ordering: char; trans: char; rows: csize_t; cols: csize_t;
                alpha: cdouble; a: ptr cdouble; lda: csize_t; stridea: csize_t;
                b: ptr cdouble; ldb: csize_t; strideb: csize_t) {.
    importc: "MKL_Domatcopy2", dynlib: libName.}

proc comatcopy2*(ordering: char; trans: char; rows: csize_t; cols: csize_t;
                alpha: MKL_Complex8; a: ptr MKL_Complex8; lda: csize_t; stridea: csize_t;
                b: ptr MKL_Complex8; ldb: csize_t; strideb: csize_t) {.
    importc: "MKL_Comatcopy2", dynlib: libName.}

proc zomatcopy2*(ordering: char; trans: char; rows: csize_t; cols: csize_t;
                alpha: MKL_Complex16; a: ptr MKL_Complex16; lda: csize_t; stridea: csize_t;
                b: ptr MKL_Complex16; ldb: csize_t; strideb: csize_t) {.
    importc: "MKL_Zomatcopy2", dynlib: libName.}

##  Out-of-place memory movement routines

proc somatadd*(ordering: char; transa: char; transb: char; rows: csize_t;
              cols: csize_t; alpha: cfloat; a: ptr cfloat; lda: csize_t; beta: cfloat;
              b: ptr cfloat; ldb: csize_t; c: ptr cfloat; ldc: csize_t) {.
    importc: "MKL_Somatadd", dynlib: libName.}

proc domatadd*(ordering: char; transa: char; transb: char; rows: csize_t;
              cols: csize_t; alpha: cdouble; a: ptr cdouble; lda: csize_t;
              beta: cdouble; b: ptr cdouble; ldb: csize_t; c: ptr cdouble; ldc: csize_t) {.
    importc: "MKL_Domatadd", dynlib: libName.}

proc comatadd*(ordering: char; transa: char; transb: char; rows: csize_t;
              cols: csize_t; alpha: MKL_Complex8; a: ptr MKL_Complex8; lda: csize_t;
              beta: MKL_Complex8; b: ptr MKL_Complex8; ldb: csize_t; c: ptr MKL_Complex8;
              ldc: csize_t) {.importc: "MKL_Comatadd", dynlib: libName.}

proc zomatadd*(ordering: char; transa: char; transb: char; rows: csize_t;
              cols: csize_t; alpha: MKL_Complex16; a: ptr MKL_Complex16; lda: csize_t;
              beta: MKL_Complex16; b: ptr MKL_Complex16; ldb: csize_t; c: ptr MKL_Complex16;
              ldc: csize_t) {.importc: "MKL_Zomatadd", dynlib: libName.}
