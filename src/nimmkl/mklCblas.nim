import mklTypes

##
##  ===========================================================================
##  Prototypes for level 1 BLAS functions (complex are recast as routines)
##  ===========================================================================
##

proc dcabs1*(z: pointer): cdouble {.importc: "cblas_dcabs1", dynlib: libName.}
proc scabs1*(c: pointer): cfloat {.importc: "cblas_scabs1", dynlib: libName.}
proc sdot*(n: Mkl_Int; x: ptr cfloat; incX: Mkl_Int; y: ptr cfloat; incY: Mkl_Int): cfloat {.
    importc: "cblas_sdot", dynlib: libName.}
proc sdoti*(n: Mkl_Int; x: ptr cfloat; indx: ptr Mkl_Int; y: ptr cfloat): cfloat {.
    importc: "cblas_sdoti", dynlib: libName.}
proc ddot*(n: Mkl_Int; x: ptr cdouble; incX: Mkl_Int; y: ptr cdouble; incY: Mkl_Int): cdouble {.
    importc: "cblas_ddot", dynlib: libName.}
proc ddoti*(n: Mkl_Int; x: ptr cdouble; indx: ptr Mkl_Int; y: ptr cdouble): cdouble {.
    importc: "cblas_ddoti", dynlib: libName.}
proc dsdot*(n: Mkl_Int; x: ptr cfloat; incX: Mkl_Int; y: ptr cfloat; incY: Mkl_Int): cdouble {.
    importc: "cblas_dsdot", dynlib: libName.}
proc sdsdot*(n: Mkl_Int; sb: cfloat; x: ptr cfloat; incX: Mkl_Int; y: ptr cfloat;
            incY: Mkl_Int): cfloat {.importc: "cblas_sdsdot", dynlib: libName.}

##
##  Functions having prefixes Z and C only
##

proc cdotuSub*(n: Mkl_Int; x: pointer; incX: Mkl_Int; y: pointer; incY: Mkl_Int;
              dotu: pointer) {.importc: "cblas_cdotu_sub", dynlib: libName.}
proc cdotuiSub*(n: Mkl_Int; x: pointer; indx: ptr Mkl_Int; y: pointer; dotui: pointer) {.
    importc: "cblas_cdotui_sub", dynlib: libName.}
proc cdotcSub*(n: Mkl_Int; x: pointer; incX: Mkl_Int; y: pointer; incY: Mkl_Int;
              dotc: pointer) {.importc: "cblas_cdotc_sub", dynlib: libName.}
proc cdotciSub*(n: Mkl_Int; x: pointer; indx: ptr Mkl_Int; y: pointer; dotui: pointer) {.
    importc: "cblas_cdotci_sub", dynlib: libName.}
proc zdotuSub*(n: Mkl_Int; x: pointer; incX: Mkl_Int; y: pointer; incY: Mkl_Int;
              dotu: pointer) {.importc: "cblas_zdotu_sub", dynlib: libName.}
proc zdotuiSub*(n: Mkl_Int; x: pointer; indx: ptr Mkl_Int; y: pointer; dotui: pointer) {.
    importc: "cblas_zdotui_sub", dynlib: libName.}
proc zdotcSub*(n: Mkl_Int; x: pointer; incX: Mkl_Int; y: pointer; incY: Mkl_Int;
              dotc: pointer) {.importc: "cblas_zdotc_sub", dynlib: libName.}
proc zdotciSub*(n: Mkl_Int; x: pointer; indx: ptr Mkl_Int; y: pointer; dotui: pointer) {.
    importc: "cblas_zdotci_sub", dynlib: libName.}

##
##  Functions having prefixes S D SC DZ
##

proc snrm2*(n: Mkl_Int; x: ptr cfloat; incX: Mkl_Int): cfloat {.importc: "cblas_snrm2",
    dynlib: libName.}
proc sasum*(n: Mkl_Int; x: ptr cfloat; incX: Mkl_Int): cfloat {.importc: "cblas_sasum",
    dynlib: libName.}
proc dnrm2*(n: Mkl_Int; x: ptr cdouble; incX: Mkl_Int): cdouble {.importc: "cblas_dnrm2",
    dynlib: libName.}
proc dasum*(n: Mkl_Int; x: ptr cdouble; incX: Mkl_Int): cdouble {.importc: "cblas_dasum",
    dynlib: libName.}
proc scnrm2*(n: Mkl_Int; x: pointer; incX: Mkl_Int): cfloat {.importc: "cblas_scnrm2",
    dynlib: libName.}
proc scasum*(n: Mkl_Int; x: pointer; incX: Mkl_Int): cfloat {.importc: "cblas_scasum",
    dynlib: libName.}
proc dznrm2*(n: Mkl_Int; x: pointer; incX: Mkl_Int): cdouble {.importc: "cblas_dznrm2",
    dynlib: libName.}
proc dzasum*(n: Mkl_Int; x: pointer; incX: Mkl_Int): cdouble {.importc: "cblas_dzasum",
    dynlib: libName.}

##
##  Functions having standard 4 prefixes (S D C Z)
##

proc isamax*(n: Mkl_Int; x: ptr cfloat; incX: Mkl_Int): Cblas_Index {.
    importc: "cblas_isamax", dynlib: libName.}
proc idamax*(n: Mkl_Int; x: ptr cdouble; incX: Mkl_Int): Cblas_Index {.
    importc: "cblas_idamax", dynlib: libName.}
proc icamax*(n: Mkl_Int; x: pointer; incX: Mkl_Int): Cblas_Index {.
    importc: "cblas_icamax", dynlib: libName.}
proc izamax*(n: Mkl_Int; x: pointer; incX: Mkl_Int): Cblas_Index {.
    importc: "cblas_izamax", dynlib: libName.}
proc isamin*(n: Mkl_Int; x: ptr cfloat; incX: Mkl_Int): Cblas_Index {.
    importc: "cblas_isamin", dynlib: libName.}
proc idamin*(n: Mkl_Int; x: ptr cdouble; incX: Mkl_Int): Cblas_Index {.
    importc: "cblas_idamin", dynlib: libName.}
proc icamin*(n: Mkl_Int; x: pointer; incX: Mkl_Int): Cblas_Index {.
    importc: "cblas_icamin", dynlib: libName.}
proc izamin*(n: Mkl_Int; x: pointer; incX: Mkl_Int): Cblas_Index {.
    importc: "cblas_izamin", dynlib: libName.}

##
##  ===========================================================================
##  Prototypes for level 1 BLAS routines
##  ===========================================================================
##
##
##  Routines with standard 4 prefixes (s, d, c, z)
##

proc sswap*(n: Mkl_Int; x: ptr cfloat; incX: Mkl_Int; y: ptr cfloat; incY: Mkl_Int) {.
    importc: "cblas_sswap", dynlib: libName.}
proc scopy*(n: Mkl_Int; x: ptr cfloat; incX: Mkl_Int; y: ptr cfloat; incY: Mkl_Int) {.
    importc: "cblas_scopy", dynlib: libName.}
proc saxpy*(n: Mkl_Int; alpha: cfloat; x: ptr cfloat; incX: Mkl_Int; y: ptr cfloat;
           incY: Mkl_Int) {.importc: "cblas_saxpy", dynlib: libName.}
proc saxpby*(n: Mkl_Int; alpha: cfloat; x: ptr cfloat; incX: Mkl_Int; beta: cfloat;
            y: ptr cfloat; incY: Mkl_Int) {.importc: "cblas_saxpby", dynlib: libName.}
proc saxpyi*(n: Mkl_Int; alpha: cfloat; x: ptr cfloat; indx: ptr Mkl_Int; y: ptr cfloat) {.
    importc: "cblas_saxpyi", dynlib: libName.}
proc sgthr*(n: Mkl_Int; y: ptr cfloat; x: ptr cfloat; indx: ptr Mkl_Int) {.
    importc: "cblas_sgthr", dynlib: libName.}
proc sgthrz*(n: Mkl_Int; y: ptr cfloat; x: ptr cfloat; indx: ptr Mkl_Int) {.
    importc: "cblas_sgthrz", dynlib: libName.}
proc ssctr*(n: Mkl_Int; x: ptr cfloat; indx: ptr Mkl_Int; y: ptr cfloat) {.
    importc: "cblas_ssctr", dynlib: libName.}
proc srotg*(a: ptr cfloat; b: ptr cfloat; c: ptr cfloat; s: ptr cfloat) {.
    importc: "cblas_srotg", dynlib: libName.}
proc dswap*(n: Mkl_Int; x: ptr cdouble; incX: Mkl_Int; y: ptr cdouble; incY: Mkl_Int) {.
    importc: "cblas_dswap", dynlib: libName.}
proc dcopy*(n: Mkl_Int; x: ptr cdouble; incX: Mkl_Int; y: ptr cdouble; incY: Mkl_Int) {.
    importc: "cblas_dcopy", dynlib: libName.}
proc daxpy*(n: Mkl_Int; alpha: cdouble; x: ptr cdouble; incX: Mkl_Int; y: ptr cdouble;
           incY: Mkl_Int) {.importc: "cblas_daxpy", dynlib: libName.}
proc daxpby*(n: Mkl_Int; alpha: cdouble; x: ptr cdouble; incX: Mkl_Int; beta: cdouble;
            y: ptr cdouble; incY: Mkl_Int) {.importc: "cblas_daxpby", dynlib: libName.}
proc daxpyi*(n: Mkl_Int; alpha: cdouble; x: ptr cdouble; indx: ptr Mkl_Int; y: ptr cdouble) {.
    importc: "cblas_daxpyi", dynlib: libName.}
proc dgthr*(n: Mkl_Int; y: ptr cdouble; x: ptr cdouble; indx: ptr Mkl_Int) {.
    importc: "cblas_dgthr", dynlib: libName.}
proc dgthrz*(n: Mkl_Int; y: ptr cdouble; x: ptr cdouble; indx: ptr Mkl_Int) {.
    importc: "cblas_dgthrz", dynlib: libName.}
proc dsctr*(n: Mkl_Int; x: ptr cdouble; indx: ptr Mkl_Int; y: ptr cdouble) {.
    importc: "cblas_dsctr", dynlib: libName.}
proc drotg*(a: ptr cdouble; b: ptr cdouble; c: ptr cdouble; s: ptr cdouble) {.
    importc: "cblas_drotg", dynlib: libName.}
proc cswap*(n: Mkl_Int; x: pointer; incX: Mkl_Int; y: pointer; incY: Mkl_Int) {.
    importc: "cblas_cswap", dynlib: libName.}
proc ccopy*(n: Mkl_Int; x: pointer; incX: Mkl_Int; y: pointer; incY: Mkl_Int) {.
    importc: "cblas_ccopy", dynlib: libName.}
proc caxpy*(n: Mkl_Int; alpha: pointer; x: pointer; incX: Mkl_Int; y: pointer;
           incY: Mkl_Int) {.importc: "cblas_caxpy", dynlib: libName.}
proc caxpby*(n: Mkl_Int; alpha: pointer; x: pointer; incX: Mkl_Int; beta: pointer;
            y: pointer; incY: Mkl_Int) {.importc: "cblas_caxpby", dynlib: libName.}
proc caxpyi*(n: Mkl_Int; alpha: pointer; x: pointer; indx: ptr Mkl_Int; y: pointer) {.
    importc: "cblas_caxpyi", dynlib: libName.}
proc cgthr*(n: Mkl_Int; y: pointer; x: pointer; indx: ptr Mkl_Int) {.
    importc: "cblas_cgthr", dynlib: libName.}
proc cgthrz*(n: Mkl_Int; y: pointer; x: pointer; indx: ptr Mkl_Int) {.
    importc: "cblas_cgthrz", dynlib: libName.}
proc csctr*(n: Mkl_Int; x: pointer; indx: ptr Mkl_Int; y: pointer) {.
    importc: "cblas_csctr", dynlib: libName.}
proc crotg*(a: pointer; b: pointer; c: ptr cfloat; s: pointer) {.importc: "cblas_crotg",
    dynlib: libName.}
proc zswap*(n: Mkl_Int; x: pointer; incX: Mkl_Int; y: pointer; incY: Mkl_Int) {.
    importc: "cblas_zswap", dynlib: libName.}
proc zcopy*(n: Mkl_Int; x: pointer; incX: Mkl_Int; y: pointer; incY: Mkl_Int) {.
    importc: "cblas_zcopy", dynlib: libName.}
proc zaxpy*(n: Mkl_Int; alpha: pointer; x: pointer; incX: Mkl_Int; y: pointer;
           incY: Mkl_Int) {.importc: "cblas_zaxpy", dynlib: libName.}
proc zaxpby*(n: Mkl_Int; alpha: pointer; x: pointer; incX: Mkl_Int; beta: pointer;
            y: pointer; incY: Mkl_Int) {.importc: "cblas_zaxpby", dynlib: libName.}
proc zaxpyi*(n: Mkl_Int; alpha: pointer; x: pointer; indx: ptr Mkl_Int; y: pointer) {.
    importc: "cblas_zaxpyi", dynlib: libName.}
proc zgthr*(n: Mkl_Int; y: pointer; x: pointer; indx: ptr Mkl_Int) {.
    importc: "cblas_zgthr", dynlib: libName.}
proc zgthrz*(n: Mkl_Int; y: pointer; x: pointer; indx: ptr Mkl_Int) {.
    importc: "cblas_zgthrz", dynlib: libName.}
proc zsctr*(n: Mkl_Int; x: pointer; indx: ptr Mkl_Int; y: pointer) {.
    importc: "cblas_zsctr", dynlib: libName.}
proc zrotg*(a: pointer; b: pointer; c: ptr cdouble; s: pointer) {.importc: "cblas_zrotg",
    dynlib: libName.}

##
##  Routines with S and D prefix only
##

proc srotmg*(d1: ptr cfloat; d2: ptr cfloat; b1: ptr cfloat; b2: cfloat; p: ptr cfloat) {.
    importc: "cblas_srotmg", dynlib: libName.}
proc srot*(n: Mkl_Int; x: ptr cfloat; incX: Mkl_Int; y: ptr cfloat; incY: Mkl_Int; c: cfloat;
          s: cfloat) {.importc: "cblas_srot", dynlib: libName.}
proc sroti*(n: Mkl_Int; x: ptr cfloat; indx: ptr Mkl_Int; y: ptr cfloat; c: cfloat; s: cfloat) {.
    importc: "cblas_sroti", dynlib: libName.}
proc srotm*(n: Mkl_Int; x: ptr cfloat; incX: Mkl_Int; y: ptr cfloat; incY: Mkl_Int;
           p: ptr cfloat) {.importc: "cblas_srotm", dynlib: libName.}
proc drotmg*(d1: ptr cdouble; d2: ptr cdouble; b1: ptr cdouble; b2: cdouble; p: ptr cdouble) {.
    importc: "cblas_drotmg", dynlib: libName.}
proc drot*(n: Mkl_Int; x: ptr cdouble; incX: Mkl_Int; y: ptr cdouble; incY: Mkl_Int;
          c: cdouble; s: cdouble) {.importc: "cblas_drot", dynlib: libName.}
proc drotm*(n: Mkl_Int; x: ptr cdouble; incX: Mkl_Int; y: ptr cdouble; incY: Mkl_Int;
           p: ptr cdouble) {.importc: "cblas_drotm", dynlib: libName.}
proc droti*(n: Mkl_Int; x: ptr cdouble; indx: ptr Mkl_Int; y: ptr cdouble; c: cdouble;
           s: cdouble) {.importc: "cblas_droti", dynlib: libName.}

##
##  Routines with CS and ZD prefix only
##

proc csrot*(n: Mkl_Int; x: pointer; incX: Mkl_Int; y: pointer; incY: Mkl_Int; c: cfloat;
           s: cfloat) {.importc: "cblas_csrot", dynlib: libName.}
proc zdrot*(n: Mkl_Int; x: pointer; incX: Mkl_Int; y: pointer; incY: Mkl_Int; c: cdouble;
           s: cdouble) {.importc: "cblas_zdrot", dynlib: libName.}

##
##  Routines with S D C Z CS and ZD prefixes
##

proc sscal*(n: Mkl_Int; alpha: cfloat; x: ptr cfloat; incX: Mkl_Int) {.
    importc: "cblas_sscal", dynlib: libName.}
proc dscal*(n: Mkl_Int; alpha: cdouble; x: ptr cdouble; incX: Mkl_Int) {.
    importc: "cblas_dscal", dynlib: libName.}
proc cscal*(n: Mkl_Int; alpha: pointer; x: pointer; incX: Mkl_Int) {.
    importc: "cblas_cscal", dynlib: libName.}
proc zscal*(n: Mkl_Int; alpha: pointer; x: pointer; incX: Mkl_Int) {.
    importc: "cblas_zscal", dynlib: libName.}
proc csscal*(n: Mkl_Int; alpha: cfloat; x: pointer; incX: Mkl_Int) {.
    importc: "cblas_csscal", dynlib: libName.}
proc zdscal*(n: Mkl_Int; alpha: cdouble; x: pointer; incX: Mkl_Int) {.
    importc: "cblas_zdscal", dynlib: libName.}

##
##  ===========================================================================
##  Prototypes for level 2 BLAS
##  ===========================================================================
##
##
##  Routines with standard 4 prefixes (S, D, C, Z)
##

proc sgemv*(layout: Cblas_Layout; transA: Cblas_Transpose; m: Mkl_Int; n: Mkl_Int;
           alpha: cfloat; a: ptr cfloat; lda: Mkl_Int; x: ptr cfloat; incX: Mkl_Int;
           beta: cfloat; y: ptr cfloat; incY: Mkl_Int) {.importc: "cblas_sgemv",
    dynlib: libName.}
proc sgbmv*(layout: Cblas_Layout; transA: Cblas_Transpose; m: Mkl_Int; n: Mkl_Int;
           kl: Mkl_Int; ku: Mkl_Int; alpha: cfloat; a: ptr cfloat; lda: Mkl_Int;
           x: ptr cfloat; incX: Mkl_Int; beta: cfloat; y: ptr cfloat; incY: Mkl_Int) {.
    importc: "cblas_sgbmv", dynlib: libName.}
proc strmv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; a: ptr cfloat; lda: Mkl_Int; x: ptr cfloat;
           incX: Mkl_Int) {.importc: "cblas_strmv", dynlib: libName.}
proc stbmv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; k: Mkl_Int; a: ptr cfloat; lda: Mkl_Int;
           x: ptr cfloat; incX: Mkl_Int) {.importc: "cblas_stbmv", dynlib: libName.}
proc stpmv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; ap: ptr cfloat; x: ptr cfloat; incX: Mkl_Int) {.
    importc: "cblas_stpmv", dynlib: libName.}
proc strsv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; a: ptr cfloat; lda: Mkl_Int; x: ptr cfloat;
           incX: Mkl_Int) {.importc: "cblas_strsv", dynlib: libName.}
proc stbsv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; k: Mkl_Int; a: ptr cfloat; lda: Mkl_Int;
           x: ptr cfloat; incX: Mkl_Int) {.importc: "cblas_stbsv", dynlib: libName.}
proc stpsv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; ap: ptr cfloat; x: ptr cfloat; incX: Mkl_Int) {.
    importc: "cblas_stpsv", dynlib: libName.}
proc dgemv*(layout: Cblas_Layout; transA: Cblas_Transpose; m: Mkl_Int; n: Mkl_Int;
           alpha: cdouble; a: ptr cdouble; lda: Mkl_Int; x: ptr cdouble; incX: Mkl_Int;
           beta: cdouble; y: ptr cdouble; incY: Mkl_Int) {.importc: "cblas_dgemv",
    dynlib: libName.}
proc dgbmv*(layout: Cblas_Layout; transA: Cblas_Transpose; m: Mkl_Int; n: Mkl_Int;
           kl: Mkl_Int; ku: Mkl_Int; alpha: cdouble; a: ptr cdouble; lda: Mkl_Int;
           x: ptr cdouble; incX: Mkl_Int; beta: cdouble; y: ptr cdouble; incY: Mkl_Int) {.
    importc: "cblas_dgbmv", dynlib: libName.}
proc dtrmv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; a: ptr cdouble; lda: Mkl_Int; x: ptr cdouble;
           incX: Mkl_Int) {.importc: "cblas_dtrmv", dynlib: libName.}
proc dtbmv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; k: Mkl_Int; a: ptr cdouble; lda: Mkl_Int;
           x: ptr cdouble; incX: Mkl_Int) {.importc: "cblas_dtbmv", dynlib: libName.}
proc dtpmv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; ap: ptr cdouble; x: ptr cdouble; incX: Mkl_Int) {.
    importc: "cblas_dtpmv", dynlib: libName.}
proc dtrsv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; a: ptr cdouble; lda: Mkl_Int; x: ptr cdouble;
           incX: Mkl_Int) {.importc: "cblas_dtrsv", dynlib: libName.}
proc dtbsv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; k: Mkl_Int; a: ptr cdouble; lda: Mkl_Int;
           x: ptr cdouble; incX: Mkl_Int) {.importc: "cblas_dtbsv", dynlib: libName.}
proc dtpsv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; ap: ptr cdouble; x: ptr cdouble; incX: Mkl_Int) {.
    importc: "cblas_dtpsv", dynlib: libName.}
proc cgemv*(layout: Cblas_Layout; transA: Cblas_Transpose; m: Mkl_Int; n: Mkl_Int;
           alpha: pointer; a: pointer; lda: Mkl_Int; x: pointer; incX: Mkl_Int;
           beta: pointer; y: pointer; incY: Mkl_Int) {.importc: "cblas_cgemv",
    dynlib: libName.}
proc cgbmv*(layout: Cblas_Layout; transA: Cblas_Transpose; m: Mkl_Int; n: Mkl_Int;
           kl: Mkl_Int; ku: Mkl_Int; alpha: pointer; a: pointer; lda: Mkl_Int; x: pointer;
           incX: Mkl_Int; beta: pointer; y: pointer; incY: Mkl_Int) {.
    importc: "cblas_cgbmv", dynlib: libName.}
proc ctrmv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; a: pointer; lda: Mkl_Int; x: pointer;
           incX: Mkl_Int) {.importc: "cblas_ctrmv", dynlib: libName.}
proc ctbmv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; k: Mkl_Int; a: pointer; lda: Mkl_Int; x: pointer;
           incX: Mkl_Int) {.importc: "cblas_ctbmv", dynlib: libName.}
proc ctpmv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; ap: pointer; x: pointer; incX: Mkl_Int) {.
    importc: "cblas_ctpmv", dynlib: libName.}
proc ctrsv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; a: pointer; lda: Mkl_Int; x: pointer;
           incX: Mkl_Int) {.importc: "cblas_ctrsv", dynlib: libName.}
proc ctbsv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; k: Mkl_Int; a: pointer; lda: Mkl_Int; x: pointer;
           incX: Mkl_Int) {.importc: "cblas_ctbsv", dynlib: libName.}
proc ctpsv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; ap: pointer; x: pointer; incX: Mkl_Int) {.
    importc: "cblas_ctpsv", dynlib: libName.}
proc zgemv*(layout: Cblas_Layout; transA: Cblas_Transpose; m: Mkl_Int; n: Mkl_Int;
           alpha: pointer; a: pointer; lda: Mkl_Int; x: pointer; incX: Mkl_Int;
           beta: pointer; y: pointer; incY: Mkl_Int) {.importc: "cblas_zgemv",
    dynlib: libName.}
proc zgbmv*(layout: Cblas_Layout; transA: Cblas_Transpose; m: Mkl_Int; n: Mkl_Int;
           kl: Mkl_Int; ku: Mkl_Int; alpha: pointer; a: pointer; lda: Mkl_Int; x: pointer;
           incX: Mkl_Int; beta: pointer; y: pointer; incY: Mkl_Int) {.
    importc: "cblas_zgbmv", dynlib: libName.}
proc ztrmv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; a: pointer; lda: Mkl_Int; x: pointer;
           incX: Mkl_Int) {.importc: "cblas_ztrmv", dynlib: libName.}
proc ztbmv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; k: Mkl_Int; a: pointer; lda: Mkl_Int; x: pointer;
           incX: Mkl_Int) {.importc: "cblas_ztbmv", dynlib: libName.}
proc ztpmv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; ap: pointer; x: pointer; incX: Mkl_Int) {.
    importc: "cblas_ztpmv", dynlib: libName.}
proc ztrsv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; a: pointer; lda: Mkl_Int; x: pointer;
           incX: Mkl_Int) {.importc: "cblas_ztrsv", dynlib: libName.}
proc ztbsv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; k: Mkl_Int; a: pointer; lda: Mkl_Int; x: pointer;
           incX: Mkl_Int) {.importc: "cblas_ztbsv", dynlib: libName.}
proc ztpsv*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
           diag: Cblas_Diag; n: Mkl_Int; ap: pointer; x: pointer; incX: Mkl_Int) {.
    importc: "cblas_ztpsv", dynlib: libName.}

##
##  Routines with S and D prefixes only
##

proc ssymv*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: cfloat;
           a: ptr cfloat; lda: Mkl_Int; x: ptr cfloat; incX: Mkl_Int; beta: cfloat;
           y: ptr cfloat; incY: Mkl_Int) {.importc: "cblas_ssymv", dynlib: libName.}
proc ssbmv*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; k: Mkl_Int;
           alpha: cfloat; a: ptr cfloat; lda: Mkl_Int; x: ptr cfloat; incX: Mkl_Int;
           beta: cfloat; y: ptr cfloat; incY: Mkl_Int) {.importc: "cblas_ssbmv",
    dynlib: libName.}
proc sspmv*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: cfloat;
           ap: ptr cfloat; x: ptr cfloat; incX: Mkl_Int; beta: cfloat; y: ptr cfloat;
           incY: Mkl_Int) {.importc: "cblas_sspmv", dynlib: libName.}
proc sger*(layout: Cblas_Layout; m: Mkl_Int; n: Mkl_Int; alpha: cfloat; x: ptr cfloat;
          incX: Mkl_Int; y: ptr cfloat; incY: Mkl_Int; a: ptr cfloat; lda: Mkl_Int) {.
    importc: "cblas_sger", dynlib: libName.}
proc ssyr*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: cfloat;
          x: ptr cfloat; incX: Mkl_Int; a: ptr cfloat; lda: Mkl_Int) {.
    importc: "cblas_ssyr", dynlib: libName.}
proc sspr*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: cfloat;
          x: ptr cfloat; incX: Mkl_Int; ap: ptr cfloat) {.importc: "cblas_sspr",
    dynlib: libName.}
proc ssyr2*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: cfloat;
           x: ptr cfloat; incX: Mkl_Int; y: ptr cfloat; incY: Mkl_Int; a: ptr cfloat;
           lda: Mkl_Int) {.importc: "cblas_ssyr2", dynlib: libName.}
proc sspr2*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: cfloat;
           x: ptr cfloat; incX: Mkl_Int; y: ptr cfloat; incY: Mkl_Int; a: ptr cfloat) {.
    importc: "cblas_sspr2", dynlib: libName.}
proc dsymv*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: cdouble;
           a: ptr cdouble; lda: Mkl_Int; x: ptr cdouble; incX: Mkl_Int; beta: cdouble;
           y: ptr cdouble; incY: Mkl_Int) {.importc: "cblas_dsymv", dynlib: libName.}
proc dsbmv*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; k: Mkl_Int;
           alpha: cdouble; a: ptr cdouble; lda: Mkl_Int; x: ptr cdouble; incX: Mkl_Int;
           beta: cdouble; y: ptr cdouble; incY: Mkl_Int) {.importc: "cblas_dsbmv",
    dynlib: libName.}
proc dspmv*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: cdouble;
           ap: ptr cdouble; x: ptr cdouble; incX: Mkl_Int; beta: cdouble; y: ptr cdouble;
           incY: Mkl_Int) {.importc: "cblas_dspmv", dynlib: libName.}
proc dger*(layout: Cblas_Layout; m: Mkl_Int; n: Mkl_Int; alpha: cdouble; x: ptr cdouble;
          incX: Mkl_Int; y: ptr cdouble; incY: Mkl_Int; a: ptr cdouble; lda: Mkl_Int) {.
    importc: "cblas_dger", dynlib: libName.}
proc dsyr*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: cdouble;
          x: ptr cdouble; incX: Mkl_Int; a: ptr cdouble; lda: Mkl_Int) {.
    importc: "cblas_dsyr", dynlib: libName.}
proc dspr*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: cdouble;
          x: ptr cdouble; incX: Mkl_Int; ap: ptr cdouble) {.importc: "cblas_dspr",
    dynlib: libName.}
proc dsyr2*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: cdouble;
           x: ptr cdouble; incX: Mkl_Int; y: ptr cdouble; incY: Mkl_Int; a: ptr cdouble;
           lda: Mkl_Int) {.importc: "cblas_dsyr2", dynlib: libName.}
proc dspr2*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: cdouble;
           x: ptr cdouble; incX: Mkl_Int; y: ptr cdouble; incY: Mkl_Int; a: ptr cdouble) {.
    importc: "cblas_dspr2", dynlib: libName.}

##
##  Routines with C and Z prefixes only
##

proc chemv*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: pointer;
           a: pointer; lda: Mkl_Int; x: pointer; incX: Mkl_Int; beta: pointer; y: pointer;
           incY: Mkl_Int) {.importc: "cblas_chemv", dynlib: libName.}
proc chbmv*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; k: Mkl_Int;
           alpha: pointer; a: pointer; lda: Mkl_Int; x: pointer; incX: Mkl_Int;
           beta: pointer; y: pointer; incY: Mkl_Int) {.importc: "cblas_chbmv",
    dynlib: libName.}
proc chpmv*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: pointer;
           ap: pointer; x: pointer; incX: Mkl_Int; beta: pointer; y: pointer;
           incY: Mkl_Int) {.importc: "cblas_chpmv", dynlib: libName.}
proc cgeru*(layout: Cblas_Layout; m: Mkl_Int; n: Mkl_Int; alpha: pointer; x: pointer;
           incX: Mkl_Int; y: pointer; incY: Mkl_Int; a: pointer; lda: Mkl_Int) {.
    importc: "cblas_cgeru", dynlib: libName.}
proc cgerc*(layout: Cblas_Layout; m: Mkl_Int; n: Mkl_Int; alpha: pointer; x: pointer;
           incX: Mkl_Int; y: pointer; incY: Mkl_Int; a: pointer; lda: Mkl_Int) {.
    importc: "cblas_cgerc", dynlib: libName.}
proc cher*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: cfloat; x: pointer;
          incX: Mkl_Int; a: pointer; lda: Mkl_Int) {.importc: "cblas_cher",
    dynlib: libName.}
proc chpr*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: cfloat; x: pointer;
          incX: Mkl_Int; a: pointer) {.importc: "cblas_chpr", dynlib: libName.}
proc cher2*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: pointer;
           x: pointer; incX: Mkl_Int; y: pointer; incY: Mkl_Int; a: pointer; lda: Mkl_Int) {.
    importc: "cblas_cher2", dynlib: libName.}
proc chpr2*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: pointer;
           x: pointer; incX: Mkl_Int; y: pointer; incY: Mkl_Int; ap: pointer) {.
    importc: "cblas_chpr2", dynlib: libName.}
proc zhemv*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: pointer;
           a: pointer; lda: Mkl_Int; x: pointer; incX: Mkl_Int; beta: pointer; y: pointer;
           incY: Mkl_Int) {.importc: "cblas_zhemv", dynlib: libName.}
proc zhbmv*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; k: Mkl_Int;
           alpha: pointer; a: pointer; lda: Mkl_Int; x: pointer; incX: Mkl_Int;
           beta: pointer; y: pointer; incY: Mkl_Int) {.importc: "cblas_zhbmv",
    dynlib: libName.}
proc zhpmv*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: pointer;
           ap: pointer; x: pointer; incX: Mkl_Int; beta: pointer; y: pointer;
           incY: Mkl_Int) {.importc: "cblas_zhpmv", dynlib: libName.}
proc zgeru*(layout: Cblas_Layout; m: Mkl_Int; n: Mkl_Int; alpha: pointer; x: pointer;
           incX: Mkl_Int; y: pointer; incY: Mkl_Int; a: pointer; lda: Mkl_Int) {.
    importc: "cblas_zgeru", dynlib: libName.}
proc zgerc*(layout: Cblas_Layout; m: Mkl_Int; n: Mkl_Int; alpha: pointer; x: pointer;
           incX: Mkl_Int; y: pointer; incY: Mkl_Int; a: pointer; lda: Mkl_Int) {.
    importc: "cblas_zgerc", dynlib: libName.}
proc zher*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: cdouble;
          x: pointer; incX: Mkl_Int; a: pointer; lda: Mkl_Int) {.importc: "cblas_zher",
    dynlib: libName.}
proc zhpr*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: cdouble;
          x: pointer; incX: Mkl_Int; a: pointer) {.importc: "cblas_zhpr",
    dynlib: libName.}
proc zher2*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: pointer;
           x: pointer; incX: Mkl_Int; y: pointer; incY: Mkl_Int; a: pointer; lda: Mkl_Int) {.
    importc: "cblas_zher2", dynlib: libName.}
proc zhpr2*(layout: Cblas_Layout; uplo: Cblas_Uplo; n: Mkl_Int; alpha: pointer;
           x: pointer; incX: Mkl_Int; y: pointer; incY: Mkl_Int; ap: pointer) {.
    importc: "cblas_zhpr2", dynlib: libName.}

##
##  ===========================================================================
##  Prototypes for level 3 BLAS
##  ===========================================================================
##
##
##  Routines with standard 4 prefixes (S, D, C, Z)
##

proc sgemm*(layout: Cblas_Layout; transA: Cblas_Transpose; transB: Cblas_Transpose;
           m: Mkl_Int; n: Mkl_Int; k: Mkl_Int; alpha: cfloat; a: ptr cfloat; lda: Mkl_Int;
           b: ptr cfloat; ldb: Mkl_Int; beta: cfloat; c: ptr cfloat; ldc: Mkl_Int) {.
    importc: "cblas_sgemm", dynlib: libName.}
proc sgemmBatch*(layout: Cblas_Layout; transA_Array: ptr Cblas_Transpose;
                transB_Array: ptr Cblas_Transpose; m_Array: ptr Mkl_Int;
                n_Array: ptr Mkl_Int; k_Array: ptr Mkl_Int; alphaArray: ptr cfloat;
                a_Array: ptr ptr cfloat; ldaArray: ptr Mkl_Int;
                b_Array: ptr ptr cfloat; ldbArray: ptr Mkl_Int; betaArray: ptr cfloat;
                c_Array: ptr ptr cfloat; ldcArray: ptr Mkl_Int; groupCount: Mkl_Int;
                groupSize: ptr Mkl_Int) {.importc: "cblas_sgemm_batch",
                                       dynlib: libName.}
proc sgemmt*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
            transB: Cblas_Transpose; n: Mkl_Int; k: Mkl_Int; alpha: cfloat;
            a: ptr cfloat; lda: Mkl_Int; b: ptr cfloat; ldb: Mkl_Int; beta: cfloat;
            c: ptr cfloat; ldc: Mkl_Int) {.importc: "cblas_sgemmt", dynlib: libName.}
proc ssymm*(layout: Cblas_Layout; side: Cblas_Side; uplo: Cblas_Uplo; m: Mkl_Int;
           n: Mkl_Int; alpha: cfloat; a: ptr cfloat; lda: Mkl_Int; b: ptr cfloat;
           ldb: Mkl_Int; beta: cfloat; c: ptr cfloat; ldc: Mkl_Int) {.
    importc: "cblas_ssymm", dynlib: libName.}
proc ssyrk*(layout: Cblas_Layout; uplo: Cblas_Uplo; trans: Cblas_Transpose; n: Mkl_Int;
           k: Mkl_Int; alpha: cfloat; a: ptr cfloat; lda: Mkl_Int; beta: cfloat;
           c: ptr cfloat; ldc: Mkl_Int) {.importc: "cblas_ssyrk", dynlib: libName.}
proc ssyr2k*(layout: Cblas_Layout; uplo: Cblas_Uplo; trans: Cblas_Transpose;
            n: Mkl_Int; k: Mkl_Int; alpha: cfloat; a: ptr cfloat; lda: Mkl_Int;
            b: ptr cfloat; ldb: Mkl_Int; beta: cfloat; c: ptr cfloat; ldc: Mkl_Int) {.
    importc: "cblas_ssyr2k", dynlib: libName.}
proc strmm*(layout: Cblas_Layout; side: Cblas_Side; uplo: Cblas_Uplo;
           transA: Cblas_Transpose; diag: Cblas_Diag; m: Mkl_Int; n: Mkl_Int;
           alpha: cfloat; a: ptr cfloat; lda: Mkl_Int; b: ptr cfloat; ldb: Mkl_Int) {.
    importc: "cblas_strmm", dynlib: libName.}
proc strsm*(layout: Cblas_Layout; side: Cblas_Side; uplo: Cblas_Uplo;
           transA: Cblas_Transpose; diag: Cblas_Diag; m: Mkl_Int; n: Mkl_Int;
           alpha: cfloat; a: ptr cfloat; lda: Mkl_Int; b: ptr cfloat; ldb: Mkl_Int) {.
    importc: "cblas_strsm", dynlib: libName.}
proc strsmBatch*(layout: Cblas_Layout; sideArray: ptr Cblas_Side;
                uploArray: ptr Cblas_Uplo; transA_Array: ptr Cblas_Transpose;
                diagArray: ptr Cblas_Diag; m_Array: ptr Mkl_Int; n_Array: ptr Mkl_Int;
                alphaArray: ptr cfloat; a_Array: ptr ptr cfloat; ldaArray: ptr Mkl_Int;
                b_Array: ptr ptr cfloat; ldbArray: ptr Mkl_Int; groupCount: Mkl_Int;
                groupSize: ptr Mkl_Int) {.importc: "cblas_strsm_batch",
                                       dynlib: libName.}
proc dgemm*(layout: Cblas_Layout; transA: Cblas_Transpose; transB: Cblas_Transpose;
           m: Mkl_Int; n: Mkl_Int; k: Mkl_Int; alpha: cdouble; a: ptr cdouble;
           lda: Mkl_Int; b: ptr cdouble; ldb: Mkl_Int; beta: cdouble; c: ptr cdouble;
           ldc: Mkl_Int) {.importc: "cblas_dgemm", dynlib: libName.}
proc dgemmBatch*(layout: Cblas_Layout; transA_Array: ptr Cblas_Transpose;
                transB_Array: ptr Cblas_Transpose; m_Array: ptr Mkl_Int;
                n_Array: ptr Mkl_Int; k_Array: ptr Mkl_Int; alphaArray: ptr cdouble;
                a_Array: ptr ptr cdouble; ldaArray: ptr Mkl_Int;
                b_Array: ptr ptr cdouble; ldbArray: ptr Mkl_Int;
                betaArray: ptr cdouble; c_Array: ptr ptr cdouble;
                ldcArray: ptr Mkl_Int; groupCount: Mkl_Int; groupSize: ptr Mkl_Int) {.
    importc: "cblas_dgemm_batch", dynlib: libName.}
proc dgemmt*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
            transB: Cblas_Transpose; n: Mkl_Int; k: Mkl_Int; alpha: cdouble;
            a: ptr cdouble; lda: Mkl_Int; b: ptr cdouble; ldb: Mkl_Int; beta: cdouble;
            c: ptr cdouble; ldc: Mkl_Int) {.importc: "cblas_dgemmt", dynlib: libName.}
proc dsymm*(layout: Cblas_Layout; side: Cblas_Side; uplo: Cblas_Uplo; m: Mkl_Int;
           n: Mkl_Int; alpha: cdouble; a: ptr cdouble; lda: Mkl_Int; b: ptr cdouble;
           ldb: Mkl_Int; beta: cdouble; c: ptr cdouble; ldc: Mkl_Int) {.
    importc: "cblas_dsymm", dynlib: libName.}
proc dsyrk*(layout: Cblas_Layout; uplo: Cblas_Uplo; trans: Cblas_Transpose; n: Mkl_Int;
           k: Mkl_Int; alpha: cdouble; a: ptr cdouble; lda: Mkl_Int; beta: cdouble;
           c: ptr cdouble; ldc: Mkl_Int) {.importc: "cblas_dsyrk", dynlib: libName.}
proc dsyr2k*(layout: Cblas_Layout; uplo: Cblas_Uplo; trans: Cblas_Transpose;
            n: Mkl_Int; k: Mkl_Int; alpha: cdouble; a: ptr cdouble; lda: Mkl_Int;
            b: ptr cdouble; ldb: Mkl_Int; beta: cdouble; c: ptr cdouble; ldc: Mkl_Int) {.
    importc: "cblas_dsyr2k", dynlib: libName.}
proc dtrmm*(layout: Cblas_Layout; side: Cblas_Side; uplo: Cblas_Uplo;
           transA: Cblas_Transpose; diag: Cblas_Diag; m: Mkl_Int; n: Mkl_Int;
           alpha: cdouble; a: ptr cdouble; lda: Mkl_Int; b: ptr cdouble; ldb: Mkl_Int) {.
    importc: "cblas_dtrmm", dynlib: libName.}
proc dtrsm*(layout: Cblas_Layout; side: Cblas_Side; uplo: Cblas_Uplo;
           transA: Cblas_Transpose; diag: Cblas_Diag; m: Mkl_Int; n: Mkl_Int;
           alpha: cdouble; a: ptr cdouble; lda: Mkl_Int; b: ptr cdouble; ldb: Mkl_Int) {.
    importc: "cblas_dtrsm", dynlib: libName.}
proc dtrsmBatch*(layout: Cblas_Layout; sideArray: ptr Cblas_Side;
                uploArray: ptr Cblas_Uplo; transaArray: ptr Cblas_Transpose;
                diagArray: ptr Cblas_Diag; m_Array: ptr Mkl_Int; n_Array: ptr Mkl_Int;
                alphaArray: ptr cdouble; a_Array: ptr ptr cdouble;
                ldaArray: ptr Mkl_Int; b_Array: ptr ptr cdouble; ldbArray: ptr Mkl_Int;
                groupCount: Mkl_Int; groupSize: ptr Mkl_Int) {.
    importc: "cblas_dtrsm_batch", dynlib: libName.}
proc cgemm*(layout: Cblas_Layout; transA: Cblas_Transpose; transB: Cblas_Transpose;
           m: Mkl_Int; n: Mkl_Int; k: Mkl_Int; alpha: pointer; a: pointer; lda: Mkl_Int;
           b: pointer; ldb: Mkl_Int; beta: pointer; c: pointer; ldc: Mkl_Int) {.
    importc: "cblas_cgemm", dynlib: libName.}
proc cgemm3m*(layout: Cblas_Layout; transA: Cblas_Transpose; transB: Cblas_Transpose;
             m: Mkl_Int; n: Mkl_Int; k: Mkl_Int; alpha: pointer; a: pointer; lda: Mkl_Int;
             b: pointer; ldb: Mkl_Int; beta: pointer; c: pointer; ldc: Mkl_Int) {.
    importc: "cblas_cgemm3m", dynlib: libName.}
proc cgemmBatch*(layout: Cblas_Layout; transA_Array: ptr Cblas_Transpose;
                transB_Array: ptr Cblas_Transpose; m_Array: ptr Mkl_Int;
                n_Array: ptr Mkl_Int; k_Array: ptr Mkl_Int; alphaArray: pointer;
                a_Array: ptr pointer; ldaArray: ptr Mkl_Int; b_Array: ptr pointer;
                ldbArray: ptr Mkl_Int; betaArray: pointer; c_Array: ptr pointer;
                ldcArray: ptr Mkl_Int; groupCount: Mkl_Int; groupSize: ptr Mkl_Int) {.
    importc: "cblas_cgemm_batch", dynlib: libName.}
proc cgemm3mBatch*(layout: Cblas_Layout; transA_Array: ptr Cblas_Transpose;
                  transB_Array: ptr Cblas_Transpose; m_Array: ptr Mkl_Int;
                  n_Array: ptr Mkl_Int; k_Array: ptr Mkl_Int; alphaArray: pointer;
                  a_Array: ptr pointer; ldaArray: ptr Mkl_Int; b_Array: ptr pointer;
                  ldbArray: ptr Mkl_Int; betaArray: pointer; c_Array: ptr pointer;
                  ldcArray: ptr Mkl_Int; groupCount: Mkl_Int; groupSize: ptr Mkl_Int) {.
    importc: "cblas_cgemm3m_batch", dynlib: libName.}
proc cgemmt*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
            transB: Cblas_Transpose; n: Mkl_Int; k: Mkl_Int; alpha: pointer; a: pointer;
            lda: Mkl_Int; b: pointer; ldb: Mkl_Int; beta: pointer; c: pointer;
            ldc: Mkl_Int) {.importc: "cblas_cgemmt", dynlib: libName.}
proc csymm*(layout: Cblas_Layout; side: Cblas_Side; uplo: Cblas_Uplo; m: Mkl_Int;
           n: Mkl_Int; alpha: pointer; a: pointer; lda: Mkl_Int; b: pointer; ldb: Mkl_Int;
           beta: pointer; c: pointer; ldc: Mkl_Int) {.importc: "cblas_csymm",
    dynlib: libName.}
proc csyrk*(layout: Cblas_Layout; uplo: Cblas_Uplo; trans: Cblas_Transpose; n: Mkl_Int;
           k: Mkl_Int; alpha: pointer; a: pointer; lda: Mkl_Int; beta: pointer;
           c: pointer; ldc: Mkl_Int) {.importc: "cblas_csyrk", dynlib: libName.}
proc csyr2k*(layout: Cblas_Layout; uplo: Cblas_Uplo; trans: Cblas_Transpose;
            n: Mkl_Int; k: Mkl_Int; alpha: pointer; a: pointer; lda: Mkl_Int; b: pointer;
            ldb: Mkl_Int; beta: pointer; c: pointer; ldc: Mkl_Int) {.
    importc: "cblas_csyr2k", dynlib: libName.}
proc ctrmm*(layout: Cblas_Layout; side: Cblas_Side; uplo: Cblas_Uplo;
           transA: Cblas_Transpose; diag: Cblas_Diag; m: Mkl_Int; n: Mkl_Int;
           alpha: pointer; a: pointer; lda: Mkl_Int; b: pointer; ldb: Mkl_Int) {.
    importc: "cblas_ctrmm", dynlib: libName.}
proc ctrsm*(layout: Cblas_Layout; side: Cblas_Side; uplo: Cblas_Uplo;
           transA: Cblas_Transpose; diag: Cblas_Diag; m: Mkl_Int; n: Mkl_Int;
           alpha: pointer; a: pointer; lda: Mkl_Int; b: pointer; ldb: Mkl_Int) {.
    importc: "cblas_ctrsm", dynlib: libName.}
proc ctrsmBatch*(layout: Cblas_Layout; sideArray: ptr Cblas_Side;
                uploArray: ptr Cblas_Uplo; transaArray: ptr Cblas_Transpose;
                diagArray: ptr Cblas_Diag; m_Array: ptr Mkl_Int; n_Array: ptr Mkl_Int;
                alphaArray: pointer; a_Array: ptr pointer; ldaArray: ptr Mkl_Int;
                b_Array: ptr pointer; ldbArray: ptr Mkl_Int; groupCount: Mkl_Int;
                groupSize: ptr Mkl_Int) {.importc: "cblas_ctrsm_batch",
                                       dynlib: libName.}
proc zgemm*(layout: Cblas_Layout; transA: Cblas_Transpose; transB: Cblas_Transpose;
           m: Mkl_Int; n: Mkl_Int; k: Mkl_Int; alpha: pointer; a: pointer; lda: Mkl_Int;
           b: pointer; ldb: Mkl_Int; beta: pointer; c: pointer; ldc: Mkl_Int) {.
    importc: "cblas_zgemm", dynlib: libName.}
proc zgemm3m*(layout: Cblas_Layout; transA: Cblas_Transpose; transB: Cblas_Transpose;
             m: Mkl_Int; n: Mkl_Int; k: Mkl_Int; alpha: pointer; a: pointer; lda: Mkl_Int;
             b: pointer; ldb: Mkl_Int; beta: pointer; c: pointer; ldc: Mkl_Int) {.
    importc: "cblas_zgemm3m", dynlib: libName.}
proc zgemmBatch*(layout: Cblas_Layout; transA_Array: ptr Cblas_Transpose;
                transB_Array: ptr Cblas_Transpose; m_Array: ptr Mkl_Int;
                n_Array: ptr Mkl_Int; k_Array: ptr Mkl_Int; alphaArray: pointer;
                a_Array: ptr pointer; ldaArray: ptr Mkl_Int; b_Array: ptr pointer;
                ldbArray: ptr Mkl_Int; betaArray: pointer; c_Array: ptr pointer;
                ldcArray: ptr Mkl_Int; groupCount: Mkl_Int; groupSize: ptr Mkl_Int) {.
    importc: "cblas_zgemm_batch", dynlib: libName.}
proc zgemm3mBatch*(layout: Cblas_Layout; transA_Array: ptr Cblas_Transpose;
                  transB_Array: ptr Cblas_Transpose; m_Array: ptr Mkl_Int;
                  n_Array: ptr Mkl_Int; k_Array: ptr Mkl_Int; alphaArray: pointer;
                  a_Array: ptr pointer; ldaArray: ptr Mkl_Int; b_Array: ptr pointer;
                  ldbArray: ptr Mkl_Int; betaArray: pointer; c_Array: ptr pointer;
                  ldcArray: ptr Mkl_Int; groupCount: Mkl_Int; groupSize: ptr Mkl_Int) {.
    importc: "cblas_zgemm3m_batch", dynlib: libName.}
proc zgemmt*(layout: Cblas_Layout; uplo: Cblas_Uplo; transA: Cblas_Transpose;
            transB: Cblas_Transpose; n: Mkl_Int; k: Mkl_Int; alpha: pointer; a: pointer;
            lda: Mkl_Int; b: pointer; ldb: Mkl_Int; beta: pointer; c: pointer;
            ldc: Mkl_Int) {.importc: "cblas_zgemmt", dynlib: libName.}
proc zsymm*(layout: Cblas_Layout; side: Cblas_Side; uplo: Cblas_Uplo; m: Mkl_Int;
           n: Mkl_Int; alpha: pointer; a: pointer; lda: Mkl_Int; b: pointer; ldb: Mkl_Int;
           beta: pointer; c: pointer; ldc: Mkl_Int) {.importc: "cblas_zsymm",
    dynlib: libName.}
proc zsyrk*(layout: Cblas_Layout; uplo: Cblas_Uplo; trans: Cblas_Transpose; n: Mkl_Int;
           k: Mkl_Int; alpha: pointer; a: pointer; lda: Mkl_Int; beta: pointer;
           c: pointer; ldc: Mkl_Int) {.importc: "cblas_zsyrk", dynlib: libName.}
proc zsyr2k*(layout: Cblas_Layout; uplo: Cblas_Uplo; trans: Cblas_Transpose;
            n: Mkl_Int; k: Mkl_Int; alpha: pointer; a: pointer; lda: Mkl_Int; b: pointer;
            ldb: Mkl_Int; beta: pointer; c: pointer; ldc: Mkl_Int) {.
    importc: "cblas_zsyr2k", dynlib: libName.}
proc ztrmm*(layout: Cblas_Layout; side: Cblas_Side; uplo: Cblas_Uplo;
           transA: Cblas_Transpose; diag: Cblas_Diag; m: Mkl_Int; n: Mkl_Int;
           alpha: pointer; a: pointer; lda: Mkl_Int; b: pointer; ldb: Mkl_Int) {.
    importc: "cblas_ztrmm", dynlib: libName.}
proc ztrsm*(layout: Cblas_Layout; side: Cblas_Side; uplo: Cblas_Uplo;
           transA: Cblas_Transpose; diag: Cblas_Diag; m: Mkl_Int; n: Mkl_Int;
           alpha: pointer; a: pointer; lda: Mkl_Int; b: pointer; ldb: Mkl_Int) {.
    importc: "cblas_ztrsm", dynlib: libName.}
proc ztrsmBatch*(layout: Cblas_Layout; sideArray: ptr Cblas_Side;
                uploArray: ptr Cblas_Uplo; transaArray: ptr Cblas_Transpose;
                diagArray: ptr Cblas_Diag; m_Array: ptr Mkl_Int; n_Array: ptr Mkl_Int;
                alphaArray: pointer; a_Array: ptr pointer; ldaArray: ptr Mkl_Int;
                b_Array: ptr pointer; ldbArray: ptr Mkl_Int; groupCount: Mkl_Int;
                groupSize: ptr Mkl_Int) {.importc: "cblas_ztrsm_batch",
                                       dynlib: libName.}

##
##  Routines with prefixes C and Z only
##

proc chemm*(layout: Cblas_Layout; side: Cblas_Side; uplo: Cblas_Uplo; m: Mkl_Int;
           n: Mkl_Int; alpha: pointer; a: pointer; lda: Mkl_Int; b: pointer; ldb: Mkl_Int;
           beta: pointer; c: pointer; ldc: Mkl_Int) {.importc: "cblas_chemm",
    dynlib: libName.}
proc cherk*(layout: Cblas_Layout; uplo: Cblas_Uplo; trans: Cblas_Transpose; n: Mkl_Int;
           k: Mkl_Int; alpha: cfloat; a: pointer; lda: Mkl_Int; beta: cfloat; c: pointer;
           ldc: Mkl_Int) {.importc: "cblas_cherk", dynlib: libName.}
proc cher2k*(layout: Cblas_Layout; uplo: Cblas_Uplo; trans: Cblas_Transpose;
            n: Mkl_Int; k: Mkl_Int; alpha: pointer; a: pointer; lda: Mkl_Int; b: pointer;
            ldb: Mkl_Int; beta: cfloat; c: pointer; ldc: Mkl_Int) {.
    importc: "cblas_cher2k", dynlib: libName.}
proc zhemm*(layout: Cblas_Layout; side: Cblas_Side; uplo: Cblas_Uplo; m: Mkl_Int;
           n: Mkl_Int; alpha: pointer; a: pointer; lda: Mkl_Int; b: pointer; ldb: Mkl_Int;
           beta: pointer; c: pointer; ldc: Mkl_Int) {.importc: "cblas_zhemm",
    dynlib: libName.}
proc zherk*(layout: Cblas_Layout; uplo: Cblas_Uplo; trans: Cblas_Transpose; n: Mkl_Int;
           k: Mkl_Int; alpha: cdouble; a: pointer; lda: Mkl_Int; beta: cdouble;
           c: pointer; ldc: Mkl_Int) {.importc: "cblas_zherk", dynlib: libName.}
proc zher2k*(layout: Cblas_Layout; uplo: Cblas_Uplo; trans: Cblas_Transpose;
            n: Mkl_Int; k: Mkl_Int; alpha: pointer; a: pointer; lda: Mkl_Int; b: pointer;
            ldb: Mkl_Int; beta: cdouble; c: pointer; ldc: Mkl_Int) {.
    importc: "cblas_zher2k", dynlib: libName.}

##
##  Routines with prefixes S and D only
##

proc sgemmAlloc*(identifier: Cblas_Identifier; m: Mkl_Int; n: Mkl_Int; k: Mkl_Int): ptr cfloat {.
    importc: "cblas_sgemm_alloc", dynlib: libName.}
proc sgemmPackGetSize*(identifier: Cblas_Identifier; m: Mkl_Int; n: Mkl_Int; k: Mkl_Int): Mkl_Int {.
    importc: "cblas_sgemm_pack_get_size", dynlib: libName.}
proc sgemmPack*(layout: Cblas_Layout; identifier: Cblas_Identifier;
               trans: Cblas_Transpose; m: Mkl_Int; n: Mkl_Int; k: Mkl_Int;
               alpha: cfloat; src: ptr cfloat; ld: Mkl_Int; dest: ptr cfloat) {.
    importc: "cblas_sgemm_pack", dynlib: libName.}
proc sgemmCompute*(layout: Cblas_Layout; transA: Mkl_Int; transB: Mkl_Int; m: Mkl_Int;
                  n: Mkl_Int; k: Mkl_Int; a: ptr cfloat; lda: Mkl_Int; b: ptr cfloat;
                  ldb: Mkl_Int; beta: cfloat; c: ptr cfloat; ldc: Mkl_Int) {.
    importc: "cblas_sgemm_compute", dynlib: libName.}
proc sgemmFree*(dest: ptr cfloat) {.importc: "cblas_sgemm_free", dynlib: libName.}
proc dgemmAlloc*(identifier: Cblas_Identifier; m: Mkl_Int; n: Mkl_Int; k: Mkl_Int): ptr cdouble {.
    importc: "cblas_dgemm_alloc", dynlib: libName.}
proc dgemmPackGetSize*(identifier: Cblas_Identifier; m: Mkl_Int; n: Mkl_Int; k: Mkl_Int): Mkl_Int {.
    importc: "cblas_dgemm_pack_get_size", dynlib: libName.}
proc dgemmPack*(layout: Cblas_Layout; identifier: Cblas_Identifier;
               trans: Cblas_Transpose; m: Mkl_Int; n: Mkl_Int; k: Mkl_Int;
               alpha: cdouble; src: ptr cdouble; ld: Mkl_Int; dest: ptr cdouble) {.
    importc: "cblas_dgemm_pack", dynlib: libName.}
proc dgemmCompute*(layout: Cblas_Layout; transA: Mkl_Int; transB: Mkl_Int; m: Mkl_Int;
                  n: Mkl_Int; k: Mkl_Int; a: ptr cdouble; lda: Mkl_Int; b: ptr cdouble;
                  ldb: Mkl_Int; beta: cdouble; c: ptr cdouble; ldc: Mkl_Int) {.
    importc: "cblas_dgemm_compute", dynlib: libName.}
proc dgemmFree*(dest: ptr cdouble) {.importc: "cblas_dgemm_free", dynlib: libName.}

##
##  Integer Routines
##

proc gemmS16s16s32*(layout: Cblas_Layout; transA: Cblas_Transpose;
                   transB: Cblas_Transpose; offsetC: Cblas_Offset; m: Mkl_Int;
                   n: Mkl_Int; k: Mkl_Int; alpha: cfloat; a: ptr Mkl_Int16; lda: Mkl_Int;
                   ao: Mkl_Int16; b: ptr Mkl_Int16; ldb: Mkl_Int; bo: Mkl_Int16;
                   beta: cfloat; c: ptr Mkl_Int32; ldc: Mkl_Int; cb: ptr Mkl_Int32) {.
    importc: "cblas_gemm_s16s16s32", dynlib: libName.}
proc gemmS8u8s32*(layout: Cblas_Layout; transA: Cblas_Transpose;
                 transB: Cblas_Transpose; offsetC: Cblas_Offset; m: Mkl_Int;
                 n: Mkl_Int; k: Mkl_Int; alpha: cfloat; a: pointer; lda: Mkl_Int;
                 ao: Mkl_Int8; b: pointer; ldb: Mkl_Int; bo: Mkl_Int8; beta: cfloat;
                 c: ptr Mkl_Int32; ldc: Mkl_Int; cb: ptr Mkl_Int32) {.
    importc: "cblas_gemm_s8u8s32", dynlib: libName.}
proc gemmS8u8s32PackGetSize*(identifier: Cblas_Identifier; m: Mkl_Int; n: Mkl_Int;
                            k: Mkl_Int): Mkl_Int {.
    importc: "cblas_gemm_s8u8s32_pack_get_size", dynlib: libName.}
proc gemmS16s16s32PackGetSize*(identifier: Cblas_Identifier; m: Mkl_Int; n: Mkl_Int;
                              k: Mkl_Int): Mkl_Int {.
    importc: "cblas_gemm_s16s16s32_pack_get_size", dynlib: libName.}
proc gemmS8u8s32Pack*(layout: Cblas_Layout; identifier: Cblas_Identifier;
                     trans: Cblas_Transpose; m: Mkl_Int; n: Mkl_Int; k: Mkl_Int;
                     src: pointer; ld: Mkl_Int; dest: pointer) {.
    importc: "cblas_gemm_s8u8s32_pack", dynlib: libName.}
proc gemmS16s16s32Pack*(layout: Cblas_Layout; identifier: Cblas_Identifier;
                       trans: Cblas_Transpose; m: Mkl_Int; n: Mkl_Int; k: Mkl_Int;
                       src: ptr Mkl_Int16; ld: Mkl_Int; dest: ptr Mkl_Int16) {.
    importc: "cblas_gemm_s16s16s32_pack", dynlib: libName.}
proc gemmS8u8s32Compute*(layout: Cblas_Layout; transA: Cblas_Transpose;
                        transB: Cblas_Transpose; offsetc: Cblas_Offset; m: Mkl_Int;
                        n: Mkl_Int; k: Mkl_Int; alpha: cfloat; a: pointer; lda: Mkl_Int;
                        ao: Mkl_Int8; b: pointer; ldb: Mkl_Int; bo: Mkl_Int8;
                        beta: cfloat; c: ptr Mkl_Int32; ldc: Mkl_Int; co: ptr Mkl_Int32) {.
    importc: "cblas_gemm_s8u8s32_compute", dynlib: libName.}
proc gemmS16s16s32Compute*(layout: Cblas_Layout; transA: Cblas_Transpose;
                          transB: Cblas_Transpose; offsetc: Cblas_Offset;
                          m: Mkl_Int; n: Mkl_Int; k: Mkl_Int; alpha: cfloat;
                          a: ptr Mkl_Int16; lda: Mkl_Int; ao: Mkl_Int16;
                          b: ptr Mkl_Int16; ldb: Mkl_Int; bo: Mkl_Int16; beta: cfloat;
                          c: ptr Mkl_Int32; ldc: Mkl_Int; co: ptr Mkl_Int32) {.
    importc: "cblas_gemm_s16s16s32_compute", dynlib: libName.}

