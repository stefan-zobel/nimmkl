import mklTypes

const
  LAPACK_ROW_MAJOR* = 101
  LAPACK_COL_MAJOR* = 102
  LAPACK_WORK_MEMORY_ERROR* = -1010
  LAPACK_TRANSPOSE_MEMORY_ERROR* = -1011

##  Callback logical functions of one, two, or three arguments are used
##  to select eigenvalues to sort to the top left of the Schur form.
##  The value is selected if function returns non-zero.

type
  Lapack_S_Select2* = proc (a1: ptr cfloat; a2: ptr cfloat): LapackLogical {.cdecl.}
  Lapack_S_Select3* = proc (a1: ptr cfloat; a2: ptr cfloat; a3: ptr cfloat): LapackLogical {.cdecl.}
  Lapack_D_Select2* = proc (a1: ptr cdouble; a2: ptr cdouble): LapackLogical {.cdecl.}
  Lapack_D_Select3* = proc (a1: ptr cdouble; a2: ptr cdouble; a3: ptr cdouble): LapackLogical {.cdecl.}
  Lapack_C_Select1* = proc (a1: ptr LapackComplexFloat): LapackLogical {.cdecl.}
  Lapack_C_Select2* = proc (a1: ptr LapackComplexFloat; a2: ptr LapackComplexFloat): LapackLogical {.cdecl.}
  Lapack_Z_Select1* = proc (a1: ptr LapackComplexDouble): LapackLogical {.cdecl.}
  Lapack_Z_Select2* = proc (a1: ptr LapackComplexDouble; a2: ptr LapackComplexDouble): LapackLogical {.cdecl.}


proc dlamch*(cmach: char): cdouble {.importc: "LAPACKE_dlamch", dynlib: libName.}
proc dlamchWork*(cmach: char): cdouble {.importc: "LAPACKE_dlamch_work",
                                     dynlib: libName.}
proc dlange*(matrixLayout: cint; norm: char; m: LapackInt; n: LapackInt; a: ptr cdouble;
            lda: LapackInt): cdouble {.importc: "LAPACKE_dlange", dynlib: libName.}
proc dlangeWork*(matrixLayout: cint; norm: char; m: LapackInt; n: LapackInt;
                a: ptr cdouble; lda: LapackInt; work: ptr cdouble): cdouble {.
    importc: "LAPACKE_dlange_work", dynlib: libName.}
proc dlansy*(matrixLayout: cint; norm: char; uplo: char; n: LapackInt; a: ptr cdouble;
            lda: LapackInt): cdouble {.importc: "LAPACKE_dlansy", dynlib: libName.}
proc dlansyWork*(matrixLayout: cint; norm: char; uplo: char; n: LapackInt;
                a: ptr cdouble; lda: LapackInt; work: ptr cdouble): cdouble {.
    importc: "LAPACKE_dlansy_work", dynlib: libName.}
proc dlantr*(matrixLayout: cint; norm: char; uplo: char; diag: char; m: LapackInt;
            n: LapackInt; a: ptr cdouble; lda: LapackInt): cdouble {.
    importc: "LAPACKE_dlantr", dynlib: libName.}
proc dlantrWork*(matrixLayout: cint; norm: char; uplo: char; diag: char; m: LapackInt;
                n: LapackInt; a: ptr cdouble; lda: LapackInt; work: ptr cdouble): cdouble {.
    importc: "LAPACKE_dlantr_work", dynlib: libName.}
proc dlapy2*(x: cdouble; y: cdouble): cdouble {.importc: "LAPACKE_dlapy2",
    dynlib: libName.}
proc dlapy2Work*(x: cdouble; y: cdouble): cdouble {.importc: "LAPACKE_dlapy2_work",
    dynlib: libName.}
proc dlapy3*(x: cdouble; y: cdouble; z: cdouble): cdouble {.importc: "LAPACKE_dlapy3",
    dynlib: libName.}
proc dlapy3Work*(x: cdouble; y: cdouble; z: cdouble): cdouble {.
    importc: "LAPACKE_dlapy3_work", dynlib: libName.}
proc zlange*(matrixLayout: cint; norm: char; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt): cdouble {.
    importc: "LAPACKE_zlange", dynlib: libName.}
proc zlangeWork*(matrixLayout: cint; norm: char; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; work: ptr cdouble): cdouble {.
    importc: "LAPACKE_zlange_work", dynlib: libName.}
proc zlanhe*(matrixLayout: cint; norm: char; uplo: char; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt): cdouble {.
    importc: "LAPACKE_zlanhe", dynlib: libName.}
proc zlanheWork*(matrixLayout: cint; norm: char; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; work: ptr cdouble): cdouble {.
    importc: "LAPACKE_zlanhe_work", dynlib: libName.}
proc zlarcm*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
            lda: LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt;
            c: ptr LapackComplexDouble; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_zlarcm", dynlib: libName.}
proc zlarcmWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt;
                c: ptr LapackComplexDouble; ldc: LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zlarcm_work", dynlib: libName.}
proc zlacrm*(matrixLayout: cint; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
            c: ptr LapackComplexDouble; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_zlacrm", dynlib: libName.}
proc zlacrmWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; b: ptr cdouble;
                ldb: LapackInt; c: ptr LapackComplexDouble; ldc: LapackInt;
                work: ptr cdouble): LapackInt {.importc: "LAPACKE_zlacrm_work",
    dynlib: libName.}
proc zlansy*(matrixLayout: cint; norm: char; uplo: char; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt): cdouble {.
    importc: "LAPACKE_zlansy", dynlib: libName.}
proc zlansyWork*(matrixLayout: cint; norm: char; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; work: ptr cdouble): cdouble {.
    importc: "LAPACKE_zlansy_work", dynlib: libName.}
proc zlantr*(matrixLayout: cint; norm: char; uplo: char; diag: char; m: LapackInt;
            n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt): cdouble {.
    importc: "LAPACKE_zlantr", dynlib: libName.}
proc zlantrWork*(matrixLayout: cint; norm: char; uplo: char; diag: char; m: LapackInt;
                n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                work: ptr cdouble): cdouble {.importc: "LAPACKE_zlantr_work",
    dynlib: libName.}
proc clange*(matrixLayout: cint; norm: char; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt): cfloat {.
    importc: "LAPACKE_clange", dynlib: libName.}
proc clangeWork*(matrixLayout: cint; norm: char; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; work: ptr cfloat): cfloat {.
    importc: "LAPACKE_clange_work", dynlib: libName.}
proc clanhe*(matrixLayout: cint; norm: char; uplo: char; n: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt): cfloat {.
    importc: "LAPACKE_clanhe", dynlib: libName.}
proc clanheWork*(matrixLayout: cint; norm: char; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; work: ptr cfloat): cfloat {.
    importc: "LAPACKE_clanhe_work", dynlib: libName.}
proc clarcm*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
            lda: LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt;
            c: ptr LapackComplexFloat; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_clarcm", dynlib: libName.}
proc clarcmWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt;
                c: ptr LapackComplexFloat; ldc: LapackInt; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_clarcm_work", dynlib: libName.}
proc clacrm*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; b: ptr cfloat; ldb: LapackInt; c: ptr LapackComplexFloat;
            ldc: LapackInt): LapackInt {.importc: "LAPACKE_clacrm", dynlib: libName.}
proc clacrmWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; b: ptr cfloat;
                ldb: LapackInt; c: ptr LapackComplexFloat; ldc: LapackInt;
                rwork: ptr cfloat): LapackInt {.importc: "LAPACKE_clacrm_work",
    dynlib: libName.}
proc clansy*(matrixLayout: cint; norm: char; uplo: char; n: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt): cfloat {.
    importc: "LAPACKE_clansy", dynlib: libName.}
proc clansyWork*(matrixLayout: cint; norm: char; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; work: ptr cfloat): cfloat {.
    importc: "LAPACKE_clansy_work", dynlib: libName.}
proc clantr*(matrixLayout: cint; norm: char; uplo: char; diag: char; m: LapackInt;
            n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt): cfloat {.
    importc: "LAPACKE_clantr", dynlib: libName.}
proc clantrWork*(matrixLayout: cint; norm: char; uplo: char; diag: char; m: LapackInt;
                n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                work: ptr cfloat): cfloat {.importc: "LAPACKE_clantr_work",
                                        dynlib: libName.}
proc slamch*(cmach: char): cfloat {.importc: "LAPACKE_slamch", dynlib: libName.}
proc slamchWork*(cmach: char): cfloat {.importc: "LAPACKE_slamch_work",
                                    dynlib: libName.}
proc slange*(matrixLayout: cint; norm: char; m: LapackInt; n: LapackInt; a: ptr cfloat;
            lda: LapackInt): cfloat {.importc: "LAPACKE_slange", dynlib: libName.}
proc slangeWork*(matrixLayout: cint; norm: char; m: LapackInt; n: LapackInt;
                a: ptr cfloat; lda: LapackInt; work: ptr cfloat): cfloat {.
    importc: "LAPACKE_slange_work", dynlib: libName.}
proc slansy*(matrixLayout: cint; norm: char; uplo: char; n: LapackInt; a: ptr cfloat;
            lda: LapackInt): cfloat {.importc: "LAPACKE_slansy", dynlib: libName.}
proc slansyWork*(matrixLayout: cint; norm: char; uplo: char; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; work: ptr cfloat): cfloat {.
    importc: "LAPACKE_slansy_work", dynlib: libName.}
proc slantr*(matrixLayout: cint; norm: char; uplo: char; diag: char; m: LapackInt;
            n: LapackInt; a: ptr cfloat; lda: LapackInt): cfloat {.
    importc: "LAPACKE_slantr", dynlib: libName.}
proc slantrWork*(matrixLayout: cint; norm: char; uplo: char; diag: char; m: LapackInt;
                n: LapackInt; a: ptr cfloat; lda: LapackInt; work: ptr cfloat): cfloat {.
    importc: "LAPACKE_slantr_work", dynlib: libName.}
proc slapy2*(x: cfloat; y: cfloat): cfloat {.importc: "LAPACKE_slapy2", dynlib: libName.}
proc slapy2Work*(x: cfloat; y: cfloat): cfloat {.importc: "LAPACKE_slapy2_work",
    dynlib: libName.}
proc slapy3*(x: cfloat; y: cfloat; z: cfloat): cfloat {.importc: "LAPACKE_slapy3",
    dynlib: libName.}
proc slapy3Work*(x: cfloat; y: cfloat; z: cfloat): cfloat {.
    importc: "LAPACKE_slapy3_work", dynlib: libName.}
proc cbbcsd*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char; jobv2t: char;
            trans: char; m: LapackInt; p: LapackInt; q: LapackInt; theta: ptr cfloat;
            phi: ptr cfloat; u1: ptr LapackComplexFloat; ldu1: LapackInt;
            u2: ptr LapackComplexFloat; ldu2: LapackInt; v1t: ptr LapackComplexFloat;
            ldv1t: LapackInt; v2t: ptr LapackComplexFloat; ldv2t: LapackInt;
            b11d: ptr cfloat; b11e: ptr cfloat; b12d: ptr cfloat; b12e: ptr cfloat;
            b21d: ptr cfloat; b21e: ptr cfloat; b22d: ptr cfloat; b22e: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cbbcsd", dynlib: libName.}
proc cbbcsdWork*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char;
                jobv2t: char; trans: char; m: LapackInt; p: LapackInt; q: LapackInt;
                theta: ptr cfloat; phi: ptr cfloat; u1: ptr LapackComplexFloat;
                ldu1: LapackInt; u2: ptr LapackComplexFloat; ldu2: LapackInt;
                v1t: ptr LapackComplexFloat; ldv1t: LapackInt;
                v2t: ptr LapackComplexFloat; ldv2t: LapackInt; b11d: ptr cfloat;
                b11e: ptr cfloat; b12d: ptr cfloat; b12e: ptr cfloat; b21d: ptr cfloat;
                b21e: ptr cfloat; b22d: ptr cfloat; b22e: ptr cfloat; rwork: ptr cfloat;
                lrwork: LapackInt): LapackInt {.importc: "LAPACKE_cbbcsd_work",
    dynlib: libName.}
proc cbdsqr*(matrixLayout: cint; uplo: char; n: LapackInt; ncvt: LapackInt;
            nru: LapackInt; ncc: LapackInt; d: ptr cfloat; e: ptr cfloat;
            vt: ptr LapackComplexFloat; ldvt: LapackInt; u: ptr LapackComplexFloat;
            ldu: LapackInt; c: ptr LapackComplexFloat; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_cbdsqr", dynlib: libName.}
proc cbdsqrWork*(matrixLayout: cint; uplo: char; n: LapackInt; ncvt: LapackInt;
                nru: LapackInt; ncc: LapackInt; d: ptr cfloat; e: ptr cfloat;
                vt: ptr LapackComplexFloat; ldvt: LapackInt;
                u: ptr LapackComplexFloat; ldu: LapackInt; c: ptr LapackComplexFloat;
                ldc: LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cbdsqr_work", dynlib: libName.}
proc cgbbrd*(matrixLayout: cint; vect: char; m: LapackInt; n: LapackInt; ncc: LapackInt;
            kl: LapackInt; ku: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
            d: ptr cfloat; e: ptr cfloat; q: ptr LapackComplexFloat; ldq: LapackInt;
            pt: ptr LapackComplexFloat; ldpt: LapackInt; c: ptr LapackComplexFloat;
            ldc: LapackInt): LapackInt {.importc: "LAPACKE_cgbbrd", dynlib: libName.}
proc cgbbrdWork*(matrixLayout: cint; vect: char; m: LapackInt; n: LapackInt;
                ncc: LapackInt; kl: LapackInt; ku: LapackInt;
                ab: ptr LapackComplexFloat; ldab: LapackInt; d: ptr cfloat;
                e: ptr cfloat; q: ptr LapackComplexFloat; ldq: LapackInt;
                pt: ptr LapackComplexFloat; ldpt: LapackInt;
                c: ptr LapackComplexFloat; ldc: LapackInt;
                work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgbbrd_work", dynlib: libName.}
proc cgbcon*(matrixLayout: cint; norm: char; n: LapackInt; kl: LapackInt; ku: LapackInt;
            ab: ptr LapackComplexFloat; ldab: LapackInt; ipiv: ptr LapackInt;
            anorm: cfloat; rcond: ptr cfloat): LapackInt {.importc: "LAPACKE_cgbcon",
    dynlib: libName.}
proc cgbconWork*(matrixLayout: cint; norm: char; n: LapackInt; kl: LapackInt;
                ku: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
                ipiv: ptr LapackInt; anorm: cfloat; rcond: ptr cfloat;
                work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgbcon_work", dynlib: libName.}
proc cgbequ*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
            ku: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt; r: ptr cfloat;
            c: ptr cfloat; rowcnd: ptr cfloat; colcnd: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgbequ", dynlib: libName.}
proc cgbequWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
                ku: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
                r: ptr cfloat; c: ptr cfloat; rowcnd: ptr cfloat; colcnd: ptr cfloat;
                amax: ptr cfloat): LapackInt {.importc: "LAPACKE_cgbequ_work",
    dynlib: libName.}
proc cgbequb*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
             ku: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt; r: ptr cfloat;
             c: ptr cfloat; rowcnd: ptr cfloat; colcnd: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgbequb", dynlib: libName.}
proc cgbequbWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
                 ku: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
                 r: ptr cfloat; c: ptr cfloat; rowcnd: ptr cfloat; colcnd: ptr cfloat;
                 amax: ptr cfloat): LapackInt {.importc: "LAPACKE_cgbequb_work",
    dynlib: libName.}
proc cgbrfs*(matrixLayout: cint; trans: char; n: LapackInt; kl: LapackInt; ku: LapackInt;
            nrhs: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
            afb: ptr LapackComplexFloat; ldafb: LapackInt; ipiv: ptr LapackInt;
            b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
            ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgbrfs", dynlib: libName.}
proc cgbrfsWork*(matrixLayout: cint; trans: char; n: LapackInt; kl: LapackInt;
                ku: LapackInt; nrhs: LapackInt; ab: ptr LapackComplexFloat;
                ldab: LapackInt; afb: ptr LapackComplexFloat; ldafb: LapackInt;
                ipiv: ptr LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt;
                x: ptr LapackComplexFloat; ldx: LapackInt; ferr: ptr cfloat;
                berr: ptr cfloat; work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgbrfs_work", dynlib: libName.}
proc cgbrfsx*(matrixLayout: cint; trans: char; equed: char; n: LapackInt; kl: LapackInt;
             ku: LapackInt; nrhs: LapackInt; ab: ptr LapackComplexFloat;
             ldab: LapackInt; afb: ptr LapackComplexFloat; ldafb: LapackInt;
             ipiv: ptr LapackInt; r: ptr cfloat; c: ptr cfloat;
             b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
             ldx: LapackInt; rcond: ptr cfloat; berr: ptr cfloat; nErrBnds: LapackInt;
             errBndsNorm: ptr cfloat; errBndsComp: ptr cfloat; nparams: LapackInt;
             params: ptr cfloat): LapackInt {.importc: "LAPACKE_cgbrfsx",
    dynlib: libName.}
proc cgbrfsxWork*(matrixLayout: cint; trans: char; equed: char; n: LapackInt;
                 kl: LapackInt; ku: LapackInt; nrhs: LapackInt;
                 ab: ptr LapackComplexFloat; ldab: LapackInt;
                 afb: ptr LapackComplexFloat; ldafb: LapackInt; ipiv: ptr LapackInt;
                 r: ptr cfloat; c: ptr cfloat; b: ptr LapackComplexFloat; ldb: LapackInt;
                 x: ptr LapackComplexFloat; ldx: LapackInt; rcond: ptr cfloat;
                 berr: ptr cfloat; nErrBnds: LapackInt; errBndsNorm: ptr cfloat;
                 errBndsComp: ptr cfloat; nparams: LapackInt; params: ptr cfloat;
                 work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgbrfsx_work", dynlib: libName.}
proc cgbsv*(matrixLayout: cint; n: LapackInt; kl: LapackInt; ku: LapackInt;
           nrhs: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
           ipiv: ptr LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cgbsv", dynlib: libName.}
proc cgbsvWork*(matrixLayout: cint; n: LapackInt; kl: LapackInt; ku: LapackInt;
               nrhs: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
               ipiv: ptr LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cgbsv_work", dynlib: libName.}
proc cgbsvx*(matrixLayout: cint; fact: char; trans: char; n: LapackInt; kl: LapackInt;
            ku: LapackInt; nrhs: LapackInt; ab: ptr LapackComplexFloat;
            ldab: LapackInt; afb: ptr LapackComplexFloat; ldafb: LapackInt;
            ipiv: ptr LapackInt; equed: cstring; r: ptr cfloat; c: ptr cfloat;
            b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
            ldx: LapackInt; rcond: ptr cfloat; ferr: ptr cfloat; berr: ptr cfloat;
            rpivot: ptr cfloat): LapackInt {.importc: "LAPACKE_cgbsvx",
    dynlib: libName.}
proc cgbsvxWork*(matrixLayout: cint; fact: char; trans: char; n: LapackInt;
                kl: LapackInt; ku: LapackInt; nrhs: LapackInt;
                ab: ptr LapackComplexFloat; ldab: LapackInt;
                afb: ptr LapackComplexFloat; ldafb: LapackInt; ipiv: ptr LapackInt;
                equed: cstring; r: ptr cfloat; c: ptr cfloat; b: ptr LapackComplexFloat;
                ldb: LapackInt; x: ptr LapackComplexFloat; ldx: LapackInt;
                rcond: ptr cfloat; ferr: ptr cfloat; berr: ptr cfloat;
                work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgbsvx_work", dynlib: libName.}
proc cgbsvxx*(matrixLayout: cint; fact: char; trans: char; n: LapackInt; kl: LapackInt;
             ku: LapackInt; nrhs: LapackInt; ab: ptr LapackComplexFloat;
             ldab: LapackInt; afb: ptr LapackComplexFloat; ldafb: LapackInt;
             ipiv: ptr LapackInt; equed: cstring; r: ptr cfloat; c: ptr cfloat;
             b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
             ldx: LapackInt; rcond: ptr cfloat; rpvgrw: ptr cfloat; berr: ptr cfloat;
             nErrBnds: LapackInt; errBndsNorm: ptr cfloat; errBndsComp: ptr cfloat;
             nparams: LapackInt; params: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgbsvxx", dynlib: libName.}
proc cgbsvxxWork*(matrixLayout: cint; fact: char; trans: char; n: LapackInt;
                 kl: LapackInt; ku: LapackInt; nrhs: LapackInt;
                 ab: ptr LapackComplexFloat; ldab: LapackInt;
                 afb: ptr LapackComplexFloat; ldafb: LapackInt; ipiv: ptr LapackInt;
                 equed: cstring; r: ptr cfloat; c: ptr cfloat;
                 b: ptr LapackComplexFloat; ldb: LapackInt;
                 x: ptr LapackComplexFloat; ldx: LapackInt; rcond: ptr cfloat;
                 rpvgrw: ptr cfloat; berr: ptr cfloat; nErrBnds: LapackInt;
                 errBndsNorm: ptr cfloat; errBndsComp: ptr cfloat; nparams: LapackInt;
                 params: ptr cfloat; work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgbsvxx_work", dynlib: libName.}
proc cgbtrf*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
            ku: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
            ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_cgbtrf",
    dynlib: libName.}
proc cgbtrfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
                ku: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
                ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_cgbtrf_work",
    dynlib: libName.}
proc cgbtrs*(matrixLayout: cint; trans: char; n: LapackInt; kl: LapackInt; ku: LapackInt;
            nrhs: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
            ipiv: ptr LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cgbtrs", dynlib: libName.}
proc cgbtrsWork*(matrixLayout: cint; trans: char; n: LapackInt; kl: LapackInt;
                ku: LapackInt; nrhs: LapackInt; ab: ptr LapackComplexFloat;
                ldab: LapackInt; ipiv: ptr LapackInt; b: ptr LapackComplexFloat;
                ldb: LapackInt): LapackInt {.importc: "LAPACKE_cgbtrs_work",
    dynlib: libName.}
proc cgebak*(matrixLayout: cint; job: char; side: char; n: LapackInt; ilo: LapackInt;
            ihi: LapackInt; scale: ptr cfloat; m: LapackInt; v: ptr LapackComplexFloat;
            ldv: LapackInt): LapackInt {.importc: "LAPACKE_cgebak", dynlib: libName.}
proc cgebakWork*(matrixLayout: cint; job: char; side: char; n: LapackInt; ilo: LapackInt;
                ihi: LapackInt; scale: ptr cfloat; m: LapackInt;
                v: ptr LapackComplexFloat; ldv: LapackInt): LapackInt {.
    importc: "LAPACKE_cgebak_work", dynlib: libName.}
proc cgebal*(matrixLayout: cint; job: char; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; ilo: ptr LapackInt; ihi: ptr LapackInt; scale: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgebal", dynlib: libName.}
proc cgebalWork*(matrixLayout: cint; job: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; ilo: ptr LapackInt;
                ihi: ptr LapackInt; scale: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgebal_work", dynlib: libName.}
proc cgebrd*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; d: ptr cfloat; e: ptr cfloat; tauq: ptr LapackComplexFloat;
            taup: ptr LapackComplexFloat): LapackInt {.importc: "LAPACKE_cgebrd",
    dynlib: libName.}
proc cgebrdWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; d: ptr cfloat; e: ptr cfloat;
                tauq: ptr LapackComplexFloat; taup: ptr LapackComplexFloat;
                work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cgebrd_work", dynlib: libName.}
proc cgecon*(matrixLayout: cint; norm: char; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; anorm: cfloat; rcond: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgecon", dynlib: libName.}
proc cgeconWork*(matrixLayout: cint; norm: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; anorm: cfloat;
                rcond: ptr cfloat; work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgecon_work", dynlib: libName.}
proc cgeequ*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; r: ptr cfloat; c: ptr cfloat; rowcnd: ptr cfloat;
            colcnd: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgeequ", dynlib: libName.}
proc cgeequWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; r: ptr cfloat; c: ptr cfloat;
                rowcnd: ptr cfloat; colcnd: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgeequ_work", dynlib: libName.}
proc cgeequb*(matrixLayout: cint; m: LapackInt; n: LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; r: ptr cfloat; c: ptr cfloat;
             rowcnd: ptr cfloat; colcnd: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgeequb", dynlib: libName.}
proc cgeequbWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt; r: ptr cfloat;
                 c: ptr cfloat; rowcnd: ptr cfloat; colcnd: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgeequb_work", dynlib: libName.}
proc cgees*(matrixLayout: cint; jobvs: char; sort: char; select: Lapack_C_Select1;
           n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
           sdim: ptr LapackInt; w: ptr LapackComplexFloat; vs: ptr LapackComplexFloat;
           ldvs: LapackInt): LapackInt {.importc: "LAPACKE_cgees", dynlib: libName.}
proc cgeesWork*(matrixLayout: cint; jobvs: char; sort: char; select: Lapack_C_Select1;
               n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
               sdim: ptr LapackInt; w: ptr LapackComplexFloat;
               vs: ptr LapackComplexFloat; ldvs: LapackInt;
               work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat;
               bwork: ptr LapackLogical): LapackInt {.importc: "LAPACKE_cgees_work",
    dynlib: libName.}
proc cgeesx*(matrixLayout: cint; jobvs: char; sort: char; select: Lapack_C_Select1;
            sense: char; n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            sdim: ptr LapackInt; w: ptr LapackComplexFloat;
            vs: ptr LapackComplexFloat; ldvs: LapackInt; rconde: ptr cfloat;
            rcondv: ptr cfloat): LapackInt {.importc: "LAPACKE_cgeesx",
    dynlib: libName.}
proc cgeesxWork*(matrixLayout: cint; jobvs: char; sort: char; select: Lapack_C_Select1;
                sense: char; n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                sdim: ptr LapackInt; w: ptr LapackComplexFloat;
                vs: ptr LapackComplexFloat; ldvs: LapackInt; rconde: ptr cfloat;
                rcondv: ptr cfloat; work: ptr LapackComplexFloat; lwork: LapackInt;
                rwork: ptr cfloat; bwork: ptr LapackLogical): LapackInt {.
    importc: "LAPACKE_cgeesx_work", dynlib: libName.}
proc cgeev*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt;
           a: ptr LapackComplexFloat; lda: LapackInt; w: ptr LapackComplexFloat;
           vl: ptr LapackComplexFloat; ldvl: LapackInt; vr: ptr LapackComplexFloat;
           ldvr: LapackInt): LapackInt {.importc: "LAPACKE_cgeev", dynlib: libName.}
proc cgeevWork*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt;
               a: ptr LapackComplexFloat; lda: LapackInt; w: ptr LapackComplexFloat;
               vl: ptr LapackComplexFloat; ldvl: LapackInt;
               vr: ptr LapackComplexFloat; ldvr: LapackInt;
               work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgeev_work", dynlib: libName.}
proc cgeevx*(matrixLayout: cint; balanc: char; jobvl: char; jobvr: char; sense: char;
            n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            w: ptr LapackComplexFloat; vl: ptr LapackComplexFloat; ldvl: LapackInt;
            vr: ptr LapackComplexFloat; ldvr: LapackInt; ilo: ptr LapackInt;
            ihi: ptr LapackInt; scale: ptr cfloat; abnrm: ptr cfloat; rconde: ptr cfloat;
            rcondv: ptr cfloat): LapackInt {.importc: "LAPACKE_cgeevx",
    dynlib: libName.}
proc cgeevxWork*(matrixLayout: cint; balanc: char; jobvl: char; jobvr: char; sense: char;
                n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                w: ptr LapackComplexFloat; vl: ptr LapackComplexFloat;
                ldvl: LapackInt; vr: ptr LapackComplexFloat; ldvr: LapackInt;
                ilo: ptr LapackInt; ihi: ptr LapackInt; scale: ptr cfloat;
                abnrm: ptr cfloat; rconde: ptr cfloat; rcondv: ptr cfloat;
                work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgeevx_work", dynlib: libName.}
proc cgehrd*(matrixLayout: cint; n: LapackInt; ilo: LapackInt; ihi: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; tau: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cgehrd", dynlib: libName.}
proc cgehrdWork*(matrixLayout: cint; n: LapackInt; ilo: LapackInt; ihi: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt;
                tau: ptr LapackComplexFloat; work: ptr LapackComplexFloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_cgehrd_work",
    dynlib: libName.}
proc cgejsv*(matrixLayout: cint; joba: char; jobu: char; jobv: char; jobr: char;
            jobt: char; jobp: char; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; sva: ptr cfloat;
            u: ptr LapackComplexFloat; ldu: LapackInt; v: ptr LapackComplexFloat;
            ldv: LapackInt; stat: ptr cfloat; istat: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_cgejsv", dynlib: libName.}
proc cgejsvWork*(matrixLayout: cint; joba: char; jobu: char; jobv: char; jobr: char;
                jobt: char; jobp: char; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; sva: ptr cfloat;
                u: ptr LapackComplexFloat; ldu: LapackInt; v: ptr LapackComplexFloat;
                ldv: LapackInt; cwork: ptr LapackComplexFloat; lwork: LapackInt;
                rwork: ptr cfloat; lrwork: LapackInt; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_cgejsv_work", dynlib: libName.}
proc cgelq2*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; tau: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cgelq2", dynlib: libName.}
proc cgelq2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt;
                tau: ptr LapackComplexFloat; work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cgelq2_work", dynlib: libName.}
proc cgelqf*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; tau: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cgelqf", dynlib: libName.}
proc cgelqfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt;
                tau: ptr LapackComplexFloat; work: ptr LapackComplexFloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_cgelqf_work",
    dynlib: libName.}
proc cgels*(matrixLayout: cint; trans: char; m: LapackInt; n: LapackInt; nrhs: LapackInt;
           a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
           ldb: LapackInt): LapackInt {.importc: "LAPACKE_cgels", dynlib: libName.}
proc cgelsWork*(matrixLayout: cint; trans: char; m: LapackInt; n: LapackInt;
               nrhs: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
               b: ptr LapackComplexFloat; ldb: LapackInt;
               work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cgels_work", dynlib: libName.}
proc cgelsd*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
            ldb: LapackInt; s: ptr cfloat; rcond: cfloat; rank: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_cgelsd", dynlib: libName.}
proc cgelsdWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
                ldb: LapackInt; s: ptr cfloat; rcond: cfloat; rank: ptr LapackInt;
                work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_cgelsd_work",
    dynlib: libName.}
proc cgelss*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
            ldb: LapackInt; s: ptr cfloat; rcond: cfloat; rank: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_cgelss", dynlib: libName.}
proc cgelssWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
                ldb: LapackInt; s: ptr cfloat; rcond: cfloat; rank: ptr LapackInt;
                work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgelss_work", dynlib: libName.}
proc cgelsy*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
            ldb: LapackInt; jpvt: ptr LapackInt; rcond: cfloat; rank: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_cgelsy", dynlib: libName.}
proc cgelsyWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
                ldb: LapackInt; jpvt: ptr LapackInt; rcond: cfloat;
                rank: ptr LapackInt; work: ptr LapackComplexFloat; lwork: LapackInt;
                rwork: ptr cfloat): LapackInt {.importc: "LAPACKE_cgelsy_work",
    dynlib: libName.}
proc cgemqrt*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
             k: LapackInt; nb: LapackInt; v: ptr LapackComplexFloat; ldv: LapackInt;
             t: ptr LapackComplexFloat; ldt: LapackInt; c: ptr LapackComplexFloat;
             ldc: LapackInt): LapackInt {.importc: "LAPACKE_cgemqrt", dynlib: libName.}
proc cgemqrtWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt;
                 n: LapackInt; k: LapackInt; nb: LapackInt; v: ptr LapackComplexFloat;
                 ldv: LapackInt; t: ptr LapackComplexFloat; ldt: LapackInt;
                 c: ptr LapackComplexFloat; ldc: LapackInt;
                 work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cgemqrt_work", dynlib: libName.}
proc cgeqlf*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; tau: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cgeqlf", dynlib: libName.}
proc cgeqlfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt;
                tau: ptr LapackComplexFloat; work: ptr LapackComplexFloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_cgeqlf_work",
    dynlib: libName.}
proc cgeqp3*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; jpvt: ptr LapackInt; tau: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cgeqp3", dynlib: libName.}
proc cgeqp3Work*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; jpvt: ptr LapackInt;
                tau: ptr LapackComplexFloat; work: ptr LapackComplexFloat;
                lwork: LapackInt; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgeqp3_work", dynlib: libName.}
proc cgeqpf*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; jpvt: ptr LapackInt; tau: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cgeqpf", dynlib: libName.}
proc cgeqpfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; jpvt: ptr LapackInt;
                tau: ptr LapackComplexFloat; work: ptr LapackComplexFloat;
                rwork: ptr cfloat): LapackInt {.importc: "LAPACKE_cgeqpf_work",
    dynlib: libName.}
proc cgeqr2*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; tau: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cgeqr2", dynlib: libName.}
proc cgeqr2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt;
                tau: ptr LapackComplexFloat; work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cgeqr2_work", dynlib: libName.}
proc cgeqrf*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; tau: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cgeqrf", dynlib: libName.}
proc cgeqrfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt;
                tau: ptr LapackComplexFloat; work: ptr LapackComplexFloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_cgeqrf_work",
    dynlib: libName.}
proc cgeqrfp*(matrixLayout: cint; m: LapackInt; n: LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; tau: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cgeqrfp", dynlib: libName.}
proc cgeqrfpWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt;
                 tau: ptr LapackComplexFloat; work: ptr LapackComplexFloat;
                 lwork: LapackInt): LapackInt {.importc: "LAPACKE_cgeqrfp_work",
    dynlib: libName.}
proc cgeqrt*(matrixLayout: cint; m: LapackInt; n: LapackInt; nb: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; t: ptr LapackComplexFloat;
            ldt: LapackInt): LapackInt {.importc: "LAPACKE_cgeqrt", dynlib: libName.}
proc cgeqrt2*(matrixLayout: cint; m: LapackInt; n: LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; t: ptr LapackComplexFloat;
             ldt: LapackInt): LapackInt {.importc: "LAPACKE_cgeqrt2", dynlib: libName.}
proc cgeqrt2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt;
                 t: ptr LapackComplexFloat; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_cgeqrt2_work", dynlib: libName.}
proc cgeqrt3*(matrixLayout: cint; m: LapackInt; n: LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; t: ptr LapackComplexFloat;
             ldt: LapackInt): LapackInt {.importc: "LAPACKE_cgeqrt3", dynlib: libName.}
proc cgeqrt3Work*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt;
                 t: ptr LapackComplexFloat; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_cgeqrt3_work", dynlib: libName.}
proc cgeqrtWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; nb: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; t: ptr LapackComplexFloat;
                ldt: LapackInt; work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cgeqrt_work", dynlib: libName.}
proc cgerfs*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; af: ptr LapackComplexFloat;
            ldaf: LapackInt; ipiv: ptr LapackInt; b: ptr LapackComplexFloat;
            ldb: LapackInt; x: ptr LapackComplexFloat; ldx: LapackInt;
            ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.importc: "LAPACKE_cgerfs",
    dynlib: libName.}
proc cgerfsWork*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt;
                af: ptr LapackComplexFloat; ldaf: LapackInt; ipiv: ptr LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
                ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat;
                work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgerfs_work", dynlib: libName.}
proc cgerfsx*(matrixLayout: cint; trans: char; equed: char; n: LapackInt;
             nrhs: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
             af: ptr LapackComplexFloat; ldaf: LapackInt; ipiv: ptr LapackInt;
             r: ptr cfloat; c: ptr cfloat; b: ptr LapackComplexFloat; ldb: LapackInt;
             x: ptr LapackComplexFloat; ldx: LapackInt; rcond: ptr cfloat;
             berr: ptr cfloat; nErrBnds: LapackInt; errBndsNorm: ptr cfloat;
             errBndsComp: ptr cfloat; nparams: LapackInt; params: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgerfsx", dynlib: libName.}
proc cgerfsxWork*(matrixLayout: cint; trans: char; equed: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                 af: ptr LapackComplexFloat; ldaf: LapackInt; ipiv: ptr LapackInt;
                 r: ptr cfloat; c: ptr cfloat; b: ptr LapackComplexFloat; ldb: LapackInt;
                 x: ptr LapackComplexFloat; ldx: LapackInt; rcond: ptr cfloat;
                 berr: ptr cfloat; nErrBnds: LapackInt; errBndsNorm: ptr cfloat;
                 errBndsComp: ptr cfloat; nparams: LapackInt; params: ptr cfloat;
                 work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgerfsx_work", dynlib: libName.}
proc cgerqf*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; tau: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cgerqf", dynlib: libName.}
proc cgerqfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt;
                tau: ptr LapackComplexFloat; work: ptr LapackComplexFloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_cgerqf_work",
    dynlib: libName.}
proc cgesdd*(matrixLayout: cint; jobz: char; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; s: ptr cfloat;
            u: ptr LapackComplexFloat; ldu: LapackInt; vt: ptr LapackComplexFloat;
            ldvt: LapackInt): LapackInt {.importc: "LAPACKE_cgesdd", dynlib: libName.}
proc cgesddWork*(matrixLayout: cint; jobz: char; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; s: ptr cfloat;
                u: ptr LapackComplexFloat; ldu: LapackInt;
                vt: ptr LapackComplexFloat; ldvt: LapackInt;
                work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_cgesdd_work",
    dynlib: libName.}
proc cgesv*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt;
           a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
           b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cgesv", dynlib: libName.}
proc cgesvWork*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt;
               a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
               b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cgesv_work", dynlib: libName.}
proc cgesvd*(matrixLayout: cint; jobu: char; jobvt: char; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; s: ptr cfloat;
            u: ptr LapackComplexFloat; ldu: LapackInt; vt: ptr LapackComplexFloat;
            ldvt: LapackInt; superb: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgesvd", dynlib: libName.}
proc cgesvdWork*(matrixLayout: cint; jobu: char; jobvt: char; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; s: ptr cfloat;
                u: ptr LapackComplexFloat; ldu: LapackInt;
                vt: ptr LapackComplexFloat; ldvt: LapackInt;
                work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgesvd_work", dynlib: libName.}
proc cgesvdx*(matrixLayout: cint; jobu: char; jobvt: char; range: char; m: LapackInt;
             n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt; vl: cfloat;
             vu: cfloat; il: LapackInt; iu: LapackInt; ns: ptr LapackInt; s: ptr cfloat;
             u: ptr LapackComplexFloat; ldu: LapackInt; vt: ptr LapackComplexFloat;
             ldvt: LapackInt; superb: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_cgesvdx", dynlib: libName.}
proc cgesvdxWork*(matrixLayout: cint; jobu: char; jobvt: char; range: char; m: LapackInt;
                 n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt; vl: cfloat;
                 vu: cfloat; il: LapackInt; iu: LapackInt; ns: ptr LapackInt;
                 s: ptr cfloat; u: ptr LapackComplexFloat; ldu: LapackInt;
                 vt: ptr LapackComplexFloat; ldvt: LapackInt;
                 work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat;
                 iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_cgesvdx_work",
    dynlib: libName.}
proc cgesvj*(matrixLayout: cint; joba: char; jobu: char; jobv: char; m: LapackInt;
            n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt; sva: ptr cfloat;
            mv: LapackInt; v: ptr LapackComplexFloat; ldv: LapackInt; stat: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgesvj", dynlib: libName.}
proc cgesvjWork*(matrixLayout: cint; joba: char; jobu: char; jobv: char; m: LapackInt;
                n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                sva: ptr cfloat; mv: LapackInt; v: ptr LapackComplexFloat;
                ldv: LapackInt; cwork: ptr LapackComplexFloat; lwork: LapackInt;
                rwork: ptr cfloat; lrwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cgesvj_work", dynlib: libName.}
proc cgesvx*(matrixLayout: cint; fact: char; trans: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; af: ptr LapackComplexFloat;
            ldaf: LapackInt; ipiv: ptr LapackInt; equed: cstring; r: ptr cfloat;
            c: ptr cfloat; b: ptr LapackComplexFloat; ldb: LapackInt;
            x: ptr LapackComplexFloat; ldx: LapackInt; rcond: ptr cfloat;
            ferr: ptr cfloat; berr: ptr cfloat; rpivot: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgesvx", dynlib: libName.}
proc cgesvxWork*(matrixLayout: cint; fact: char; trans: char; n: LapackInt;
                nrhs: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                af: ptr LapackComplexFloat; ldaf: LapackInt; ipiv: ptr LapackInt;
                equed: cstring; r: ptr cfloat; c: ptr cfloat; b: ptr LapackComplexFloat;
                ldb: LapackInt; x: ptr LapackComplexFloat; ldx: LapackInt;
                rcond: ptr cfloat; ferr: ptr cfloat; berr: ptr cfloat;
                work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgesvx_work", dynlib: libName.}
proc cgesvxx*(matrixLayout: cint; fact: char; trans: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; af: ptr LapackComplexFloat;
             ldaf: LapackInt; ipiv: ptr LapackInt; equed: cstring; r: ptr cfloat;
             c: ptr cfloat; b: ptr LapackComplexFloat; ldb: LapackInt;
             x: ptr LapackComplexFloat; ldx: LapackInt; rcond: ptr cfloat;
             rpvgrw: ptr cfloat; berr: ptr cfloat; nErrBnds: LapackInt;
             errBndsNorm: ptr cfloat; errBndsComp: ptr cfloat; nparams: LapackInt;
             params: ptr cfloat): LapackInt {.importc: "LAPACKE_cgesvxx",
    dynlib: libName.}
proc cgesvxxWork*(matrixLayout: cint; fact: char; trans: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                 af: ptr LapackComplexFloat; ldaf: LapackInt; ipiv: ptr LapackInt;
                 equed: cstring; r: ptr cfloat; c: ptr cfloat;
                 b: ptr LapackComplexFloat; ldb: LapackInt;
                 x: ptr LapackComplexFloat; ldx: LapackInt; rcond: ptr cfloat;
                 rpvgrw: ptr cfloat; berr: ptr cfloat; nErrBnds: LapackInt;
                 errBndsNorm: ptr cfloat; errBndsComp: ptr cfloat; nparams: LapackInt;
                 params: ptr cfloat; work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgesvxx_work", dynlib: libName.}
proc cgetf2*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_cgetf2", dynlib: libName.}
proc cgetf2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_cgetf2_work", dynlib: libName.}
proc cgetrf*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_cgetrf", dynlib: libName.}
proc cgetrf2*(matrixLayout: cint; m: LapackInt; n: LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_cgetrf2", dynlib: libName.}
proc cgetrf2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_cgetrf2_work", dynlib: libName.}
proc cgetrfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_cgetrf_work", dynlib: libName.}
proc cgetri*(matrixLayout: cint; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_cgetri", dynlib: libName.}
proc cgetriWork*(matrixLayout: cint; n: LapackInt; a: ptr LapackComplexFloat;
                lda: LapackInt; ipiv: ptr LapackInt; work: ptr LapackComplexFloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_cgetri_work",
    dynlib: libName.}
proc cgetrs*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
            b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cgetrs", dynlib: libName.}
proc cgetrsWork*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cgetrs_work", dynlib: libName.}
proc cggbak*(matrixLayout: cint; job: char; side: char; n: LapackInt; ilo: LapackInt;
            ihi: LapackInt; lscale: ptr cfloat; rscale: ptr cfloat; m: LapackInt;
            v: ptr LapackComplexFloat; ldv: LapackInt): LapackInt {.
    importc: "LAPACKE_cggbak", dynlib: libName.}
proc cggbakWork*(matrixLayout: cint; job: char; side: char; n: LapackInt; ilo: LapackInt;
                ihi: LapackInt; lscale: ptr cfloat; rscale: ptr cfloat; m: LapackInt;
                v: ptr LapackComplexFloat; ldv: LapackInt): LapackInt {.
    importc: "LAPACKE_cggbak_work", dynlib: libName.}
proc cggbal*(matrixLayout: cint; job: char; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt;
            ilo: ptr LapackInt; ihi: ptr LapackInt; lscale: ptr cfloat;
            rscale: ptr cfloat): LapackInt {.importc: "LAPACKE_cggbal",
    dynlib: libName.}
proc cggbalWork*(matrixLayout: cint; job: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
                ldb: LapackInt; ilo: ptr LapackInt; ihi: ptr LapackInt;
                lscale: ptr cfloat; rscale: ptr cfloat; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cggbal_work", dynlib: libName.}
proc cgges*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
           selctg: Lapack_C_Select2; n: LapackInt; a: ptr LapackComplexFloat;
           lda: LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt;
           sdim: ptr LapackInt; alpha: ptr LapackComplexFloat;
           beta: ptr LapackComplexFloat; vsl: ptr LapackComplexFloat;
           ldvsl: LapackInt; vsr: ptr LapackComplexFloat; ldvsr: LapackInt): LapackInt {.
    importc: "LAPACKE_cgges", dynlib: libName.}
proc cgges3*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
            selctg: Lapack_C_Select2; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt;
            sdim: ptr LapackInt; alpha: ptr LapackComplexFloat;
            beta: ptr LapackComplexFloat; vsl: ptr LapackComplexFloat;
            ldvsl: LapackInt; vsr: ptr LapackComplexFloat; ldvsr: LapackInt): LapackInt {.
    importc: "LAPACKE_cgges3", dynlib: libName.}
proc cgges3Work*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
                selctg: Lapack_C_Select2; n: LapackInt; a: ptr LapackComplexFloat;
                lda: LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt;
                sdim: ptr LapackInt; alpha: ptr LapackComplexFloat;
                beta: ptr LapackComplexFloat; vsl: ptr LapackComplexFloat;
                ldvsl: LapackInt; vsr: ptr LapackComplexFloat; ldvsr: LapackInt;
                work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat;
                bwork: ptr LapackLogical): LapackInt {.
    importc: "LAPACKE_cgges3_work", dynlib: libName.}
proc cggesWork*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
               selctg: Lapack_C_Select2; n: LapackInt; a: ptr LapackComplexFloat;
               lda: LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt;
               sdim: ptr LapackInt; alpha: ptr LapackComplexFloat;
               beta: ptr LapackComplexFloat; vsl: ptr LapackComplexFloat;
               ldvsl: LapackInt; vsr: ptr LapackComplexFloat; ldvsr: LapackInt;
               work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat;
               bwork: ptr LapackLogical): LapackInt {.importc: "LAPACKE_cgges_work",
    dynlib: libName.}
proc cggesx*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
            selctg: Lapack_C_Select2; sense: char; n: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
            ldb: LapackInt; sdim: ptr LapackInt; alpha: ptr LapackComplexFloat;
            beta: ptr LapackComplexFloat; vsl: ptr LapackComplexFloat;
            ldvsl: LapackInt; vsr: ptr LapackComplexFloat; ldvsr: LapackInt;
            rconde: ptr cfloat; rcondv: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cggesx", dynlib: libName.}
proc cggesxWork*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
                selctg: Lapack_C_Select2; sense: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
                ldb: LapackInt; sdim: ptr LapackInt; alpha: ptr LapackComplexFloat;
                beta: ptr LapackComplexFloat; vsl: ptr LapackComplexFloat;
                ldvsl: LapackInt; vsr: ptr LapackComplexFloat; ldvsr: LapackInt;
                rconde: ptr cfloat; rcondv: ptr cfloat; work: ptr LapackComplexFloat;
                lwork: LapackInt; rwork: ptr cfloat; iwork: ptr LapackInt;
                liwork: LapackInt; bwork: ptr LapackLogical): LapackInt {.
    importc: "LAPACKE_cggesx_work", dynlib: libName.}
proc cggev*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt;
           a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
           ldb: LapackInt; alpha: ptr LapackComplexFloat;
           beta: ptr LapackComplexFloat; vl: ptr LapackComplexFloat; ldvl: LapackInt;
           vr: ptr LapackComplexFloat; ldvr: LapackInt): LapackInt {.
    importc: "LAPACKE_cggev", dynlib: libName.}
proc cggev3*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
            ldb: LapackInt; alpha: ptr LapackComplexFloat;
            beta: ptr LapackComplexFloat; vl: ptr LapackComplexFloat; ldvl: LapackInt;
            vr: ptr LapackComplexFloat; ldvr: LapackInt): LapackInt {.
    importc: "LAPACKE_cggev3", dynlib: libName.}
proc cggev3Work*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
                ldb: LapackInt; alpha: ptr LapackComplexFloat;
                beta: ptr LapackComplexFloat; vl: ptr LapackComplexFloat;
                ldvl: LapackInt; vr: ptr LapackComplexFloat; ldvr: LapackInt;
                work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cggev3_work", dynlib: libName.}
proc cggevWork*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt;
               a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
               ldb: LapackInt; alpha: ptr LapackComplexFloat;
               beta: ptr LapackComplexFloat; vl: ptr LapackComplexFloat;
               ldvl: LapackInt; vr: ptr LapackComplexFloat; ldvr: LapackInt;
               work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cggev_work", dynlib: libName.}
proc cggevx*(matrixLayout: cint; balanc: char; jobvl: char; jobvr: char; sense: char;
            n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            b: ptr LapackComplexFloat; ldb: LapackInt; alpha: ptr LapackComplexFloat;
            beta: ptr LapackComplexFloat; vl: ptr LapackComplexFloat; ldvl: LapackInt;
            vr: ptr LapackComplexFloat; ldvr: LapackInt; ilo: ptr LapackInt;
            ihi: ptr LapackInt; lscale: ptr cfloat; rscale: ptr cfloat;
            abnrm: ptr cfloat; bbnrm: ptr cfloat; rconde: ptr cfloat; rcondv: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cggevx", dynlib: libName.}
proc cggevxWork*(matrixLayout: cint; balanc: char; jobvl: char; jobvr: char; sense: char;
                n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt;
                alpha: ptr LapackComplexFloat; beta: ptr LapackComplexFloat;
                vl: ptr LapackComplexFloat; ldvl: LapackInt;
                vr: ptr LapackComplexFloat; ldvr: LapackInt; ilo: ptr LapackInt;
                ihi: ptr LapackInt; lscale: ptr cfloat; rscale: ptr cfloat;
                abnrm: ptr cfloat; bbnrm: ptr cfloat; rconde: ptr cfloat;
                rcondv: ptr cfloat; work: ptr LapackComplexFloat; lwork: LapackInt;
                rwork: ptr cfloat; iwork: ptr LapackInt; bwork: ptr LapackLogical): LapackInt {.
    importc: "LAPACKE_cggevx_work", dynlib: libName.}
proc cggglm*(matrixLayout: cint; n: LapackInt; m: LapackInt; p: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
            ldb: LapackInt; d: ptr LapackComplexFloat; x: ptr LapackComplexFloat;
            y: ptr LapackComplexFloat): LapackInt {.importc: "LAPACKE_cggglm",
    dynlib: libName.}
proc cggglmWork*(matrixLayout: cint; n: LapackInt; m: LapackInt; p: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
                ldb: LapackInt; d: ptr LapackComplexFloat; x: ptr LapackComplexFloat;
                y: ptr LapackComplexFloat; work: ptr LapackComplexFloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_cggglm_work",
    dynlib: libName.}
proc cgghd3*(matrixLayout: cint; compq: char; compz: char; n: LapackInt; ilo: LapackInt;
            ihi: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            b: ptr LapackComplexFloat; ldb: LapackInt; q: ptr LapackComplexFloat;
            ldq: LapackInt; z: ptr LapackComplexFloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_cgghd3", dynlib: libName.}
proc cgghd3Work*(matrixLayout: cint; compq: char; compz: char; n: LapackInt;
                ilo: LapackInt; ihi: LapackInt; a: ptr LapackComplexFloat;
                lda: LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt;
                q: ptr LapackComplexFloat; ldq: LapackInt; z: ptr LapackComplexFloat;
                ldz: LapackInt; work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cgghd3_work", dynlib: libName.}
proc cgghrd*(matrixLayout: cint; compq: char; compz: char; n: LapackInt; ilo: LapackInt;
            ihi: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            b: ptr LapackComplexFloat; ldb: LapackInt; q: ptr LapackComplexFloat;
            ldq: LapackInt; z: ptr LapackComplexFloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_cgghrd", dynlib: libName.}
proc cgghrdWork*(matrixLayout: cint; compq: char; compz: char; n: LapackInt;
                ilo: LapackInt; ihi: LapackInt; a: ptr LapackComplexFloat;
                lda: LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt;
                q: ptr LapackComplexFloat; ldq: LapackInt; z: ptr LapackComplexFloat;
                ldz: LapackInt): LapackInt {.importc: "LAPACKE_cgghrd_work",
    dynlib: libName.}
proc cgglse*(matrixLayout: cint; m: LapackInt; n: LapackInt; p: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
            ldb: LapackInt; c: ptr LapackComplexFloat; d: ptr LapackComplexFloat;
            x: ptr LapackComplexFloat): LapackInt {.importc: "LAPACKE_cgglse",
    dynlib: libName.}
proc cgglseWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; p: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
                ldb: LapackInt; c: ptr LapackComplexFloat; d: ptr LapackComplexFloat;
                x: ptr LapackComplexFloat; work: ptr LapackComplexFloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_cgglse_work",
    dynlib: libName.}
proc cggqrf*(matrixLayout: cint; n: LapackInt; m: LapackInt; p: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; taua: ptr LapackComplexFloat;
            b: ptr LapackComplexFloat; ldb: LapackInt; taub: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cggqrf", dynlib: libName.}
proc cggqrfWork*(matrixLayout: cint; n: LapackInt; m: LapackInt; p: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt;
                taua: ptr LapackComplexFloat; b: ptr LapackComplexFloat;
                ldb: LapackInt; taub: ptr LapackComplexFloat;
                work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cggqrf_work", dynlib: libName.}
proc cggrqf*(matrixLayout: cint; m: LapackInt; p: LapackInt; n: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; taua: ptr LapackComplexFloat;
            b: ptr LapackComplexFloat; ldb: LapackInt; taub: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cggrqf", dynlib: libName.}
proc cggrqfWork*(matrixLayout: cint; m: LapackInt; p: LapackInt; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt;
                taua: ptr LapackComplexFloat; b: ptr LapackComplexFloat;
                ldb: LapackInt; taub: ptr LapackComplexFloat;
                work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cggrqf_work", dynlib: libName.}
proc cggsvd*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
            n: LapackInt; p: LapackInt; k: ptr LapackInt; l: ptr LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
            ldb: LapackInt; alpha: ptr cfloat; beta: ptr cfloat;
            u: ptr LapackComplexFloat; ldu: LapackInt; v: ptr LapackComplexFloat;
            ldv: LapackInt; q: ptr LapackComplexFloat; ldq: LapackInt;
            iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_cggsvd",
    dynlib: libName.}
proc cggsvd3*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
             n: LapackInt; p: LapackInt; k: ptr LapackInt; l: ptr LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
             ldb: LapackInt; alpha: ptr cfloat; beta: ptr cfloat;
             u: ptr LapackComplexFloat; ldu: LapackInt; v: ptr LapackComplexFloat;
             ldv: LapackInt; q: ptr LapackComplexFloat; ldq: LapackInt;
             iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_cggsvd3",
    dynlib: libName.}
proc cggsvd3Work*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
                 n: LapackInt; p: LapackInt; k: ptr LapackInt; l: ptr LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt;
                 b: ptr LapackComplexFloat; ldb: LapackInt; alpha: ptr cfloat;
                 beta: ptr cfloat; u: ptr LapackComplexFloat; ldu: LapackInt;
                 v: ptr LapackComplexFloat; ldv: LapackInt;
                 q: ptr LapackComplexFloat; ldq: LapackInt;
                 work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat;
                 iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_cggsvd3_work",
    dynlib: libName.}
proc cggsvdWork*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
                n: LapackInt; p: LapackInt; k: ptr LapackInt; l: ptr LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
                ldb: LapackInt; alpha: ptr cfloat; beta: ptr cfloat;
                u: ptr LapackComplexFloat; ldu: LapackInt; v: ptr LapackComplexFloat;
                ldv: LapackInt; q: ptr LapackComplexFloat; ldq: LapackInt;
                work: ptr LapackComplexFloat; rwork: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_cggsvd_work", dynlib: libName.}
proc cggsvp*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
            p: LapackInt; n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            b: ptr LapackComplexFloat; ldb: LapackInt; tola: cfloat; tolb: cfloat;
            k: ptr LapackInt; l: ptr LapackInt; u: ptr LapackComplexFloat;
            ldu: LapackInt; v: ptr LapackComplexFloat; ldv: LapackInt;
            q: ptr LapackComplexFloat; ldq: LapackInt): LapackInt {.
    importc: "LAPACKE_cggsvp", dynlib: libName.}
proc cggsvp3*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
             p: LapackInt; n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
             b: ptr LapackComplexFloat; ldb: LapackInt; tola: cfloat; tolb: cfloat;
             k: ptr LapackInt; l: ptr LapackInt; u: ptr LapackComplexFloat;
             ldu: LapackInt; v: ptr LapackComplexFloat; ldv: LapackInt;
             q: ptr LapackComplexFloat; ldq: LapackInt): LapackInt {.
    importc: "LAPACKE_cggsvp3", dynlib: libName.}
proc cggsvp3Work*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
                 p: LapackInt; n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                 b: ptr LapackComplexFloat; ldb: LapackInt; tola: cfloat; tolb: cfloat;
                 k: ptr LapackInt; l: ptr LapackInt; u: ptr LapackComplexFloat;
                 ldu: LapackInt; v: ptr LapackComplexFloat; ldv: LapackInt;
                 q: ptr LapackComplexFloat; ldq: LapackInt; iwork: ptr LapackInt;
                 rwork: ptr cfloat; tau: ptr LapackComplexFloat;
                 work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cggsvp3_work", dynlib: libName.}
proc cggsvpWork*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
                p: LapackInt; n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt; tola: cfloat; tolb: cfloat;
                k: ptr LapackInt; l: ptr LapackInt; u: ptr LapackComplexFloat;
                ldu: LapackInt; v: ptr LapackComplexFloat; ldv: LapackInt;
                q: ptr LapackComplexFloat; ldq: LapackInt; iwork: ptr LapackInt;
                rwork: ptr cfloat; tau: ptr LapackComplexFloat;
                work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cggsvp_work", dynlib: libName.}
proc cgtcon*(norm: char; n: LapackInt; dl: ptr LapackComplexFloat;
            d: ptr LapackComplexFloat; du: ptr LapackComplexFloat;
            du2: ptr LapackComplexFloat; ipiv: ptr LapackInt; anorm: cfloat;
            rcond: ptr cfloat): LapackInt {.importc: "LAPACKE_cgtcon", dynlib: libName.}
proc cgtconWork*(norm: char; n: LapackInt; dl: ptr LapackComplexFloat;
                d: ptr LapackComplexFloat; du: ptr LapackComplexFloat;
                du2: ptr LapackComplexFloat; ipiv: ptr LapackInt; anorm: cfloat;
                rcond: ptr cfloat; work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cgtcon_work", dynlib: libName.}
proc cgtrfs*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
            dl: ptr LapackComplexFloat; d: ptr LapackComplexFloat;
            du: ptr LapackComplexFloat; dlf: ptr LapackComplexFloat;
            df: ptr LapackComplexFloat; duf: ptr LapackComplexFloat;
            du2: ptr LapackComplexFloat; ipiv: ptr LapackInt;
            b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
            ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgtrfs", dynlib: libName.}
proc cgtrfsWork*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
                dl: ptr LapackComplexFloat; d: ptr LapackComplexFloat;
                du: ptr LapackComplexFloat; dlf: ptr LapackComplexFloat;
                df: ptr LapackComplexFloat; duf: ptr LapackComplexFloat;
                du2: ptr LapackComplexFloat; ipiv: ptr LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
                ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat;
                work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgtrfs_work", dynlib: libName.}
proc cgtsv*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt;
           dl: ptr LapackComplexFloat; d: ptr LapackComplexFloat;
           du: ptr LapackComplexFloat; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cgtsv", dynlib: libName.}
proc cgtsvWork*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt;
               dl: ptr LapackComplexFloat; d: ptr LapackComplexFloat;
               du: ptr LapackComplexFloat; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cgtsv_work", dynlib: libName.}
proc cgtsvx*(matrixLayout: cint; fact: char; trans: char; n: LapackInt; nrhs: LapackInt;
            dl: ptr LapackComplexFloat; d: ptr LapackComplexFloat;
            du: ptr LapackComplexFloat; dlf: ptr LapackComplexFloat;
            df: ptr LapackComplexFloat; duf: ptr LapackComplexFloat;
            du2: ptr LapackComplexFloat; ipiv: ptr LapackInt;
            b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
            ldx: LapackInt; rcond: ptr cfloat; ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cgtsvx", dynlib: libName.}
proc cgtsvxWork*(matrixLayout: cint; fact: char; trans: char; n: LapackInt;
                nrhs: LapackInt; dl: ptr LapackComplexFloat;
                d: ptr LapackComplexFloat; du: ptr LapackComplexFloat;
                dlf: ptr LapackComplexFloat; df: ptr LapackComplexFloat;
                duf: ptr LapackComplexFloat; du2: ptr LapackComplexFloat;
                ipiv: ptr LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt;
                x: ptr LapackComplexFloat; ldx: LapackInt; rcond: ptr cfloat;
                ferr: ptr cfloat; berr: ptr cfloat; work: ptr LapackComplexFloat;
                rwork: ptr cfloat): LapackInt {.importc: "LAPACKE_cgtsvx_work",
    dynlib: libName.}
proc cgttrf*(n: LapackInt; dl: ptr LapackComplexFloat; d: ptr LapackComplexFloat;
            du: ptr LapackComplexFloat; du2: ptr LapackComplexFloat;
            ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_cgttrf",
    dynlib: libName.}
proc cgttrfWork*(n: LapackInt; dl: ptr LapackComplexFloat; d: ptr LapackComplexFloat;
                du: ptr LapackComplexFloat; du2: ptr LapackComplexFloat;
                ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_cgttrf_work",
    dynlib: libName.}
proc cgttrs*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
            dl: ptr LapackComplexFloat; d: ptr LapackComplexFloat;
            du: ptr LapackComplexFloat; du2: ptr LapackComplexFloat;
            ipiv: ptr LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cgttrs", dynlib: libName.}
proc cgttrsWork*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
                dl: ptr LapackComplexFloat; d: ptr LapackComplexFloat;
                du: ptr LapackComplexFloat; du2: ptr LapackComplexFloat;
                ipiv: ptr LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cgttrs_work", dynlib: libName.}
proc chbev*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; kd: LapackInt;
           ab: ptr LapackComplexFloat; ldab: LapackInt; w: ptr cfloat;
           z: ptr LapackComplexFloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_chbev", dynlib: libName.}
proc chbevWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; kd: LapackInt;
               ab: ptr LapackComplexFloat; ldab: LapackInt; w: ptr cfloat;
               z: ptr LapackComplexFloat; ldz: LapackInt;
               work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_chbev_work", dynlib: libName.}
proc chbevd*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; kd: LapackInt;
            ab: ptr LapackComplexFloat; ldab: LapackInt; w: ptr cfloat;
            z: ptr LapackComplexFloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_chbevd", dynlib: libName.}
proc chbevdWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; kd: LapackInt;
                ab: ptr LapackComplexFloat; ldab: LapackInt; w: ptr cfloat;
                z: ptr LapackComplexFloat; ldz: LapackInt;
                work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat;
                lrwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_chbevd_work", dynlib: libName.}
proc chbevx*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
            kd: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
            q: ptr LapackComplexFloat; ldq: LapackInt; vl: cfloat; vu: cfloat;
            il: LapackInt; iu: LapackInt; abstol: cfloat; m: ptr LapackInt;
            w: ptr cfloat; z: ptr LapackComplexFloat; ldz: LapackInt;
            ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_chbevx",
    dynlib: libName.}
proc chbevxWork*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                kd: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
                q: ptr LapackComplexFloat; ldq: LapackInt; vl: cfloat; vu: cfloat;
                il: LapackInt; iu: LapackInt; abstol: cfloat; m: ptr LapackInt;
                w: ptr cfloat; z: ptr LapackComplexFloat; ldz: LapackInt;
                work: ptr LapackComplexFloat; rwork: ptr cfloat; iwork: ptr LapackInt;
                ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_chbevx_work",
    dynlib: libName.}
proc chbgst*(matrixLayout: cint; vect: char; uplo: char; n: LapackInt; ka: LapackInt;
            kb: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
            bb: ptr LapackComplexFloat; ldbb: LapackInt; x: ptr LapackComplexFloat;
            ldx: LapackInt): LapackInt {.importc: "LAPACKE_chbgst", dynlib: libName.}
proc chbgstWork*(matrixLayout: cint; vect: char; uplo: char; n: LapackInt; ka: LapackInt;
                kb: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
                bb: ptr LapackComplexFloat; ldbb: LapackInt;
                x: ptr LapackComplexFloat; ldx: LapackInt;
                work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_chbgst_work", dynlib: libName.}
proc chbgv*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; ka: LapackInt;
           kb: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
           bb: ptr LapackComplexFloat; ldbb: LapackInt; w: ptr cfloat;
           z: ptr LapackComplexFloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_chbgv", dynlib: libName.}
proc chbgvWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; ka: LapackInt;
               kb: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
               bb: ptr LapackComplexFloat; ldbb: LapackInt; w: ptr cfloat;
               z: ptr LapackComplexFloat; ldz: LapackInt;
               work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_chbgv_work", dynlib: libName.}
proc chbgvd*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; ka: LapackInt;
            kb: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
            bb: ptr LapackComplexFloat; ldbb: LapackInt; w: ptr cfloat;
            z: ptr LapackComplexFloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_chbgvd", dynlib: libName.}
proc chbgvdWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; ka: LapackInt;
                kb: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
                bb: ptr LapackComplexFloat; ldbb: LapackInt; w: ptr cfloat;
                z: ptr LapackComplexFloat; ldz: LapackInt;
                work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat;
                lrwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_chbgvd_work", dynlib: libName.}
proc chbgvx*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
            ka: LapackInt; kb: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
            bb: ptr LapackComplexFloat; ldbb: LapackInt; q: ptr LapackComplexFloat;
            ldq: LapackInt; vl: cfloat; vu: cfloat; il: LapackInt; iu: LapackInt;
            abstol: cfloat; m: ptr LapackInt; w: ptr cfloat; z: ptr LapackComplexFloat;
            ldz: LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_chbgvx", dynlib: libName.}
proc chbgvxWork*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                ka: LapackInt; kb: LapackInt; ab: ptr LapackComplexFloat;
                ldab: LapackInt; bb: ptr LapackComplexFloat; ldbb: LapackInt;
                q: ptr LapackComplexFloat; ldq: LapackInt; vl: cfloat; vu: cfloat;
                il: LapackInt; iu: LapackInt; abstol: cfloat; m: ptr LapackInt;
                w: ptr cfloat; z: ptr LapackComplexFloat; ldz: LapackInt;
                work: ptr LapackComplexFloat; rwork: ptr cfloat; iwork: ptr LapackInt;
                ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_chbgvx_work",
    dynlib: libName.}
proc chbtrd*(matrixLayout: cint; vect: char; uplo: char; n: LapackInt; kd: LapackInt;
            ab: ptr LapackComplexFloat; ldab: LapackInt; d: ptr cfloat; e: ptr cfloat;
            q: ptr LapackComplexFloat; ldq: LapackInt): LapackInt {.
    importc: "LAPACKE_chbtrd", dynlib: libName.}
proc chbtrdWork*(matrixLayout: cint; vect: char; uplo: char; n: LapackInt; kd: LapackInt;
                ab: ptr LapackComplexFloat; ldab: LapackInt; d: ptr cfloat;
                e: ptr cfloat; q: ptr LapackComplexFloat; ldq: LapackInt;
                work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_chbtrd_work", dynlib: libName.}
proc checon*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; ipiv: ptr LapackInt; anorm: cfloat; rcond: ptr cfloat): LapackInt {.
    importc: "LAPACKE_checon", dynlib: libName.}
proc checonWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                anorm: cfloat; rcond: ptr cfloat; work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_checon_work", dynlib: libName.}
proc cheequb*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
             lda: LapackInt; s: ptr cfloat; scond: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cheequb", dynlib: libName.}
proc cheequbWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt; s: ptr cfloat;
                 scond: ptr cfloat; amax: ptr cfloat; work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cheequb_work", dynlib: libName.}
proc cheev*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
           a: ptr LapackComplexFloat; lda: LapackInt; w: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cheev", dynlib: libName.}
proc cheevWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
               a: ptr LapackComplexFloat; lda: LapackInt; w: ptr cfloat;
               work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cheev_work", dynlib: libName.}
proc cheevd*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; w: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cheevd", dynlib: libName.}
proc cheevdWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; w: ptr cfloat;
                work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat;
                lrwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cheevd_work", dynlib: libName.}
proc cheevr*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; vl: cfloat; vu: cfloat;
            il: LapackInt; iu: LapackInt; abstol: cfloat; m: ptr LapackInt;
            w: ptr cfloat; z: ptr LapackComplexFloat; ldz: LapackInt;
            isuppz: ptr LapackInt): LapackInt {.importc: "LAPACKE_cheevr",
    dynlib: libName.}
proc cheevrWork*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; vl: cfloat; vu: cfloat;
                il: LapackInt; iu: LapackInt; abstol: cfloat; m: ptr LapackInt;
                w: ptr cfloat; z: ptr LapackComplexFloat; ldz: LapackInt;
                isuppz: ptr LapackInt; work: ptr LapackComplexFloat; lwork: LapackInt;
                rwork: ptr cfloat; lrwork: LapackInt; iwork: ptr LapackInt;
                liwork: LapackInt): LapackInt {.importc: "LAPACKE_cheevr_work",
    dynlib: libName.}
proc cheevx*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; vl: cfloat; vu: cfloat;
            il: LapackInt; iu: LapackInt; abstol: cfloat; m: ptr LapackInt;
            w: ptr cfloat; z: ptr LapackComplexFloat; ldz: LapackInt;
            ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_cheevx",
    dynlib: libName.}
proc cheevxWork*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; vl: cfloat; vu: cfloat;
                il: LapackInt; iu: LapackInt; abstol: cfloat; m: ptr LapackInt;
                w: ptr cfloat; z: ptr LapackComplexFloat; ldz: LapackInt;
                work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat;
                iwork: ptr LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_cheevx_work", dynlib: libName.}
proc chegst*(matrixLayout: cint; itype: LapackInt; uplo: char; n: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
            ldb: LapackInt): LapackInt {.importc: "LAPACKE_chegst", dynlib: libName.}
proc chegstWork*(matrixLayout: cint; itype: LapackInt; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
                ldb: LapackInt): LapackInt {.importc: "LAPACKE_chegst_work",
    dynlib: libName.}
proc chegv*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char; n: LapackInt;
           a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
           ldb: LapackInt; w: ptr cfloat): LapackInt {.importc: "LAPACKE_chegv",
    dynlib: libName.}
proc chegvWork*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
               n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
               b: ptr LapackComplexFloat; ldb: LapackInt; w: ptr cfloat;
               work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_chegv_work", dynlib: libName.}
proc chegvd*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char; n: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
            ldb: LapackInt; w: ptr cfloat): LapackInt {.importc: "LAPACKE_chegvd",
    dynlib: libName.}
proc chegvdWork*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
                n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt; w: ptr cfloat;
                work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat;
                lrwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_chegvd_work", dynlib: libName.}
proc chegvx*(matrixLayout: cint; itype: LapackInt; jobz: char; range: char; uplo: char;
            n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            b: ptr LapackComplexFloat; ldb: LapackInt; vl: cfloat; vu: cfloat;
            il: LapackInt; iu: LapackInt; abstol: cfloat; m: ptr LapackInt;
            w: ptr cfloat; z: ptr LapackComplexFloat; ldz: LapackInt;
            ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_chegvx",
    dynlib: libName.}
proc chegvxWork*(matrixLayout: cint; itype: LapackInt; jobz: char; range: char;
                uplo: char; n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt; vl: cfloat; vu: cfloat;
                il: LapackInt; iu: LapackInt; abstol: cfloat; m: ptr LapackInt;
                w: ptr cfloat; z: ptr LapackComplexFloat; ldz: LapackInt;
                work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat;
                iwork: ptr LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_chegvx_work", dynlib: libName.}
proc cherfs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; af: ptr LapackComplexFloat;
            ldaf: LapackInt; ipiv: ptr LapackInt; b: ptr LapackComplexFloat;
            ldb: LapackInt; x: ptr LapackComplexFloat; ldx: LapackInt;
            ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.importc: "LAPACKE_cherfs",
    dynlib: libName.}
proc cherfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt;
                af: ptr LapackComplexFloat; ldaf: LapackInt; ipiv: ptr LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
                ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat;
                work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cherfs_work", dynlib: libName.}
proc cherfsx*(matrixLayout: cint; uplo: char; equed: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; af: ptr LapackComplexFloat;
             ldaf: LapackInt; ipiv: ptr LapackInt; s: ptr cfloat;
             b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
             ldx: LapackInt; rcond: ptr cfloat; berr: ptr cfloat; nErrBnds: LapackInt;
             errBndsNorm: ptr cfloat; errBndsComp: ptr cfloat; nparams: LapackInt;
             params: ptr cfloat): LapackInt {.importc: "LAPACKE_cherfsx",
    dynlib: libName.}
proc cherfsxWork*(matrixLayout: cint; uplo: char; equed: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                 af: ptr LapackComplexFloat; ldaf: LapackInt; ipiv: ptr LapackInt;
                 s: ptr cfloat; b: ptr LapackComplexFloat; ldb: LapackInt;
                 x: ptr LapackComplexFloat; ldx: LapackInt; rcond: ptr cfloat;
                 berr: ptr cfloat; nErrBnds: LapackInt; errBndsNorm: ptr cfloat;
                 errBndsComp: ptr cfloat; nparams: LapackInt; params: ptr cfloat;
                 work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cherfsx_work", dynlib: libName.}
proc chesv*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
           a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
           b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_chesv", dynlib: libName.}
proc chesvWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
               b: ptr LapackComplexFloat; ldb: LapackInt;
               work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_chesv_work", dynlib: libName.}
proc chesvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; af: ptr LapackComplexFloat;
            ldaf: LapackInt; ipiv: ptr LapackInt; b: ptr LapackComplexFloat;
            ldb: LapackInt; x: ptr LapackComplexFloat; ldx: LapackInt;
            rcond: ptr cfloat; ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.
    importc: "LAPACKE_chesvx", dynlib: libName.}
proc chesvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                af: ptr LapackComplexFloat; ldaf: LapackInt; ipiv: ptr LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
                ldx: LapackInt; rcond: ptr cfloat; ferr: ptr cfloat; berr: ptr cfloat;
                work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_chesvx_work", dynlib: libName.}
proc chesvxx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; af: ptr LapackComplexFloat;
             ldaf: LapackInt; ipiv: ptr LapackInt; equed: cstring; s: ptr cfloat;
             b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
             ldx: LapackInt; rcond: ptr cfloat; rpvgrw: ptr cfloat; berr: ptr cfloat;
             nErrBnds: LapackInt; errBndsNorm: ptr cfloat; errBndsComp: ptr cfloat;
             nparams: LapackInt; params: ptr cfloat): LapackInt {.
    importc: "LAPACKE_chesvxx", dynlib: libName.}
proc chesvxxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                 af: ptr LapackComplexFloat; ldaf: LapackInt; ipiv: ptr LapackInt;
                 equed: cstring; s: ptr cfloat; b: ptr LapackComplexFloat;
                 ldb: LapackInt; x: ptr LapackComplexFloat; ldx: LapackInt;
                 rcond: ptr cfloat; rpvgrw: ptr cfloat; berr: ptr cfloat;
                 nErrBnds: LapackInt; errBndsNorm: ptr cfloat;
                 errBndsComp: ptr cfloat; nparams: LapackInt; params: ptr cfloat;
                 work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_chesvxx_work", dynlib: libName.}
proc cheswapr*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
              lda: LapackInt; i1: LapackInt; i2: LapackInt): LapackInt {.
    importc: "LAPACKE_cheswapr", dynlib: libName.}
proc cheswaprWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                  a: ptr LapackComplexFloat; lda: LapackInt; i1: LapackInt;
                  i2: LapackInt): LapackInt {.importc: "LAPACKE_cheswapr_work",
    dynlib: libName.}
proc chetrd*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; d: ptr cfloat; e: ptr cfloat; tau: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_chetrd", dynlib: libName.}
proc chetrdWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; d: ptr cfloat; e: ptr cfloat;
                tau: ptr LapackComplexFloat; work: ptr LapackComplexFloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_chetrd_work",
    dynlib: libName.}
proc chetrf*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_chetrf", dynlib: libName.}
proc chetrfRook*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_chetrf_rook", dynlib: libName.}
proc chetrfRookWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                    a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                    work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_chetrf_rook_work", dynlib: libName.}
proc chetrfWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_chetrf_work", dynlib: libName.}
proc chetri*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_chetri", dynlib: libName.}
proc chetri2*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
             lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_chetri2", dynlib: libName.}
proc chetri2Work*(matrixLayout: cint; uplo: char; n: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                 work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_chetri2_work", dynlib: libName.}
proc chetri2x*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
              lda: LapackInt; ipiv: ptr LapackInt; nb: LapackInt): LapackInt {.
    importc: "LAPACKE_chetri2x", dynlib: libName.}
proc chetri2xWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                  a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                  work: ptr LapackComplexFloat; nb: LapackInt): LapackInt {.
    importc: "LAPACKE_chetri2x_work", dynlib: libName.}
proc chetriWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_chetri_work", dynlib: libName.}
proc chetrs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
            b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_chetrs", dynlib: libName.}
proc chetrs2*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
             b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_chetrs2", dynlib: libName.}
proc chetrs2Work*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                 b: ptr LapackComplexFloat; ldb: LapackInt;
                 work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_chetrs2_work", dynlib: libName.}
proc chetrsRook*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_chetrs_rook", dynlib: libName.}
proc chetrsRookWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                    a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                    b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_chetrs_rook_work", dynlib: libName.}
proc chetrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_chetrs_work", dynlib: libName.}
proc chfrk*(matrixLayout: cint; transr: char; uplo: char; trans: char; n: LapackInt;
           k: LapackInt; alpha: cfloat; a: ptr LapackComplexFloat; lda: LapackInt;
           beta: cfloat; c: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_chfrk", dynlib: libName.}
proc chfrkWork*(matrixLayout: cint; transr: char; uplo: char; trans: char; n: LapackInt;
               k: LapackInt; alpha: cfloat; a: ptr LapackComplexFloat; lda: LapackInt;
               beta: cfloat; c: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_chfrk_work", dynlib: libName.}
proc chgeqz*(matrixLayout: cint; job: char; compq: char; compz: char; n: LapackInt;
            ilo: LapackInt; ihi: LapackInt; h: ptr LapackComplexFloat; ldh: LapackInt;
            t: ptr LapackComplexFloat; ldt: LapackInt; alpha: ptr LapackComplexFloat;
            beta: ptr LapackComplexFloat; q: ptr LapackComplexFloat; ldq: LapackInt;
            z: ptr LapackComplexFloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_chgeqz", dynlib: libName.}
proc chgeqzWork*(matrixLayout: cint; job: char; compq: char; compz: char; n: LapackInt;
                ilo: LapackInt; ihi: LapackInt; h: ptr LapackComplexFloat;
                ldh: LapackInt; t: ptr LapackComplexFloat; ldt: LapackInt;
                alpha: ptr LapackComplexFloat; beta: ptr LapackComplexFloat;
                q: ptr LapackComplexFloat; ldq: LapackInt; z: ptr LapackComplexFloat;
                ldz: LapackInt; work: ptr LapackComplexFloat; lwork: LapackInt;
                rwork: ptr cfloat): LapackInt {.importc: "LAPACKE_chgeqz_work",
    dynlib: libName.}
proc chpcon*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexFloat;
            ipiv: ptr LapackInt; anorm: cfloat; rcond: ptr cfloat): LapackInt {.
    importc: "LAPACKE_chpcon", dynlib: libName.}
proc chpconWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexFloat; ipiv: ptr LapackInt; anorm: cfloat;
                rcond: ptr cfloat; work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_chpcon_work", dynlib: libName.}
proc chpev*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
           ap: ptr LapackComplexFloat; w: ptr cfloat; z: ptr LapackComplexFloat;
           ldz: LapackInt): LapackInt {.importc: "LAPACKE_chpev", dynlib: libName.}
proc chpevWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
               ap: ptr LapackComplexFloat; w: ptr cfloat; z: ptr LapackComplexFloat;
               ldz: LapackInt; work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_chpev_work", dynlib: libName.}
proc chpevd*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
            ap: ptr LapackComplexFloat; w: ptr cfloat; z: ptr LapackComplexFloat;
            ldz: LapackInt): LapackInt {.importc: "LAPACKE_chpevd", dynlib: libName.}
proc chpevdWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                ap: ptr LapackComplexFloat; w: ptr cfloat; z: ptr LapackComplexFloat;
                ldz: LapackInt; work: ptr LapackComplexFloat; lwork: LapackInt;
                rwork: ptr cfloat; lrwork: LapackInt; iwork: ptr LapackInt;
                liwork: LapackInt): LapackInt {.importc: "LAPACKE_chpevd_work",
    dynlib: libName.}
proc chpevx*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
            ap: ptr LapackComplexFloat; vl: cfloat; vu: cfloat; il: LapackInt;
            iu: LapackInt; abstol: cfloat; m: ptr LapackInt; w: ptr cfloat;
            z: ptr LapackComplexFloat; ldz: LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_chpevx", dynlib: libName.}
proc chpevxWork*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                ap: ptr LapackComplexFloat; vl: cfloat; vu: cfloat; il: LapackInt;
                iu: LapackInt; abstol: cfloat; m: ptr LapackInt; w: ptr cfloat;
                z: ptr LapackComplexFloat; ldz: LapackInt;
                work: ptr LapackComplexFloat; rwork: ptr cfloat; iwork: ptr LapackInt;
                ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_chpevx_work",
    dynlib: libName.}
proc chpgst*(matrixLayout: cint; itype: LapackInt; uplo: char; n: LapackInt;
            ap: ptr LapackComplexFloat; bp: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_chpgst", dynlib: libName.}
proc chpgstWork*(matrixLayout: cint; itype: LapackInt; uplo: char; n: LapackInt;
                ap: ptr LapackComplexFloat; bp: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_chpgst_work", dynlib: libName.}
proc chpgv*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char; n: LapackInt;
           ap: ptr LapackComplexFloat; bp: ptr LapackComplexFloat; w: ptr cfloat;
           z: ptr LapackComplexFloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_chpgv", dynlib: libName.}
proc chpgvWork*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
               n: LapackInt; ap: ptr LapackComplexFloat; bp: ptr LapackComplexFloat;
               w: ptr cfloat; z: ptr LapackComplexFloat; ldz: LapackInt;
               work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_chpgv_work", dynlib: libName.}
proc chpgvd*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char; n: LapackInt;
            ap: ptr LapackComplexFloat; bp: ptr LapackComplexFloat; w: ptr cfloat;
            z: ptr LapackComplexFloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_chpgvd", dynlib: libName.}
proc chpgvdWork*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
                n: LapackInt; ap: ptr LapackComplexFloat; bp: ptr LapackComplexFloat;
                w: ptr cfloat; z: ptr LapackComplexFloat; ldz: LapackInt;
                work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat;
                lrwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_chpgvd_work", dynlib: libName.}
proc chpgvx*(matrixLayout: cint; itype: LapackInt; jobz: char; range: char; uplo: char;
            n: LapackInt; ap: ptr LapackComplexFloat; bp: ptr LapackComplexFloat;
            vl: cfloat; vu: cfloat; il: LapackInt; iu: LapackInt; abstol: cfloat;
            m: ptr LapackInt; w: ptr cfloat; z: ptr LapackComplexFloat; ldz: LapackInt;
            ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_chpgvx",
    dynlib: libName.}
proc chpgvxWork*(matrixLayout: cint; itype: LapackInt; jobz: char; range: char;
                uplo: char; n: LapackInt; ap: ptr LapackComplexFloat;
                bp: ptr LapackComplexFloat; vl: cfloat; vu: cfloat; il: LapackInt;
                iu: LapackInt; abstol: cfloat; m: ptr LapackInt; w: ptr cfloat;
                z: ptr LapackComplexFloat; ldz: LapackInt;
                work: ptr LapackComplexFloat; rwork: ptr cfloat; iwork: ptr LapackInt;
                ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_chpgvx_work",
    dynlib: libName.}
proc chprfs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr LapackComplexFloat; afp: ptr LapackComplexFloat;
            ipiv: ptr LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt;
            x: ptr LapackComplexFloat; ldx: LapackInt; ferr: ptr cfloat;
            berr: ptr cfloat): LapackInt {.importc: "LAPACKE_chprfs", dynlib: libName.}
proc chprfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                ap: ptr LapackComplexFloat; afp: ptr LapackComplexFloat;
                ipiv: ptr LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt;
                x: ptr LapackComplexFloat; ldx: LapackInt; ferr: ptr cfloat;
                berr: ptr cfloat; work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_chprfs_work", dynlib: libName.}
proc chpsv*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
           ap: ptr LapackComplexFloat; ipiv: ptr LapackInt; b: ptr LapackComplexFloat;
           ldb: LapackInt): LapackInt {.importc: "LAPACKE_chpsv", dynlib: libName.}
proc chpsvWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               ap: ptr LapackComplexFloat; ipiv: ptr LapackInt;
               b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_chpsv_work", dynlib: libName.}
proc chpsvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr LapackComplexFloat; afp: ptr LapackComplexFloat;
            ipiv: ptr LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt;
            x: ptr LapackComplexFloat; ldx: LapackInt; rcond: ptr cfloat;
            ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.importc: "LAPACKE_chpsvx",
    dynlib: libName.}
proc chpsvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; ap: ptr LapackComplexFloat;
                afp: ptr LapackComplexFloat; ipiv: ptr LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
                ldx: LapackInt; rcond: ptr cfloat; ferr: ptr cfloat; berr: ptr cfloat;
                work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_chpsvx_work", dynlib: libName.}
proc chptrd*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexFloat;
            d: ptr cfloat; e: ptr cfloat; tau: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_chptrd", dynlib: libName.}
proc chptrdWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexFloat; d: ptr cfloat; e: ptr cfloat;
                tau: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_chptrd_work", dynlib: libName.}
proc chptrf*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexFloat;
            ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_chptrf",
    dynlib: libName.}
proc chptrfWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexFloat; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_chptrf_work", dynlib: libName.}
proc chptri*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexFloat;
            ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_chptri",
    dynlib: libName.}
proc chptriWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexFloat; ipiv: ptr LapackInt;
                work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_chptri_work", dynlib: libName.}
proc chptrs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr LapackComplexFloat; ipiv: ptr LapackInt;
            b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_chptrs", dynlib: libName.}
proc chptrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                ap: ptr LapackComplexFloat; ipiv: ptr LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_chptrs_work", dynlib: libName.}
proc chsein*(matrixLayout: cint; job: char; eigsrc: char; initv: char;
            select: ptr LapackLogical; n: LapackInt; h: ptr LapackComplexFloat;
            ldh: LapackInt; w: ptr LapackComplexFloat; vl: ptr LapackComplexFloat;
            ldvl: LapackInt; vr: ptr LapackComplexFloat; ldvr: LapackInt;
            mm: LapackInt; m: ptr LapackInt; ifaill: ptr LapackInt;
            ifailr: ptr LapackInt): LapackInt {.importc: "LAPACKE_chsein",
    dynlib: libName.}
proc chseinWork*(matrixLayout: cint; job: char; eigsrc: char; initv: char;
                select: ptr LapackLogical; n: LapackInt; h: ptr LapackComplexFloat;
                ldh: LapackInt; w: ptr LapackComplexFloat;
                vl: ptr LapackComplexFloat; ldvl: LapackInt;
                vr: ptr LapackComplexFloat; ldvr: LapackInt; mm: LapackInt;
                m: ptr LapackInt; work: ptr LapackComplexFloat; rwork: ptr cfloat;
                ifaill: ptr LapackInt; ifailr: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_chsein_work", dynlib: libName.}
proc chseqr*(matrixLayout: cint; job: char; compz: char; n: LapackInt; ilo: LapackInt;
            ihi: LapackInt; h: ptr LapackComplexFloat; ldh: LapackInt;
            w: ptr LapackComplexFloat; z: ptr LapackComplexFloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_chseqr", dynlib: libName.}
proc chseqrWork*(matrixLayout: cint; job: char; compz: char; n: LapackInt;
                ilo: LapackInt; ihi: LapackInt; h: ptr LapackComplexFloat;
                ldh: LapackInt; w: ptr LapackComplexFloat; z: ptr LapackComplexFloat;
                ldz: LapackInt; work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_chseqr_work", dynlib: libName.}
proc clacgv*(n: LapackInt; x: ptr LapackComplexFloat; incx: LapackInt): LapackInt {.
    importc: "LAPACKE_clacgv", dynlib: libName.}
proc clacgvWork*(n: LapackInt; x: ptr LapackComplexFloat; incx: LapackInt): LapackInt {.
    importc: "LAPACKE_clacgv_work", dynlib: libName.}
proc clacn2*(n: LapackInt; v: ptr LapackComplexFloat; x: ptr LapackComplexFloat;
            est: ptr cfloat; kase: ptr LapackInt; isave: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_clacn2", dynlib: libName.}
proc clacn2Work*(n: LapackInt; v: ptr LapackComplexFloat; x: ptr LapackComplexFloat;
                est: ptr cfloat; kase: ptr LapackInt; isave: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_clacn2_work", dynlib: libName.}
proc clacp2*(matrixLayout: cint; uplo: char; m: LapackInt; n: LapackInt; a: ptr cfloat;
            lda: LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_clacp2", dynlib: libName.}
proc clacp2Work*(matrixLayout: cint; uplo: char; m: LapackInt; n: LapackInt;
                a: ptr cfloat; lda: LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_clacp2_work", dynlib: libName.}
proc clacpy*(matrixLayout: cint; uplo: char; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
            ldb: LapackInt): LapackInt {.importc: "LAPACKE_clacpy", dynlib: libName.}
proc clacpyWork*(matrixLayout: cint; uplo: char; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
                ldb: LapackInt): LapackInt {.importc: "LAPACKE_clacpy_work",
    dynlib: libName.}
proc clag2z*(matrixLayout: cint; m: LapackInt; n: LapackInt;
            sa: ptr LapackComplexFloat; ldsa: LapackInt; a: ptr LapackComplexDouble;
            lda: LapackInt): LapackInt {.importc: "LAPACKE_clag2z", dynlib: libName.}
proc clag2zWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                sa: ptr LapackComplexFloat; ldsa: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_clag2z_work", dynlib: libName.}
proc clagge*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
            ku: LapackInt; d: ptr cfloat; a: ptr LapackComplexFloat; lda: LapackInt;
            iseed: ptr LapackInt): LapackInt {.importc: "LAPACKE_clagge",
    dynlib: libName.}
proc claggeWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
                ku: LapackInt; d: ptr cfloat; a: ptr LapackComplexFloat; lda: LapackInt;
                iseed: ptr LapackInt; work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_clagge_work", dynlib: libName.}
proc claghe*(matrixLayout: cint; n: LapackInt; k: LapackInt; d: ptr cfloat;
            a: ptr LapackComplexFloat; lda: LapackInt; iseed: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_claghe", dynlib: libName.}
proc clagheWork*(matrixLayout: cint; n: LapackInt; k: LapackInt; d: ptr cfloat;
                a: ptr LapackComplexFloat; lda: LapackInt; iseed: ptr LapackInt;
                work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_claghe_work", dynlib: libName.}
proc clagsy*(matrixLayout: cint; n: LapackInt; k: LapackInt; d: ptr cfloat;
            a: ptr LapackComplexFloat; lda: LapackInt; iseed: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_clagsy", dynlib: libName.}
proc clagsyWork*(matrixLayout: cint; n: LapackInt; k: LapackInt; d: ptr cfloat;
                a: ptr LapackComplexFloat; lda: LapackInt; iseed: ptr LapackInt;
                work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_clagsy_work", dynlib: libName.}
proc clapmr*(matrixLayout: cint; forwrd: LapackLogical; m: LapackInt; n: LapackInt;
            x: ptr LapackComplexFloat; ldx: LapackInt; k: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_clapmr", dynlib: libName.}
proc clapmrWork*(matrixLayout: cint; forwrd: LapackLogical; m: LapackInt; n: LapackInt;
                x: ptr LapackComplexFloat; ldx: LapackInt; k: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_clapmr_work", dynlib: libName.}
proc clapmt*(matrixLayout: cint; forwrd: LapackLogical; m: LapackInt; n: LapackInt;
            x: ptr LapackComplexFloat; ldx: LapackInt; k: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_clapmt", dynlib: libName.}
proc clapmtWork*(matrixLayout: cint; forwrd: LapackLogical; m: LapackInt; n: LapackInt;
                x: ptr LapackComplexFloat; ldx: LapackInt; k: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_clapmt_work", dynlib: libName.}
proc clarfb*(matrixLayout: cint; side: char; trans: char; direct: char; storev: char;
            m: LapackInt; n: LapackInt; k: LapackInt; v: ptr LapackComplexFloat;
            ldv: LapackInt; t: ptr LapackComplexFloat; ldt: LapackInt;
            c: ptr LapackComplexFloat; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_clarfb", dynlib: libName.}
proc clarfbWork*(matrixLayout: cint; side: char; trans: char; direct: char; storev: char;
                m: LapackInt; n: LapackInt; k: LapackInt; v: ptr LapackComplexFloat;
                ldv: LapackInt; t: ptr LapackComplexFloat; ldt: LapackInt;
                c: ptr LapackComplexFloat; ldc: LapackInt;
                work: ptr LapackComplexFloat; ldwork: LapackInt): LapackInt {.
    importc: "LAPACKE_clarfb_work", dynlib: libName.}
proc clarfg*(n: LapackInt; alpha: ptr LapackComplexFloat; x: ptr LapackComplexFloat;
            incx: LapackInt; tau: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_clarfg", dynlib: libName.}
proc clarfgWork*(n: LapackInt; alpha: ptr LapackComplexFloat;
                x: ptr LapackComplexFloat; incx: LapackInt;
                tau: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_clarfg_work", dynlib: libName.}
proc clarft*(matrixLayout: cint; direct: char; storev: char; n: LapackInt; k: LapackInt;
            v: ptr LapackComplexFloat; ldv: LapackInt; tau: ptr LapackComplexFloat;
            t: ptr LapackComplexFloat; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_clarft", dynlib: libName.}
proc clarftWork*(matrixLayout: cint; direct: char; storev: char; n: LapackInt;
                k: LapackInt; v: ptr LapackComplexFloat; ldv: LapackInt;
                tau: ptr LapackComplexFloat; t: ptr LapackComplexFloat; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_clarft_work", dynlib: libName.}
proc clarfx*(matrixLayout: cint; side: char; m: LapackInt; n: LapackInt;
            v: ptr LapackComplexFloat; tau: LapackComplexFloat;
            c: ptr LapackComplexFloat; ldc: LapackInt; work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_clarfx", dynlib: libName.}
proc clarfxWork*(matrixLayout: cint; side: char; m: LapackInt; n: LapackInt;
                v: ptr LapackComplexFloat; tau: LapackComplexFloat;
                c: ptr LapackComplexFloat; ldc: LapackInt;
                work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_clarfx_work", dynlib: libName.}
proc clarnv*(idist: LapackInt; iseed: ptr LapackInt; n: LapackInt;
            x: ptr LapackComplexFloat): LapackInt {.importc: "LAPACKE_clarnv",
    dynlib: libName.}
proc clarnvWork*(idist: LapackInt; iseed: ptr LapackInt; n: LapackInt;
                x: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_clarnv_work", dynlib: libName.}
proc clascl*(matrixLayout: cint; `type`: char; kl: LapackInt; ku: LapackInt;
            cfrom: cfloat; cto: cfloat; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_clascl", dynlib: libName.}
proc clasclWork*(matrixLayout: cint; `type`: char; kl: LapackInt; ku: LapackInt;
                cfrom: cfloat; cto: cfloat; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_clascl_work", dynlib: libName.}
proc claset*(matrixLayout: cint; uplo: char; m: LapackInt; n: LapackInt;
            alpha: LapackComplexFloat; beta: LapackComplexFloat;
            a: ptr LapackComplexFloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_claset", dynlib: libName.}
proc clasetWork*(matrixLayout: cint; uplo: char; m: LapackInt; n: LapackInt;
                alpha: LapackComplexFloat; beta: LapackComplexFloat;
                a: ptr LapackComplexFloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_claset_work", dynlib: libName.}
proc claswp*(matrixLayout: cint; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; k1: LapackInt; k2: LapackInt; ipiv: ptr LapackInt;
            incx: LapackInt): LapackInt {.importc: "LAPACKE_claswp", dynlib: libName.}
proc claswpWork*(matrixLayout: cint; n: LapackInt; a: ptr LapackComplexFloat;
                lda: LapackInt; k1: LapackInt; k2: LapackInt; ipiv: ptr LapackInt;
                incx: LapackInt): LapackInt {.importc: "LAPACKE_claswp_work",
    dynlib: libName.}
proc clatms*(matrixLayout: cint; m: LapackInt; n: LapackInt; dist: char;
            iseed: ptr LapackInt; sym: char; d: ptr cfloat; mode: LapackInt; cond: cfloat;
            dmax: cfloat; kl: LapackInt; ku: LapackInt; pack: char;
            a: ptr LapackComplexFloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_clatms", dynlib: libName.}
proc clatmsWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; dist: char;
                iseed: ptr LapackInt; sym: char; d: ptr cfloat; mode: LapackInt;
                cond: cfloat; dmax: cfloat; kl: LapackInt; ku: LapackInt; pack: char;
                a: ptr LapackComplexFloat; lda: LapackInt;
                work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_clatms_work", dynlib: libName.}
proc clauum*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt): LapackInt {.importc: "LAPACKE_clauum", dynlib: libName.}
proc clauumWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_clauum_work", dynlib: libName.}
proc cpbcon*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
            ab: ptr LapackComplexFloat; ldab: LapackInt; anorm: cfloat;
            rcond: ptr cfloat): LapackInt {.importc: "LAPACKE_cpbcon", dynlib: libName.}
proc cpbconWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
                ab: ptr LapackComplexFloat; ldab: LapackInt; anorm: cfloat;
                rcond: ptr cfloat; work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cpbcon_work", dynlib: libName.}
proc cpbequ*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
            ab: ptr LapackComplexFloat; ldab: LapackInt; s: ptr cfloat;
            scond: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cpbequ", dynlib: libName.}
proc cpbequWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
                ab: ptr LapackComplexFloat; ldab: LapackInt; s: ptr cfloat;
                scond: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cpbequ_work", dynlib: libName.}
proc cpbrfs*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
            nrhs: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
            afb: ptr LapackComplexFloat; ldafb: LapackInt; b: ptr LapackComplexFloat;
            ldb: LapackInt; x: ptr LapackComplexFloat; ldx: LapackInt;
            ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.importc: "LAPACKE_cpbrfs",
    dynlib: libName.}
proc cpbrfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
                nrhs: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
                afb: ptr LapackComplexFloat; ldafb: LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
                ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat;
                work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cpbrfs_work", dynlib: libName.}
proc cpbstf*(matrixLayout: cint; uplo: char; n: LapackInt; kb: LapackInt;
            bb: ptr LapackComplexFloat; ldbb: LapackInt): LapackInt {.
    importc: "LAPACKE_cpbstf", dynlib: libName.}
proc cpbstfWork*(matrixLayout: cint; uplo: char; n: LapackInt; kb: LapackInt;
                bb: ptr LapackComplexFloat; ldbb: LapackInt): LapackInt {.
    importc: "LAPACKE_cpbstf_work", dynlib: libName.}
proc cpbsv*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt; nrhs: LapackInt;
           ab: ptr LapackComplexFloat; ldab: LapackInt; b: ptr LapackComplexFloat;
           ldb: LapackInt): LapackInt {.importc: "LAPACKE_cpbsv", dynlib: libName.}
proc cpbsvWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
               nrhs: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
               b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cpbsv_work", dynlib: libName.}
proc cpbsvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; kd: LapackInt;
            nrhs: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
            afb: ptr LapackComplexFloat; ldafb: LapackInt; equed: cstring;
            s: ptr cfloat; b: ptr LapackComplexFloat; ldb: LapackInt;
            x: ptr LapackComplexFloat; ldx: LapackInt; rcond: ptr cfloat;
            ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.importc: "LAPACKE_cpbsvx",
    dynlib: libName.}
proc cpbsvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; kd: LapackInt;
                nrhs: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
                afb: ptr LapackComplexFloat; ldafb: LapackInt; equed: cstring;
                s: ptr cfloat; b: ptr LapackComplexFloat; ldb: LapackInt;
                x: ptr LapackComplexFloat; ldx: LapackInt; rcond: ptr cfloat;
                ferr: ptr cfloat; berr: ptr cfloat; work: ptr LapackComplexFloat;
                rwork: ptr cfloat): LapackInt {.importc: "LAPACKE_cpbsvx_work",
    dynlib: libName.}
proc cpbtrf*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
            ab: ptr LapackComplexFloat; ldab: LapackInt): LapackInt {.
    importc: "LAPACKE_cpbtrf", dynlib: libName.}
proc cpbtrfWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
                ab: ptr LapackComplexFloat; ldab: LapackInt): LapackInt {.
    importc: "LAPACKE_cpbtrf_work", dynlib: libName.}
proc cpbtrs*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
            nrhs: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
            b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cpbtrs", dynlib: libName.}
proc cpbtrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
                nrhs: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cpbtrs_work", dynlib: libName.}
proc cpftrf*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
            a: ptr LapackComplexFloat): LapackInt {.importc: "LAPACKE_cpftrf",
    dynlib: libName.}
proc cpftrfWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cpftrf_work", dynlib: libName.}
proc cpftri*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
            a: ptr LapackComplexFloat): LapackInt {.importc: "LAPACKE_cpftri",
    dynlib: libName.}
proc cpftriWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cpftri_work", dynlib: libName.}
proc cpftrs*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexFloat; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cpftrs", dynlib: libName.}
proc cpftrsWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; a: ptr LapackComplexFloat;
                b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cpftrs_work", dynlib: libName.}
proc cpocon*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; anorm: cfloat; rcond: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cpocon", dynlib: libName.}
proc cpoconWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; anorm: cfloat;
                rcond: ptr cfloat; work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cpocon_work", dynlib: libName.}
proc cpoequ*(matrixLayout: cint; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; s: ptr cfloat; scond: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cpoequ", dynlib: libName.}
proc cpoequWork*(matrixLayout: cint; n: LapackInt; a: ptr LapackComplexFloat;
                lda: LapackInt; s: ptr cfloat; scond: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cpoequ_work", dynlib: libName.}
proc cpoequb*(matrixLayout: cint; n: LapackInt; a: ptr LapackComplexFloat;
             lda: LapackInt; s: ptr cfloat; scond: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cpoequb", dynlib: libName.}
proc cpoequbWork*(matrixLayout: cint; n: LapackInt; a: ptr LapackComplexFloat;
                 lda: LapackInt; s: ptr cfloat; scond: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cpoequb_work", dynlib: libName.}
proc cporfs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; af: ptr LapackComplexFloat;
            ldaf: LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt;
            x: ptr LapackComplexFloat; ldx: LapackInt; ferr: ptr cfloat;
            berr: ptr cfloat): LapackInt {.importc: "LAPACKE_cporfs", dynlib: libName.}
proc cporfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt;
                af: ptr LapackComplexFloat; ldaf: LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
                ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat;
                work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cporfs_work", dynlib: libName.}
proc cporfsx*(matrixLayout: cint; uplo: char; equed: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; af: ptr LapackComplexFloat;
             ldaf: LapackInt; s: ptr cfloat; b: ptr LapackComplexFloat; ldb: LapackInt;
             x: ptr LapackComplexFloat; ldx: LapackInt; rcond: ptr cfloat;
             berr: ptr cfloat; nErrBnds: LapackInt; errBndsNorm: ptr cfloat;
             errBndsComp: ptr cfloat; nparams: LapackInt; params: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cporfsx", dynlib: libName.}
proc cporfsxWork*(matrixLayout: cint; uplo: char; equed: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                 af: ptr LapackComplexFloat; ldaf: LapackInt; s: ptr cfloat;
                 b: ptr LapackComplexFloat; ldb: LapackInt;
                 x: ptr LapackComplexFloat; ldx: LapackInt; rcond: ptr cfloat;
                 berr: ptr cfloat; nErrBnds: LapackInt; errBndsNorm: ptr cfloat;
                 errBndsComp: ptr cfloat; nparams: LapackInt; params: ptr cfloat;
                 work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cporfsx_work", dynlib: libName.}
proc cposv*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
           a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
           ldb: LapackInt): LapackInt {.importc: "LAPACKE_cposv", dynlib: libName.}
proc cposvWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
               ldb: LapackInt): LapackInt {.importc: "LAPACKE_cposv_work",
    dynlib: libName.}
proc cposvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; af: ptr LapackComplexFloat;
            ldaf: LapackInt; equed: cstring; s: ptr cfloat; b: ptr LapackComplexFloat;
            ldb: LapackInt; x: ptr LapackComplexFloat; ldx: LapackInt;
            rcond: ptr cfloat; ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cposvx", dynlib: libName.}
proc cposvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                af: ptr LapackComplexFloat; ldaf: LapackInt; equed: cstring;
                s: ptr cfloat; b: ptr LapackComplexFloat; ldb: LapackInt;
                x: ptr LapackComplexFloat; ldx: LapackInt; rcond: ptr cfloat;
                ferr: ptr cfloat; berr: ptr cfloat; work: ptr LapackComplexFloat;
                rwork: ptr cfloat): LapackInt {.importc: "LAPACKE_cposvx_work",
    dynlib: libName.}
proc cposvxx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; af: ptr LapackComplexFloat;
             ldaf: LapackInt; equed: cstring; s: ptr cfloat; b: ptr LapackComplexFloat;
             ldb: LapackInt; x: ptr LapackComplexFloat; ldx: LapackInt;
             rcond: ptr cfloat; rpvgrw: ptr cfloat; berr: ptr cfloat;
             nErrBnds: LapackInt; errBndsNorm: ptr cfloat; errBndsComp: ptr cfloat;
             nparams: LapackInt; params: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cposvxx", dynlib: libName.}
proc cposvxxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                 af: ptr LapackComplexFloat; ldaf: LapackInt; equed: cstring;
                 s: ptr cfloat; b: ptr LapackComplexFloat; ldb: LapackInt;
                 x: ptr LapackComplexFloat; ldx: LapackInt; rcond: ptr cfloat;
                 rpvgrw: ptr cfloat; berr: ptr cfloat; nErrBnds: LapackInt;
                 errBndsNorm: ptr cfloat; errBndsComp: ptr cfloat; nparams: LapackInt;
                 params: ptr cfloat; work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cposvxx_work", dynlib: libName.}
proc cpotrf*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt): LapackInt {.importc: "LAPACKE_cpotrf", dynlib: libName.}
proc cpotrf2*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
             lda: LapackInt): LapackInt {.importc: "LAPACKE_cpotrf2", dynlib: libName.}
proc cpotrf2Work*(matrixLayout: cint; uplo: char; n: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_cpotrf2_work", dynlib: libName.}
proc cpotrfWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_cpotrf_work", dynlib: libName.}
proc cpotri*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt): LapackInt {.importc: "LAPACKE_cpotri", dynlib: libName.}
proc cpotriWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_cpotri_work", dynlib: libName.}
proc cpotrs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
            ldb: LapackInt): LapackInt {.importc: "LAPACKE_cpotrs", dynlib: libName.}
proc cpotrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
                ldb: LapackInt): LapackInt {.importc: "LAPACKE_cpotrs_work",
    dynlib: libName.}
proc cppcon*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexFloat;
            anorm: cfloat; rcond: ptr cfloat): LapackInt {.importc: "LAPACKE_cppcon",
    dynlib: libName.}
proc cppconWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexFloat; anorm: cfloat; rcond: ptr cfloat;
                work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cppcon_work", dynlib: libName.}
proc cppequ*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexFloat;
            s: ptr cfloat; scond: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cppequ", dynlib: libName.}
proc cppequWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexFloat; s: ptr cfloat; scond: ptr cfloat;
                amax: ptr cfloat): LapackInt {.importc: "LAPACKE_cppequ_work",
    dynlib: libName.}
proc cpprfs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr LapackComplexFloat; afp: ptr LapackComplexFloat;
            b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
            ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cpprfs", dynlib: libName.}
proc cpprfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                ap: ptr LapackComplexFloat; afp: ptr LapackComplexFloat;
                b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
                ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat;
                work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cpprfs_work", dynlib: libName.}
proc cppsv*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
           ap: ptr LapackComplexFloat; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cppsv", dynlib: libName.}
proc cppsvWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               ap: ptr LapackComplexFloat; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cppsv_work", dynlib: libName.}
proc cppsvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr LapackComplexFloat; afp: ptr LapackComplexFloat; equed: cstring;
            s: ptr cfloat; b: ptr LapackComplexFloat; ldb: LapackInt;
            x: ptr LapackComplexFloat; ldx: LapackInt; rcond: ptr cfloat;
            ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.importc: "LAPACKE_cppsvx",
    dynlib: libName.}
proc cppsvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; ap: ptr LapackComplexFloat;
                afp: ptr LapackComplexFloat; equed: cstring; s: ptr cfloat;
                b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
                ldx: LapackInt; rcond: ptr cfloat; ferr: ptr cfloat; berr: ptr cfloat;
                work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cppsvx_work", dynlib: libName.}
proc cpptrf*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cpptrf", dynlib: libName.}
proc cpptrfWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cpptrf_work", dynlib: libName.}
proc cpptri*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cpptri", dynlib: libName.}
proc cpptriWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cpptri_work", dynlib: libName.}
proc cpptrs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr LapackComplexFloat; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cpptrs", dynlib: libName.}
proc cpptrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                ap: ptr LapackComplexFloat; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cpptrs_work", dynlib: libName.}
proc cpstrf*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; piv: ptr LapackInt; rank: ptr LapackInt; tol: cfloat): LapackInt {.
    importc: "LAPACKE_cpstrf", dynlib: libName.}
proc cpstrfWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; piv: ptr LapackInt;
                rank: ptr LapackInt; tol: cfloat; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cpstrf_work", dynlib: libName.}
proc cptcon*(n: LapackInt; d: ptr cfloat; e: ptr LapackComplexFloat; anorm: cfloat;
            rcond: ptr cfloat): LapackInt {.importc: "LAPACKE_cptcon", dynlib: libName.}
proc cptconWork*(n: LapackInt; d: ptr cfloat; e: ptr LapackComplexFloat; anorm: cfloat;
                rcond: ptr cfloat; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cptcon_work", dynlib: libName.}
proc cpteqr*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cfloat; e: ptr cfloat;
            z: ptr LapackComplexFloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_cpteqr", dynlib: libName.}
proc cpteqrWork*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cfloat;
                e: ptr cfloat; z: ptr LapackComplexFloat; ldz: LapackInt;
                work: ptr cfloat): LapackInt {.importc: "LAPACKE_cpteqr_work",
    dynlib: libName.}
proc cptrfs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            d: ptr cfloat; e: ptr LapackComplexFloat; df: ptr cfloat;
            ef: ptr LapackComplexFloat; b: ptr LapackComplexFloat; ldb: LapackInt;
            x: ptr LapackComplexFloat; ldx: LapackInt; ferr: ptr cfloat;
            berr: ptr cfloat): LapackInt {.importc: "LAPACKE_cptrfs", dynlib: libName.}
proc cptrfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                d: ptr cfloat; e: ptr LapackComplexFloat; df: ptr cfloat;
                ef: ptr LapackComplexFloat; b: ptr LapackComplexFloat; ldb: LapackInt;
                x: ptr LapackComplexFloat; ldx: LapackInt; ferr: ptr cfloat;
                berr: ptr cfloat; work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cptrfs_work", dynlib: libName.}
proc cptsv*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; d: ptr cfloat;
           e: ptr LapackComplexFloat; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cptsv", dynlib: libName.}
proc cptsvWork*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; d: ptr cfloat;
               e: ptr LapackComplexFloat; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cptsv_work", dynlib: libName.}
proc cptsvx*(matrixLayout: cint; fact: char; n: LapackInt; nrhs: LapackInt;
            d: ptr cfloat; e: ptr LapackComplexFloat; df: ptr cfloat;
            ef: ptr LapackComplexFloat; b: ptr LapackComplexFloat; ldb: LapackInt;
            x: ptr LapackComplexFloat; ldx: LapackInt; rcond: ptr cfloat;
            ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.importc: "LAPACKE_cptsvx",
    dynlib: libName.}
proc cptsvxWork*(matrixLayout: cint; fact: char; n: LapackInt; nrhs: LapackInt;
                d: ptr cfloat; e: ptr LapackComplexFloat; df: ptr cfloat;
                ef: ptr LapackComplexFloat; b: ptr LapackComplexFloat; ldb: LapackInt;
                x: ptr LapackComplexFloat; ldx: LapackInt; rcond: ptr cfloat;
                ferr: ptr cfloat; berr: ptr cfloat; work: ptr LapackComplexFloat;
                rwork: ptr cfloat): LapackInt {.importc: "LAPACKE_cptsvx_work",
    dynlib: libName.}
proc cpttrf*(n: LapackInt; d: ptr cfloat; e: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cpttrf", dynlib: libName.}
proc cpttrfWork*(n: LapackInt; d: ptr cfloat; e: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cpttrf_work", dynlib: libName.}
proc cpttrs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            d: ptr cfloat; e: ptr LapackComplexFloat; b: ptr LapackComplexFloat;
            ldb: LapackInt): LapackInt {.importc: "LAPACKE_cpttrs", dynlib: libName.}
proc cpttrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                d: ptr cfloat; e: ptr LapackComplexFloat; b: ptr LapackComplexFloat;
                ldb: LapackInt): LapackInt {.importc: "LAPACKE_cpttrs_work",
    dynlib: libName.}
proc cspcon*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexFloat;
            ipiv: ptr LapackInt; anorm: cfloat; rcond: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cspcon", dynlib: libName.}
proc cspconWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexFloat; ipiv: ptr LapackInt; anorm: cfloat;
                rcond: ptr cfloat; work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cspcon_work", dynlib: libName.}
proc csprfs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr LapackComplexFloat; afp: ptr LapackComplexFloat;
            ipiv: ptr LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt;
            x: ptr LapackComplexFloat; ldx: LapackInt; ferr: ptr cfloat;
            berr: ptr cfloat): LapackInt {.importc: "LAPACKE_csprfs", dynlib: libName.}
proc csprfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                ap: ptr LapackComplexFloat; afp: ptr LapackComplexFloat;
                ipiv: ptr LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt;
                x: ptr LapackComplexFloat; ldx: LapackInt; ferr: ptr cfloat;
                berr: ptr cfloat; work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_csprfs_work", dynlib: libName.}
proc cspsv*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
           ap: ptr LapackComplexFloat; ipiv: ptr LapackInt; b: ptr LapackComplexFloat;
           ldb: LapackInt): LapackInt {.importc: "LAPACKE_cspsv", dynlib: libName.}
proc cspsvWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               ap: ptr LapackComplexFloat; ipiv: ptr LapackInt;
               b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cspsv_work", dynlib: libName.}
proc cspsvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr LapackComplexFloat; afp: ptr LapackComplexFloat;
            ipiv: ptr LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt;
            x: ptr LapackComplexFloat; ldx: LapackInt; rcond: ptr cfloat;
            ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.importc: "LAPACKE_cspsvx",
    dynlib: libName.}
proc cspsvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; ap: ptr LapackComplexFloat;
                afp: ptr LapackComplexFloat; ipiv: ptr LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
                ldx: LapackInt; rcond: ptr cfloat; ferr: ptr cfloat; berr: ptr cfloat;
                work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cspsvx_work", dynlib: libName.}
proc csptrf*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexFloat;
            ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_csptrf",
    dynlib: libName.}
proc csptrfWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexFloat; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_csptrf_work", dynlib: libName.}
proc csptri*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexFloat;
            ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_csptri",
    dynlib: libName.}
proc csptriWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexFloat; ipiv: ptr LapackInt;
                work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_csptri_work", dynlib: libName.}
proc csptrs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr LapackComplexFloat; ipiv: ptr LapackInt;
            b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_csptrs", dynlib: libName.}
proc csptrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                ap: ptr LapackComplexFloat; ipiv: ptr LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_csptrs_work", dynlib: libName.}
proc cstedc*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cfloat; e: ptr cfloat;
            z: ptr LapackComplexFloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_cstedc", dynlib: libName.}
proc cstedcWork*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cfloat;
                e: ptr cfloat; z: ptr LapackComplexFloat; ldz: LapackInt;
                work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat;
                lrwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cstedc_work", dynlib: libName.}
proc cstegr*(matrixLayout: cint; jobz: char; range: char; n: LapackInt; d: ptr cfloat;
            e: ptr cfloat; vl: cfloat; vu: cfloat; il: LapackInt; iu: LapackInt;
            abstol: cfloat; m: ptr LapackInt; w: ptr cfloat; z: ptr LapackComplexFloat;
            ldz: LapackInt; isuppz: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_cstegr", dynlib: libName.}
proc cstegrWork*(matrixLayout: cint; jobz: char; range: char; n: LapackInt;
                d: ptr cfloat; e: ptr cfloat; vl: cfloat; vu: cfloat; il: LapackInt;
                iu: LapackInt; abstol: cfloat; m: ptr LapackInt; w: ptr cfloat;
                z: ptr LapackComplexFloat; ldz: LapackInt; isuppz: ptr LapackInt;
                work: ptr cfloat; lwork: LapackInt; iwork: ptr LapackInt;
                liwork: LapackInt): LapackInt {.importc: "LAPACKE_cstegr_work",
    dynlib: libName.}
proc cstein*(matrixLayout: cint; n: LapackInt; d: ptr cfloat; e: ptr cfloat; m: LapackInt;
            w: ptr cfloat; iblock: ptr LapackInt; isplit: ptr LapackInt;
            z: ptr LapackComplexFloat; ldz: LapackInt; ifailv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_cstein", dynlib: libName.}
proc csteinWork*(matrixLayout: cint; n: LapackInt; d: ptr cfloat; e: ptr cfloat;
                m: LapackInt; w: ptr cfloat; iblock: ptr LapackInt;
                isplit: ptr LapackInt; z: ptr LapackComplexFloat; ldz: LapackInt;
                work: ptr cfloat; iwork: ptr LapackInt; ifailv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_cstein_work", dynlib: libName.}
proc cstemr*(matrixLayout: cint; jobz: char; range: char; n: LapackInt; d: ptr cfloat;
            e: ptr cfloat; vl: cfloat; vu: cfloat; il: LapackInt; iu: LapackInt;
            m: ptr LapackInt; w: ptr cfloat; z: ptr LapackComplexFloat; ldz: LapackInt;
            nzc: LapackInt; isuppz: ptr LapackInt; tryrac: ptr LapackLogical): LapackInt {.
    importc: "LAPACKE_cstemr", dynlib: libName.}
proc cstemrWork*(matrixLayout: cint; jobz: char; range: char; n: LapackInt;
                d: ptr cfloat; e: ptr cfloat; vl: cfloat; vu: cfloat; il: LapackInt;
                iu: LapackInt; m: ptr LapackInt; w: ptr cfloat;
                z: ptr LapackComplexFloat; ldz: LapackInt; nzc: LapackInt;
                isuppz: ptr LapackInt; tryrac: ptr LapackLogical; work: ptr cfloat;
                lwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cstemr_work", dynlib: libName.}
proc csteqr*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cfloat; e: ptr cfloat;
            z: ptr LapackComplexFloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_csteqr", dynlib: libName.}
proc csteqrWork*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cfloat;
                e: ptr cfloat; z: ptr LapackComplexFloat; ldz: LapackInt;
                work: ptr cfloat): LapackInt {.importc: "LAPACKE_csteqr_work",
    dynlib: libName.}
proc csycon*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; ipiv: ptr LapackInt; anorm: cfloat; rcond: ptr cfloat): LapackInt {.
    importc: "LAPACKE_csycon", dynlib: libName.}
proc csyconWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                anorm: cfloat; rcond: ptr cfloat; work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_csycon_work", dynlib: libName.}
proc csyconv*(matrixLayout: cint; uplo: char; way: char; n: LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
             e: ptr LapackComplexFloat): LapackInt {.importc: "LAPACKE_csyconv",
    dynlib: libName.}
proc csyconvWork*(matrixLayout: cint; uplo: char; way: char; n: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                 e: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_csyconv_work", dynlib: libName.}
proc csyequb*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
             lda: LapackInt; s: ptr cfloat; scond: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_csyequb", dynlib: libName.}
proc csyequbWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt; s: ptr cfloat;
                 scond: ptr cfloat; amax: ptr cfloat; work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_csyequb_work", dynlib: libName.}
proc csyr*(matrixLayout: cint; uplo: char; n: LapackInt; alpha: LapackComplexFloat;
          x: ptr LapackComplexFloat; incx: LapackInt; a: ptr LapackComplexFloat;
          lda: LapackInt): LapackInt {.importc: "LAPACKE_csyr", dynlib: libName.}
proc csyrWork*(matrixLayout: cint; uplo: char; n: LapackInt; alpha: LapackComplexFloat;
              x: ptr LapackComplexFloat; incx: LapackInt; a: ptr LapackComplexFloat;
              lda: LapackInt): LapackInt {.importc: "LAPACKE_csyr_work",
                                        dynlib: libName.}
proc csyrfs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; af: ptr LapackComplexFloat;
            ldaf: LapackInt; ipiv: ptr LapackInt; b: ptr LapackComplexFloat;
            ldb: LapackInt; x: ptr LapackComplexFloat; ldx: LapackInt;
            ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.importc: "LAPACKE_csyrfs",
    dynlib: libName.}
proc csyrfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt;
                af: ptr LapackComplexFloat; ldaf: LapackInt; ipiv: ptr LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
                ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat;
                work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_csyrfs_work", dynlib: libName.}
proc csyrfsx*(matrixLayout: cint; uplo: char; equed: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; af: ptr LapackComplexFloat;
             ldaf: LapackInt; ipiv: ptr LapackInt; s: ptr cfloat;
             b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
             ldx: LapackInt; rcond: ptr cfloat; berr: ptr cfloat; nErrBnds: LapackInt;
             errBndsNorm: ptr cfloat; errBndsComp: ptr cfloat; nparams: LapackInt;
             params: ptr cfloat): LapackInt {.importc: "LAPACKE_csyrfsx",
    dynlib: libName.}
proc csyrfsxWork*(matrixLayout: cint; uplo: char; equed: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                 af: ptr LapackComplexFloat; ldaf: LapackInt; ipiv: ptr LapackInt;
                 s: ptr cfloat; b: ptr LapackComplexFloat; ldb: LapackInt;
                 x: ptr LapackComplexFloat; ldx: LapackInt; rcond: ptr cfloat;
                 berr: ptr cfloat; nErrBnds: LapackInt; errBndsNorm: ptr cfloat;
                 errBndsComp: ptr cfloat; nparams: LapackInt; params: ptr cfloat;
                 work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_csyrfsx_work", dynlib: libName.}
proc csysv*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
           a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
           b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_csysv", dynlib: libName.}
proc csysvRook*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
               b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_csysv_rook", dynlib: libName.}
proc csysvRookWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                   a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                   b: ptr LapackComplexFloat; ldb: LapackInt;
                   work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_csysv_rook_work", dynlib: libName.}
proc csysvWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
               b: ptr LapackComplexFloat; ldb: LapackInt;
               work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_csysv_work", dynlib: libName.}
proc csysvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; af: ptr LapackComplexFloat;
            ldaf: LapackInt; ipiv: ptr LapackInt; b: ptr LapackComplexFloat;
            ldb: LapackInt; x: ptr LapackComplexFloat; ldx: LapackInt;
            rcond: ptr cfloat; ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.
    importc: "LAPACKE_csysvx", dynlib: libName.}
proc csysvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                af: ptr LapackComplexFloat; ldaf: LapackInt; ipiv: ptr LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
                ldx: LapackInt; rcond: ptr cfloat; ferr: ptr cfloat; berr: ptr cfloat;
                work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_csysvx_work", dynlib: libName.}
proc csysvxx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; af: ptr LapackComplexFloat;
             ldaf: LapackInt; ipiv: ptr LapackInt; equed: cstring; s: ptr cfloat;
             b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
             ldx: LapackInt; rcond: ptr cfloat; rpvgrw: ptr cfloat; berr: ptr cfloat;
             nErrBnds: LapackInt; errBndsNorm: ptr cfloat; errBndsComp: ptr cfloat;
             nparams: LapackInt; params: ptr cfloat): LapackInt {.
    importc: "LAPACKE_csysvxx", dynlib: libName.}
proc csysvxxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                 af: ptr LapackComplexFloat; ldaf: LapackInt; ipiv: ptr LapackInt;
                 equed: cstring; s: ptr cfloat; b: ptr LapackComplexFloat;
                 ldb: LapackInt; x: ptr LapackComplexFloat; ldx: LapackInt;
                 rcond: ptr cfloat; rpvgrw: ptr cfloat; berr: ptr cfloat;
                 nErrBnds: LapackInt; errBndsNorm: ptr cfloat;
                 errBndsComp: ptr cfloat; nparams: LapackInt; params: ptr cfloat;
                 work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_csysvxx_work", dynlib: libName.}
proc csyswapr*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
              lda: LapackInt; i1: LapackInt; i2: LapackInt): LapackInt {.
    importc: "LAPACKE_csyswapr", dynlib: libName.}
proc csyswaprWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                  a: ptr LapackComplexFloat; lda: LapackInt; i1: LapackInt;
                  i2: LapackInt): LapackInt {.importc: "LAPACKE_csyswapr_work",
    dynlib: libName.}
proc csytrf*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_csytrf", dynlib: libName.}
proc csytrfRook*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_csytrf_rook", dynlib: libName.}
proc csytrfRookWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                    a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                    work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_csytrf_rook_work", dynlib: libName.}
proc csytrfWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_csytrf_work", dynlib: libName.}
proc csytri*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_csytri", dynlib: libName.}
proc csytri2*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
             lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_csytri2", dynlib: libName.}
proc csytri2Work*(matrixLayout: cint; uplo: char; n: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                 work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_csytri2_work", dynlib: libName.}
proc csytri2x*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
              lda: LapackInt; ipiv: ptr LapackInt; nb: LapackInt): LapackInt {.
    importc: "LAPACKE_csytri2x", dynlib: libName.}
proc csytri2xWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                  a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                  work: ptr LapackComplexFloat; nb: LapackInt): LapackInt {.
    importc: "LAPACKE_csytri2x_work", dynlib: libName.}
proc csytriWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_csytri_work", dynlib: libName.}
proc csytrs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
            b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_csytrs", dynlib: libName.}
proc csytrs2*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
             b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_csytrs2", dynlib: libName.}
proc csytrs2Work*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                 b: ptr LapackComplexFloat; ldb: LapackInt;
                 work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_csytrs2_work", dynlib: libName.}
proc csytrsRook*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_csytrs_rook", dynlib: libName.}
proc csytrsRookWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                    a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                    b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_csytrs_rook_work", dynlib: libName.}
proc csytrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_csytrs_work", dynlib: libName.}
proc ctbcon*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
            kd: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
            rcond: ptr cfloat): LapackInt {.importc: "LAPACKE_ctbcon", dynlib: libName.}
proc ctbconWork*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
                kd: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
                rcond: ptr cfloat; work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ctbcon_work", dynlib: libName.}
proc ctbrfs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            kd: LapackInt; nrhs: LapackInt; ab: ptr LapackComplexFloat;
            ldab: LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt;
            x: ptr LapackComplexFloat; ldx: LapackInt; ferr: ptr cfloat;
            berr: ptr cfloat): LapackInt {.importc: "LAPACKE_ctbrfs", dynlib: libName.}
proc ctbrfsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                kd: LapackInt; nrhs: LapackInt; ab: ptr LapackComplexFloat;
                ldab: LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt;
                x: ptr LapackComplexFloat; ldx: LapackInt; ferr: ptr cfloat;
                berr: ptr cfloat; work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ctbrfs_work", dynlib: libName.}
proc ctbtrs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            kd: LapackInt; nrhs: LapackInt; ab: ptr LapackComplexFloat;
            ldab: LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ctbtrs", dynlib: libName.}
proc ctbtrsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                kd: LapackInt; nrhs: LapackInt; ab: ptr LapackComplexFloat;
                ldab: LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ctbtrs_work", dynlib: libName.}
proc ctfsm*(matrixLayout: cint; transr: char; side: char; uplo: char; trans: char;
           diag: char; m: LapackInt; n: LapackInt; alpha: LapackComplexFloat;
           a: ptr LapackComplexFloat; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ctfsm", dynlib: libName.}
proc ctfsmWork*(matrixLayout: cint; transr: char; side: char; uplo: char; trans: char;
               diag: char; m: LapackInt; n: LapackInt; alpha: LapackComplexFloat;
               a: ptr LapackComplexFloat; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ctfsm_work", dynlib: libName.}
proc ctftri*(matrixLayout: cint; transr: char; uplo: char; diag: char; n: LapackInt;
            a: ptr LapackComplexFloat): LapackInt {.importc: "LAPACKE_ctftri",
    dynlib: libName.}
proc ctftriWork*(matrixLayout: cint; transr: char; uplo: char; diag: char; n: LapackInt;
                a: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_ctftri_work", dynlib: libName.}
proc ctfttp*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
            arf: ptr LapackComplexFloat; ap: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_ctfttp", dynlib: libName.}
proc ctfttpWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                arf: ptr LapackComplexFloat; ap: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_ctfttp_work", dynlib: libName.}
proc ctfttr*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
            arf: ptr LapackComplexFloat; a: ptr LapackComplexFloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_ctfttr", dynlib: libName.}
proc ctfttrWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                arf: ptr LapackComplexFloat; a: ptr LapackComplexFloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_ctfttr_work", dynlib: libName.}
proc ctgevc*(matrixLayout: cint; side: char; howmny: char; select: ptr LapackLogical;
            n: LapackInt; s: ptr LapackComplexFloat; lds: LapackInt;
            p: ptr LapackComplexFloat; ldp: LapackInt; vl: ptr LapackComplexFloat;
            ldvl: LapackInt; vr: ptr LapackComplexFloat; ldvr: LapackInt;
            mm: LapackInt; m: ptr LapackInt): LapackInt {.importc: "LAPACKE_ctgevc",
    dynlib: libName.}
proc ctgevcWork*(matrixLayout: cint; side: char; howmny: char;
                select: ptr LapackLogical; n: LapackInt; s: ptr LapackComplexFloat;
                lds: LapackInt; p: ptr LapackComplexFloat; ldp: LapackInt;
                vl: ptr LapackComplexFloat; ldvl: LapackInt;
                vr: ptr LapackComplexFloat; ldvr: LapackInt; mm: LapackInt;
                m: ptr LapackInt; work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ctgevc_work", dynlib: libName.}
proc ctgexc*(matrixLayout: cint; wantq: LapackLogical; wantz: LapackLogical;
            n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            b: ptr LapackComplexFloat; ldb: LapackInt; q: ptr LapackComplexFloat;
            ldq: LapackInt; z: ptr LapackComplexFloat; ldz: LapackInt; ifst: LapackInt;
            ilst: LapackInt): LapackInt {.importc: "LAPACKE_ctgexc", dynlib: libName.}
proc ctgexcWork*(matrixLayout: cint; wantq: LapackLogical; wantz: LapackLogical;
                n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt; q: ptr LapackComplexFloat;
                ldq: LapackInt; z: ptr LapackComplexFloat; ldz: LapackInt;
                ifst: LapackInt; ilst: LapackInt): LapackInt {.
    importc: "LAPACKE_ctgexc_work", dynlib: libName.}
proc ctgsen*(matrixLayout: cint; ijob: LapackInt; wantq: LapackLogical;
            wantz: LapackLogical; select: ptr LapackLogical; n: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
            ldb: LapackInt; alpha: ptr LapackComplexFloat;
            beta: ptr LapackComplexFloat; q: ptr LapackComplexFloat; ldq: LapackInt;
            z: ptr LapackComplexFloat; ldz: LapackInt; m: ptr LapackInt; pl: ptr cfloat;
            pr: ptr cfloat; dif: ptr cfloat): LapackInt {.importc: "LAPACKE_ctgsen",
    dynlib: libName.}
proc ctgsenWork*(matrixLayout: cint; ijob: LapackInt; wantq: LapackLogical;
                wantz: LapackLogical; select: ptr LapackLogical; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
                ldb: LapackInt; alpha: ptr LapackComplexFloat;
                beta: ptr LapackComplexFloat; q: ptr LapackComplexFloat;
                ldq: LapackInt; z: ptr LapackComplexFloat; ldz: LapackInt;
                m: ptr LapackInt; pl: ptr cfloat; pr: ptr cfloat; dif: ptr cfloat;
                work: ptr LapackComplexFloat; lwork: LapackInt; iwork: ptr LapackInt;
                liwork: LapackInt): LapackInt {.importc: "LAPACKE_ctgsen_work",
    dynlib: libName.}
proc ctgsja*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
            p: LapackInt; n: LapackInt; k: LapackInt; l: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
            ldb: LapackInt; tola: cfloat; tolb: cfloat; alpha: ptr cfloat;
            beta: ptr cfloat; u: ptr LapackComplexFloat; ldu: LapackInt;
            v: ptr LapackComplexFloat; ldv: LapackInt; q: ptr LapackComplexFloat;
            ldq: LapackInt; ncycle: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ctgsja", dynlib: libName.}
proc ctgsjaWork*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
                p: LapackInt; n: LapackInt; k: LapackInt; l: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
                ldb: LapackInt; tola: cfloat; tolb: cfloat; alpha: ptr cfloat;
                beta: ptr cfloat; u: ptr LapackComplexFloat; ldu: LapackInt;
                v: ptr LapackComplexFloat; ldv: LapackInt; q: ptr LapackComplexFloat;
                ldq: LapackInt; work: ptr LapackComplexFloat; ncycle: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ctgsja_work", dynlib: libName.}
proc ctgsna*(matrixLayout: cint; job: char; howmny: char; select: ptr LapackLogical;
            n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            b: ptr LapackComplexFloat; ldb: LapackInt; vl: ptr LapackComplexFloat;
            ldvl: LapackInt; vr: ptr LapackComplexFloat; ldvr: LapackInt;
            s: ptr cfloat; dif: ptr cfloat; mm: LapackInt; m: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ctgsna", dynlib: libName.}
proc ctgsnaWork*(matrixLayout: cint; job: char; howmny: char;
                select: ptr LapackLogical; n: LapackInt; a: ptr LapackComplexFloat;
                lda: LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt;
                vl: ptr LapackComplexFloat; ldvl: LapackInt;
                vr: ptr LapackComplexFloat; ldvr: LapackInt; s: ptr cfloat;
                dif: ptr cfloat; mm: LapackInt; m: ptr LapackInt;
                work: ptr LapackComplexFloat; lwork: LapackInt; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ctgsna_work", dynlib: libName.}
proc ctgsyl*(matrixLayout: cint; trans: char; ijob: LapackInt; m: LapackInt;
            n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            b: ptr LapackComplexFloat; ldb: LapackInt; c: ptr LapackComplexFloat;
            ldc: LapackInt; d: ptr LapackComplexFloat; ldd: LapackInt;
            e: ptr LapackComplexFloat; lde: LapackInt; f: ptr LapackComplexFloat;
            ldf: LapackInt; scale: ptr cfloat; dif: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ctgsyl", dynlib: libName.}
proc ctgsylWork*(matrixLayout: cint; trans: char; ijob: LapackInt; m: LapackInt;
                n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt; c: ptr LapackComplexFloat;
                ldc: LapackInt; d: ptr LapackComplexFloat; ldd: LapackInt;
                e: ptr LapackComplexFloat; lde: LapackInt; f: ptr LapackComplexFloat;
                ldf: LapackInt; scale: ptr cfloat; dif: ptr cfloat;
                work: ptr LapackComplexFloat; lwork: LapackInt; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ctgsyl_work", dynlib: libName.}
proc ctpcon*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
            ap: ptr LapackComplexFloat; rcond: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ctpcon", dynlib: libName.}
proc ctpconWork*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
                ap: ptr LapackComplexFloat; rcond: ptr cfloat;
                work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ctpcon_work", dynlib: libName.}
proc ctpmqrt*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
             k: LapackInt; l: LapackInt; nb: LapackInt; v: ptr LapackComplexFloat;
             ldv: LapackInt; t: ptr LapackComplexFloat; ldt: LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
             ldb: LapackInt): LapackInt {.importc: "LAPACKE_ctpmqrt", dynlib: libName.}
proc ctpmqrtWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt;
                 n: LapackInt; k: LapackInt; l: LapackInt; nb: LapackInt;
                 v: ptr LapackComplexFloat; ldv: LapackInt;
                 t: ptr LapackComplexFloat; ldt: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt;
                 b: ptr LapackComplexFloat; ldb: LapackInt;
                 work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_ctpmqrt_work", dynlib: libName.}
proc ctpqrt*(matrixLayout: cint; m: LapackInt; n: LapackInt; l: LapackInt; nb: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
            ldb: LapackInt; t: ptr LapackComplexFloat; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_ctpqrt", dynlib: libName.}
proc ctpqrt2*(matrixLayout: cint; m: LapackInt; n: LapackInt; l: LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; b: ptr LapackComplexFloat;
             ldb: LapackInt; t: ptr LapackComplexFloat; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_ctpqrt2", dynlib: libName.}
proc ctpqrt2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt; l: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt;
                 b: ptr LapackComplexFloat; ldb: LapackInt;
                 t: ptr LapackComplexFloat; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_ctpqrt2_work", dynlib: libName.}
proc ctpqrtWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; l: LapackInt;
                nb: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt; t: ptr LapackComplexFloat;
                ldt: LapackInt; work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_ctpqrt_work", dynlib: libName.}
proc ctprfb*(matrixLayout: cint; side: char; trans: char; direct: char; storev: char;
            m: LapackInt; n: LapackInt; k: LapackInt; l: LapackInt;
            v: ptr LapackComplexFloat; ldv: LapackInt; t: ptr LapackComplexFloat;
            ldt: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ctprfb", dynlib: libName.}
proc ctprfbWork*(matrixLayout: cint; side: char; trans: char; direct: char; storev: char;
                m: LapackInt; n: LapackInt; k: LapackInt; l: LapackInt;
                v: ptr LapackComplexFloat; ldv: LapackInt; t: ptr LapackComplexFloat;
                ldt: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt;
                work: ptr LapackComplexFloat; ldwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ctprfb_work", dynlib: libName.}
proc ctprfs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            nrhs: LapackInt; ap: ptr LapackComplexFloat; b: ptr LapackComplexFloat;
            ldb: LapackInt; x: ptr LapackComplexFloat; ldx: LapackInt;
            ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.importc: "LAPACKE_ctprfs",
    dynlib: libName.}
proc ctprfsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                nrhs: LapackInt; ap: ptr LapackComplexFloat;
                b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
                ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat;
                work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ctprfs_work", dynlib: libName.}
proc ctptri*(matrixLayout: cint; uplo: char; diag: char; n: LapackInt;
            ap: ptr LapackComplexFloat): LapackInt {.importc: "LAPACKE_ctptri",
    dynlib: libName.}
proc ctptriWork*(matrixLayout: cint; uplo: char; diag: char; n: LapackInt;
                ap: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_ctptri_work", dynlib: libName.}
proc ctptrs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            nrhs: LapackInt; ap: ptr LapackComplexFloat; b: ptr LapackComplexFloat;
            ldb: LapackInt): LapackInt {.importc: "LAPACKE_ctptrs", dynlib: libName.}
proc ctptrsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                nrhs: LapackInt; ap: ptr LapackComplexFloat;
                b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ctptrs_work", dynlib: libName.}
proc ctpttf*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
            ap: ptr LapackComplexFloat; arf: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_ctpttf", dynlib: libName.}
proc ctpttfWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                ap: ptr LapackComplexFloat; arf: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_ctpttf_work", dynlib: libName.}
proc ctpttr*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexFloat;
            a: ptr LapackComplexFloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_ctpttr", dynlib: libName.}
proc ctpttrWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexFloat; a: ptr LapackComplexFloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_ctpttr_work", dynlib: libName.}
proc ctrcon*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; rcond: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ctrcon", dynlib: libName.}
proc ctrconWork*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; rcond: ptr cfloat;
                work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ctrcon_work", dynlib: libName.}
proc ctrevc*(matrixLayout: cint; side: char; howmny: char; select: ptr LapackLogical;
            n: LapackInt; t: ptr LapackComplexFloat; ldt: LapackInt;
            vl: ptr LapackComplexFloat; ldvl: LapackInt; vr: ptr LapackComplexFloat;
            ldvr: LapackInt; mm: LapackInt; m: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ctrevc", dynlib: libName.}
proc ctrevcWork*(matrixLayout: cint; side: char; howmny: char;
                select: ptr LapackLogical; n: LapackInt; t: ptr LapackComplexFloat;
                ldt: LapackInt; vl: ptr LapackComplexFloat; ldvl: LapackInt;
                vr: ptr LapackComplexFloat; ldvr: LapackInt; mm: LapackInt;
                m: ptr LapackInt; work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ctrevc_work", dynlib: libName.}
proc ctrexc*(matrixLayout: cint; compq: char; n: LapackInt; t: ptr LapackComplexFloat;
            ldt: LapackInt; q: ptr LapackComplexFloat; ldq: LapackInt; ifst: LapackInt;
            ilst: LapackInt): LapackInt {.importc: "LAPACKE_ctrexc", dynlib: libName.}
proc ctrexcWork*(matrixLayout: cint; compq: char; n: LapackInt;
                t: ptr LapackComplexFloat; ldt: LapackInt; q: ptr LapackComplexFloat;
                ldq: LapackInt; ifst: LapackInt; ilst: LapackInt): LapackInt {.
    importc: "LAPACKE_ctrexc_work", dynlib: libName.}
proc ctrrfs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            nrhs: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
            ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ctrrfs", dynlib: libName.}
proc ctrrfsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                nrhs: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt; x: ptr LapackComplexFloat;
                ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat;
                work: ptr LapackComplexFloat; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ctrrfs_work", dynlib: libName.}
proc ctrsen*(matrixLayout: cint; job: char; compq: char; select: ptr LapackLogical;
            n: LapackInt; t: ptr LapackComplexFloat; ldt: LapackInt;
            q: ptr LapackComplexFloat; ldq: LapackInt; w: ptr LapackComplexFloat;
            m: ptr LapackInt; s: ptr cfloat; sep: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ctrsen", dynlib: libName.}
proc ctrsenWork*(matrixLayout: cint; job: char; compq: char; select: ptr LapackLogical;
                n: LapackInt; t: ptr LapackComplexFloat; ldt: LapackInt;
                q: ptr LapackComplexFloat; ldq: LapackInt; w: ptr LapackComplexFloat;
                m: ptr LapackInt; s: ptr cfloat; sep: ptr cfloat;
                work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ctrsen_work", dynlib: libName.}
proc ctrsna*(matrixLayout: cint; job: char; howmny: char; select: ptr LapackLogical;
            n: LapackInt; t: ptr LapackComplexFloat; ldt: LapackInt;
            vl: ptr LapackComplexFloat; ldvl: LapackInt; vr: ptr LapackComplexFloat;
            ldvr: LapackInt; s: ptr cfloat; sep: ptr cfloat; mm: LapackInt;
            m: ptr LapackInt): LapackInt {.importc: "LAPACKE_ctrsna", dynlib: libName.}
proc ctrsnaWork*(matrixLayout: cint; job: char; howmny: char;
                select: ptr LapackLogical; n: LapackInt; t: ptr LapackComplexFloat;
                ldt: LapackInt; vl: ptr LapackComplexFloat; ldvl: LapackInt;
                vr: ptr LapackComplexFloat; ldvr: LapackInt; s: ptr cfloat;
                sep: ptr cfloat; mm: LapackInt; m: ptr LapackInt;
                work: ptr LapackComplexFloat; ldwork: LapackInt; rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ctrsna_work", dynlib: libName.}
proc ctrsyl*(matrixLayout: cint; trana: char; tranb: char; isgn: LapackInt; m: LapackInt;
            n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            b: ptr LapackComplexFloat; ldb: LapackInt; c: ptr LapackComplexFloat;
            ldc: LapackInt; scale: ptr cfloat): LapackInt {.importc: "LAPACKE_ctrsyl",
    dynlib: libName.}
proc ctrsylWork*(matrixLayout: cint; trana: char; tranb: char; isgn: LapackInt;
                m: LapackInt; n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt; c: ptr LapackComplexFloat;
                ldc: LapackInt; scale: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ctrsyl_work", dynlib: libName.}
proc ctrtri*(matrixLayout: cint; uplo: char; diag: char; n: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_ctrtri", dynlib: libName.}
proc ctrtriWork*(matrixLayout: cint; uplo: char; diag: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_ctrtri_work", dynlib: libName.}
proc ctrtrs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            nrhs: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ctrtrs", dynlib: libName.}
proc ctrtrsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                nrhs: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ctrtrs_work", dynlib: libName.}
proc ctrttf*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; arf: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_ctrttf", dynlib: libName.}
proc ctrttfWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt;
                arf: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_ctrttf_work", dynlib: libName.}
proc ctrttp*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; ap: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_ctrttp", dynlib: libName.}
proc ctrttpWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt; ap: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_ctrttp_work", dynlib: libName.}
proc ctzrzf*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; tau: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_ctzrzf", dynlib: libName.}
proc ctzrzfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt;
                tau: ptr LapackComplexFloat; work: ptr LapackComplexFloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_ctzrzf_work",
    dynlib: libName.}
proc cunbdb*(matrixLayout: cint; trans: char; signs: char; m: LapackInt; p: LapackInt;
            q: LapackInt; x11: ptr LapackComplexFloat; ldx11: LapackInt;
            x12: ptr LapackComplexFloat; ldx12: LapackInt;
            x21: ptr LapackComplexFloat; ldx21: LapackInt;
            x22: ptr LapackComplexFloat; ldx22: LapackInt; theta: ptr cfloat;
            phi: ptr cfloat; taup1: ptr LapackComplexFloat;
            taup2: ptr LapackComplexFloat; tauq1: ptr LapackComplexFloat;
            tauq2: ptr LapackComplexFloat): LapackInt {.importc: "LAPACKE_cunbdb",
    dynlib: libName.}
proc cunbdbWork*(matrixLayout: cint; trans: char; signs: char; m: LapackInt;
                p: LapackInt; q: LapackInt; x11: ptr LapackComplexFloat;
                ldx11: LapackInt; x12: ptr LapackComplexFloat; ldx12: LapackInt;
                x21: ptr LapackComplexFloat; ldx21: LapackInt;
                x22: ptr LapackComplexFloat; ldx22: LapackInt; theta: ptr cfloat;
                phi: ptr cfloat; taup1: ptr LapackComplexFloat;
                taup2: ptr LapackComplexFloat; tauq1: ptr LapackComplexFloat;
                tauq2: ptr LapackComplexFloat; work: ptr LapackComplexFloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_cunbdb_work",
    dynlib: libName.}
proc cuncsd*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char; jobv2t: char;
            trans: char; signs: char; m: LapackInt; p: LapackInt; q: LapackInt;
            x11: ptr LapackComplexFloat; ldx11: LapackInt;
            x12: ptr LapackComplexFloat; ldx12: LapackInt;
            x21: ptr LapackComplexFloat; ldx21: LapackInt;
            x22: ptr LapackComplexFloat; ldx22: LapackInt; theta: ptr cfloat;
            u1: ptr LapackComplexFloat; ldu1: LapackInt; u2: ptr LapackComplexFloat;
            ldu2: LapackInt; v1t: ptr LapackComplexFloat; ldv1t: LapackInt;
            v2t: ptr LapackComplexFloat; ldv2t: LapackInt): LapackInt {.
    importc: "LAPACKE_cuncsd", dynlib: libName.}
proc cuncsd2by1*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char;
                m: LapackInt; p: LapackInt; q: LapackInt; x11: ptr LapackComplexFloat;
                ldx11: LapackInt; x21: ptr LapackComplexFloat; ldx21: LapackInt;
                theta: ptr cfloat; u1: ptr LapackComplexFloat; ldu1: LapackInt;
                u2: ptr LapackComplexFloat; ldu2: LapackInt;
                v1t: ptr LapackComplexFloat; ldv1t: LapackInt): LapackInt {.
    importc: "LAPACKE_cuncsd2by1", dynlib: libName.}
proc cuncsd2by1Work*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char;
                    m: LapackInt; p: LapackInt; q: LapackInt;
                    x11: ptr LapackComplexFloat; ldx11: LapackInt;
                    x21: ptr LapackComplexFloat; ldx21: LapackInt; theta: ptr cfloat;
                    u1: ptr LapackComplexFloat; ldu1: LapackInt;
                    u2: ptr LapackComplexFloat; ldu2: LapackInt;
                    v1t: ptr LapackComplexFloat; ldv1t: LapackInt;
                    work: ptr LapackComplexFloat; lwork: LapackInt;
                    rwork: ptr cfloat; lrwork: LapackInt; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_cuncsd2by1_work", dynlib: libName.}
proc cuncsdWork*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char;
                jobv2t: char; trans: char; signs: char; m: LapackInt; p: LapackInt;
                q: LapackInt; x11: ptr LapackComplexFloat; ldx11: LapackInt;
                x12: ptr LapackComplexFloat; ldx12: LapackInt;
                x21: ptr LapackComplexFloat; ldx21: LapackInt;
                x22: ptr LapackComplexFloat; ldx22: LapackInt; theta: ptr cfloat;
                u1: ptr LapackComplexFloat; ldu1: LapackInt;
                u2: ptr LapackComplexFloat; ldu2: LapackInt;
                v1t: ptr LapackComplexFloat; ldv1t: LapackInt;
                v2t: ptr LapackComplexFloat; ldv2t: LapackInt;
                work: ptr LapackComplexFloat; lwork: LapackInt; rwork: ptr cfloat;
                lrwork: LapackInt; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_cuncsd_work", dynlib: libName.}
proc cungbr*(matrixLayout: cint; vect: char; m: LapackInt; n: LapackInt; k: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; tau: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cungbr", dynlib: libName.}
proc cungbrWork*(matrixLayout: cint; vect: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                tau: ptr LapackComplexFloat; work: ptr LapackComplexFloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_cungbr_work",
    dynlib: libName.}
proc cunghr*(matrixLayout: cint; n: LapackInt; ilo: LapackInt; ihi: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; tau: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cunghr", dynlib: libName.}
proc cunghrWork*(matrixLayout: cint; n: LapackInt; ilo: LapackInt; ihi: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt;
                tau: ptr LapackComplexFloat; work: ptr LapackComplexFloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_cunghr_work",
    dynlib: libName.}
proc cunglq*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; tau: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cunglq", dynlib: libName.}
proc cunglqWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt;
                tau: ptr LapackComplexFloat; work: ptr LapackComplexFloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_cunglq_work",
    dynlib: libName.}
proc cungql*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; tau: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cungql", dynlib: libName.}
proc cungqlWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt;
                tau: ptr LapackComplexFloat; work: ptr LapackComplexFloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_cungql_work",
    dynlib: libName.}
proc cungqr*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; tau: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cungqr", dynlib: libName.}
proc cungqrWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt;
                tau: ptr LapackComplexFloat; work: ptr LapackComplexFloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_cungqr_work",
    dynlib: libName.}
proc cungrq*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
            a: ptr LapackComplexFloat; lda: LapackInt; tau: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cungrq", dynlib: libName.}
proc cungrqWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt;
                tau: ptr LapackComplexFloat; work: ptr LapackComplexFloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_cungrq_work",
    dynlib: libName.}
proc cungtr*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
            lda: LapackInt; tau: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cungtr", dynlib: libName.}
proc cungtrWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexFloat; lda: LapackInt;
                tau: ptr LapackComplexFloat; work: ptr LapackComplexFloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_cungtr_work",
    dynlib: libName.}
proc cunmbr*(matrixLayout: cint; vect: char; side: char; trans: char; m: LapackInt;
            n: LapackInt; k: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            tau: ptr LapackComplexFloat; c: ptr LapackComplexFloat; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_cunmbr", dynlib: libName.}
proc cunmbrWork*(matrixLayout: cint; vect: char; side: char; trans: char; m: LapackInt;
                n: LapackInt; k: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                tau: ptr LapackComplexFloat; c: ptr LapackComplexFloat;
                ldc: LapackInt; work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cunmbr_work", dynlib: libName.}
proc cunmhr*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            ilo: LapackInt; ihi: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            tau: ptr LapackComplexFloat; c: ptr LapackComplexFloat; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_cunmhr", dynlib: libName.}
proc cunmhrWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                ilo: LapackInt; ihi: LapackInt; a: ptr LapackComplexFloat;
                lda: LapackInt; tau: ptr LapackComplexFloat;
                c: ptr LapackComplexFloat; ldc: LapackInt;
                work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cunmhr_work", dynlib: libName.}
proc cunmlq*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            tau: ptr LapackComplexFloat; c: ptr LapackComplexFloat; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_cunmlq", dynlib: libName.}
proc cunmlqWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                tau: ptr LapackComplexFloat; c: ptr LapackComplexFloat;
                ldc: LapackInt; work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cunmlq_work", dynlib: libName.}
proc cunmql*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            tau: ptr LapackComplexFloat; c: ptr LapackComplexFloat; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_cunmql", dynlib: libName.}
proc cunmqlWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                tau: ptr LapackComplexFloat; c: ptr LapackComplexFloat;
                ldc: LapackInt; work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cunmql_work", dynlib: libName.}
proc cunmqr*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            tau: ptr LapackComplexFloat; c: ptr LapackComplexFloat; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_cunmqr", dynlib: libName.}
proc cunmqrWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                tau: ptr LapackComplexFloat; c: ptr LapackComplexFloat;
                ldc: LapackInt; work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cunmqr_work", dynlib: libName.}
proc cunmrq*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            tau: ptr LapackComplexFloat; c: ptr LapackComplexFloat; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_cunmrq", dynlib: libName.}
proc cunmrqWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                tau: ptr LapackComplexFloat; c: ptr LapackComplexFloat;
                ldc: LapackInt; work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cunmrq_work", dynlib: libName.}
proc cunmrz*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; l: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            tau: ptr LapackComplexFloat; c: ptr LapackComplexFloat; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_cunmrz", dynlib: libName.}
proc cunmrzWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; l: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                tau: ptr LapackComplexFloat; c: ptr LapackComplexFloat;
                ldc: LapackInt; work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cunmrz_work", dynlib: libName.}
proc cunmtr*(matrixLayout: cint; side: char; uplo: char; trans: char; m: LapackInt;
            n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            tau: ptr LapackComplexFloat; c: ptr LapackComplexFloat; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_cunmtr", dynlib: libName.}
proc cunmtrWork*(matrixLayout: cint; side: char; uplo: char; trans: char; m: LapackInt;
                n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                tau: ptr LapackComplexFloat; c: ptr LapackComplexFloat;
                ldc: LapackInt; work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cunmtr_work", dynlib: libName.}
proc cupgtr*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexFloat;
            tau: ptr LapackComplexFloat; q: ptr LapackComplexFloat; ldq: LapackInt): LapackInt {.
    importc: "LAPACKE_cupgtr", dynlib: libName.}
proc cupgtrWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexFloat; tau: ptr LapackComplexFloat;
                q: ptr LapackComplexFloat; ldq: LapackInt;
                work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cupgtr_work", dynlib: libName.}
proc cupmtr*(matrixLayout: cint; side: char; uplo: char; trans: char; m: LapackInt;
            n: LapackInt; ap: ptr LapackComplexFloat; tau: ptr LapackComplexFloat;
            c: ptr LapackComplexFloat; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_cupmtr", dynlib: libName.}
proc cupmtrWork*(matrixLayout: cint; side: char; uplo: char; trans: char; m: LapackInt;
                n: LapackInt; ap: ptr LapackComplexFloat;
                tau: ptr LapackComplexFloat; c: ptr LapackComplexFloat;
                ldc: LapackInt; work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_cupmtr_work", dynlib: libName.}
proc dbbcsd*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char; jobv2t: char;
            trans: char; m: LapackInt; p: LapackInt; q: LapackInt; theta: ptr cdouble;
            phi: ptr cdouble; u1: ptr cdouble; ldu1: LapackInt; u2: ptr cdouble;
            ldu2: LapackInt; v1t: ptr cdouble; ldv1t: LapackInt; v2t: ptr cdouble;
            ldv2t: LapackInt; b11d: ptr cdouble; b11e: ptr cdouble; b12d: ptr cdouble;
            b12e: ptr cdouble; b21d: ptr cdouble; b21e: ptr cdouble; b22d: ptr cdouble;
            b22e: ptr cdouble): LapackInt {.importc: "LAPACKE_dbbcsd", dynlib: libName.}
proc dbbcsdWork*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char;
                jobv2t: char; trans: char; m: LapackInt; p: LapackInt; q: LapackInt;
                theta: ptr cdouble; phi: ptr cdouble; u1: ptr cdouble; ldu1: LapackInt;
                u2: ptr cdouble; ldu2: LapackInt; v1t: ptr cdouble; ldv1t: LapackInt;
                v2t: ptr cdouble; ldv2t: LapackInt; b11d: ptr cdouble;
                b11e: ptr cdouble; b12d: ptr cdouble; b12e: ptr cdouble;
                b21d: ptr cdouble; b21e: ptr cdouble; b22d: ptr cdouble;
                b22e: ptr cdouble; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dbbcsd_work", dynlib: libName.}
proc dbdsdc*(matrixLayout: cint; uplo: char; compq: char; n: LapackInt; d: ptr cdouble;
            e: ptr cdouble; u: ptr cdouble; ldu: LapackInt; vt: ptr cdouble;
            ldvt: LapackInt; q: ptr cdouble; iq: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dbdsdc", dynlib: libName.}
proc dbdsdcWork*(matrixLayout: cint; uplo: char; compq: char; n: LapackInt;
                d: ptr cdouble; e: ptr cdouble; u: ptr cdouble; ldu: LapackInt;
                vt: ptr cdouble; ldvt: LapackInt; q: ptr cdouble; iq: ptr LapackInt;
                work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dbdsdc_work", dynlib: libName.}
proc dbdsqr*(matrixLayout: cint; uplo: char; n: LapackInt; ncvt: LapackInt;
            nru: LapackInt; ncc: LapackInt; d: ptr cdouble; e: ptr cdouble;
            vt: ptr cdouble; ldvt: LapackInt; u: ptr cdouble; ldu: LapackInt;
            c: ptr cdouble; ldc: LapackInt): LapackInt {.importc: "LAPACKE_dbdsqr",
    dynlib: libName.}
proc dbdsqrWork*(matrixLayout: cint; uplo: char; n: LapackInt; ncvt: LapackInt;
                nru: LapackInt; ncc: LapackInt; d: ptr cdouble; e: ptr cdouble;
                vt: ptr cdouble; ldvt: LapackInt; u: ptr cdouble; ldu: LapackInt;
                c: ptr cdouble; ldc: LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dbdsqr_work", dynlib: libName.}
proc dbdsvdx*(matrixLayout: cint; uplo: char; jobz: char; range: char; n: LapackInt;
             d: ptr cdouble; e: ptr cdouble; vl: cdouble; vu: cdouble; il: LapackInt;
             iu: LapackInt; ns: ptr LapackInt; s: ptr cdouble; z: ptr cdouble;
             ldz: LapackInt; superb: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dbdsvdx", dynlib: libName.}
proc dbdsvdxWork*(matrixLayout: cint; uplo: char; jobz: char; range: char; n: LapackInt;
                 d: ptr cdouble; e: ptr cdouble; vl: cdouble; vu: cdouble; il: LapackInt;
                 iu: LapackInt; ns: ptr LapackInt; s: ptr cdouble; z: ptr cdouble;
                 ldz: LapackInt; work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dbdsvdx_work", dynlib: libName.}
proc ddisna*(job: char; m: LapackInt; n: LapackInt; d: ptr cdouble; sep: ptr cdouble): LapackInt {.
    importc: "LAPACKE_ddisna", dynlib: libName.}
proc ddisnaWork*(job: char; m: LapackInt; n: LapackInt; d: ptr cdouble; sep: ptr cdouble): LapackInt {.
    importc: "LAPACKE_ddisna_work", dynlib: libName.}
proc dgbbrd*(matrixLayout: cint; vect: char; m: LapackInt; n: LapackInt; ncc: LapackInt;
            kl: LapackInt; ku: LapackInt; ab: ptr cdouble; ldab: LapackInt;
            d: ptr cdouble; e: ptr cdouble; q: ptr cdouble; ldq: LapackInt;
            pt: ptr cdouble; ldpt: LapackInt; c: ptr cdouble; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_dgbbrd", dynlib: libName.}
proc dgbbrdWork*(matrixLayout: cint; vect: char; m: LapackInt; n: LapackInt;
                ncc: LapackInt; kl: LapackInt; ku: LapackInt; ab: ptr cdouble;
                ldab: LapackInt; d: ptr cdouble; e: ptr cdouble; q: ptr cdouble;
                ldq: LapackInt; pt: ptr cdouble; ldpt: LapackInt; c: ptr cdouble;
                ldc: LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgbbrd_work", dynlib: libName.}
proc dgbcon*(matrixLayout: cint; norm: char; n: LapackInt; kl: LapackInt; ku: LapackInt;
            ab: ptr cdouble; ldab: LapackInt; ipiv: ptr LapackInt; anorm: cdouble;
            rcond: ptr cdouble): LapackInt {.importc: "LAPACKE_dgbcon",
    dynlib: libName.}
proc dgbconWork*(matrixLayout: cint; norm: char; n: LapackInt; kl: LapackInt;
                ku: LapackInt; ab: ptr cdouble; ldab: LapackInt; ipiv: ptr LapackInt;
                anorm: cdouble; rcond: ptr cdouble; work: ptr cdouble;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dgbcon_work",
    dynlib: libName.}
proc dgbequ*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
            ku: LapackInt; ab: ptr cdouble; ldab: LapackInt; r: ptr cdouble;
            c: ptr cdouble; rowcnd: ptr cdouble; colcnd: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgbequ", dynlib: libName.}
proc dgbequWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
                ku: LapackInt; ab: ptr cdouble; ldab: LapackInt; r: ptr cdouble;
                c: ptr cdouble; rowcnd: ptr cdouble; colcnd: ptr cdouble;
                amax: ptr cdouble): LapackInt {.importc: "LAPACKE_dgbequ_work",
    dynlib: libName.}
proc dgbequb*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
             ku: LapackInt; ab: ptr cdouble; ldab: LapackInt; r: ptr cdouble;
             c: ptr cdouble; rowcnd: ptr cdouble; colcnd: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgbequb", dynlib: libName.}
proc dgbequbWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
                 ku: LapackInt; ab: ptr cdouble; ldab: LapackInt; r: ptr cdouble;
                 c: ptr cdouble; rowcnd: ptr cdouble; colcnd: ptr cdouble;
                 amax: ptr cdouble): LapackInt {.importc: "LAPACKE_dgbequb_work",
    dynlib: libName.}
proc dgbrfs*(matrixLayout: cint; trans: char; n: LapackInt; kl: LapackInt; ku: LapackInt;
            nrhs: LapackInt; ab: ptr cdouble; ldab: LapackInt; afb: ptr cdouble;
            ldafb: LapackInt; ipiv: ptr LapackInt; b: ptr cdouble; ldb: LapackInt;
            x: ptr cdouble; ldx: LapackInt; ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgbrfs", dynlib: libName.}
proc dgbrfsWork*(matrixLayout: cint; trans: char; n: LapackInt; kl: LapackInt;
                ku: LapackInt; nrhs: LapackInt; ab: ptr cdouble; ldab: LapackInt;
                afb: ptr cdouble; ldafb: LapackInt; ipiv: ptr LapackInt; b: ptr cdouble;
                ldb: LapackInt; x: ptr cdouble; ldx: LapackInt; ferr: ptr cdouble;
                berr: ptr cdouble; work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgbrfs_work", dynlib: libName.}
proc dgbrfsx*(matrixLayout: cint; trans: char; equed: char; n: LapackInt; kl: LapackInt;
             ku: LapackInt; nrhs: LapackInt; ab: ptr cdouble; ldab: LapackInt;
             afb: ptr cdouble; ldafb: LapackInt; ipiv: ptr LapackInt; r: ptr cdouble;
             c: ptr cdouble; b: ptr cdouble; ldb: LapackInt; x: ptr cdouble;
             ldx: LapackInt; rcond: ptr cdouble; berr: ptr cdouble; nErrBnds: LapackInt;
             errBndsNorm: ptr cdouble; errBndsComp: ptr cdouble; nparams: LapackInt;
             params: ptr cdouble): LapackInt {.importc: "LAPACKE_dgbrfsx",
    dynlib: libName.}
proc dgbrfsxWork*(matrixLayout: cint; trans: char; equed: char; n: LapackInt;
                 kl: LapackInt; ku: LapackInt; nrhs: LapackInt; ab: ptr cdouble;
                 ldab: LapackInt; afb: ptr cdouble; ldafb: LapackInt;
                 ipiv: ptr LapackInt; r: ptr cdouble; c: ptr cdouble; b: ptr cdouble;
                 ldb: LapackInt; x: ptr cdouble; ldx: LapackInt; rcond: ptr cdouble;
                 berr: ptr cdouble; nErrBnds: LapackInt; errBndsNorm: ptr cdouble;
                 errBndsComp: ptr cdouble; nparams: LapackInt; params: ptr cdouble;
                 work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgbrfsx_work", dynlib: libName.}
proc dgbsv*(matrixLayout: cint; n: LapackInt; kl: LapackInt; ku: LapackInt;
           nrhs: LapackInt; ab: ptr cdouble; ldab: LapackInt; ipiv: ptr LapackInt;
           b: ptr cdouble; ldb: LapackInt): LapackInt {.importc: "LAPACKE_dgbsv",
    dynlib: libName.}
proc dgbsvWork*(matrixLayout: cint; n: LapackInt; kl: LapackInt; ku: LapackInt;
               nrhs: LapackInt; ab: ptr cdouble; ldab: LapackInt; ipiv: ptr LapackInt;
               b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dgbsv_work", dynlib: libName.}
proc dgbsvx*(matrixLayout: cint; fact: char; trans: char; n: LapackInt; kl: LapackInt;
            ku: LapackInt; nrhs: LapackInt; ab: ptr cdouble; ldab: LapackInt;
            afb: ptr cdouble; ldafb: LapackInt; ipiv: ptr LapackInt; equed: cstring;
            r: ptr cdouble; c: ptr cdouble; b: ptr cdouble; ldb: LapackInt; x: ptr cdouble;
            ldx: LapackInt; rcond: ptr cdouble; ferr: ptr cdouble; berr: ptr cdouble;
            rpivot: ptr cdouble): LapackInt {.importc: "LAPACKE_dgbsvx",
    dynlib: libName.}
proc dgbsvxWork*(matrixLayout: cint; fact: char; trans: char; n: LapackInt;
                kl: LapackInt; ku: LapackInt; nrhs: LapackInt; ab: ptr cdouble;
                ldab: LapackInt; afb: ptr cdouble; ldafb: LapackInt;
                ipiv: ptr LapackInt; equed: cstring; r: ptr cdouble; c: ptr cdouble;
                b: ptr cdouble; ldb: LapackInt; x: ptr cdouble; ldx: LapackInt;
                rcond: ptr cdouble; ferr: ptr cdouble; berr: ptr cdouble;
                work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgbsvx_work", dynlib: libName.}
proc dgbsvxx*(matrixLayout: cint; fact: char; trans: char; n: LapackInt; kl: LapackInt;
             ku: LapackInt; nrhs: LapackInt; ab: ptr cdouble; ldab: LapackInt;
             afb: ptr cdouble; ldafb: LapackInt; ipiv: ptr LapackInt; equed: cstring;
             r: ptr cdouble; c: ptr cdouble; b: ptr cdouble; ldb: LapackInt;
             x: ptr cdouble; ldx: LapackInt; rcond: ptr cdouble; rpvgrw: ptr cdouble;
             berr: ptr cdouble; nErrBnds: LapackInt; errBndsNorm: ptr cdouble;
             errBndsComp: ptr cdouble; nparams: LapackInt; params: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgbsvxx", dynlib: libName.}
proc dgbsvxxWork*(matrixLayout: cint; fact: char; trans: char; n: LapackInt;
                 kl: LapackInt; ku: LapackInt; nrhs: LapackInt; ab: ptr cdouble;
                 ldab: LapackInt; afb: ptr cdouble; ldafb: LapackInt;
                 ipiv: ptr LapackInt; equed: cstring; r: ptr cdouble; c: ptr cdouble;
                 b: ptr cdouble; ldb: LapackInt; x: ptr cdouble; ldx: LapackInt;
                 rcond: ptr cdouble; rpvgrw: ptr cdouble; berr: ptr cdouble;
                 nErrBnds: LapackInt; errBndsNorm: ptr cdouble;
                 errBndsComp: ptr cdouble; nparams: LapackInt; params: ptr cdouble;
                 work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgbsvxx_work", dynlib: libName.}
proc dgbtrf*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
            ku: LapackInt; ab: ptr cdouble; ldab: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgbtrf", dynlib: libName.}
proc dgbtrfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
                ku: LapackInt; ab: ptr cdouble; ldab: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgbtrf_work", dynlib: libName.}
proc dgbtrs*(matrixLayout: cint; trans: char; n: LapackInt; kl: LapackInt; ku: LapackInt;
            nrhs: LapackInt; ab: ptr cdouble; ldab: LapackInt; ipiv: ptr LapackInt;
            b: ptr cdouble; ldb: LapackInt): LapackInt {.importc: "LAPACKE_dgbtrs",
    dynlib: libName.}
proc dgbtrsWork*(matrixLayout: cint; trans: char; n: LapackInt; kl: LapackInt;
                ku: LapackInt; nrhs: LapackInt; ab: ptr cdouble; ldab: LapackInt;
                ipiv: ptr LapackInt; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dgbtrs_work", dynlib: libName.}
proc dgebak*(matrixLayout: cint; job: char; side: char; n: LapackInt; ilo: LapackInt;
            ihi: LapackInt; scale: ptr cdouble; m: LapackInt; v: ptr cdouble;
            ldv: LapackInt): LapackInt {.importc: "LAPACKE_dgebak", dynlib: libName.}
proc dgebakWork*(matrixLayout: cint; job: char; side: char; n: LapackInt; ilo: LapackInt;
                ihi: LapackInt; scale: ptr cdouble; m: LapackInt; v: ptr cdouble;
                ldv: LapackInt): LapackInt {.importc: "LAPACKE_dgebak_work",
    dynlib: libName.}
proc dgebal*(matrixLayout: cint; job: char; n: LapackInt; a: ptr cdouble; lda: LapackInt;
            ilo: ptr LapackInt; ihi: ptr LapackInt; scale: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgebal", dynlib: libName.}
proc dgebalWork*(matrixLayout: cint; job: char; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; ilo: ptr LapackInt; ihi: ptr LapackInt;
                scale: ptr cdouble): LapackInt {.importc: "LAPACKE_dgebal_work",
    dynlib: libName.}
proc dgebrd*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
            lda: LapackInt; d: ptr cdouble; e: ptr cdouble; tauq: ptr cdouble;
            taup: ptr cdouble): LapackInt {.importc: "LAPACKE_dgebrd", dynlib: libName.}
proc dgebrdWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; d: ptr cdouble; e: ptr cdouble; tauq: ptr cdouble;
                taup: ptr cdouble; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dgebrd_work", dynlib: libName.}
proc dgecon*(matrixLayout: cint; norm: char; n: LapackInt; a: ptr cdouble; lda: LapackInt;
            anorm: cdouble; rcond: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgecon", dynlib: libName.}
proc dgeconWork*(matrixLayout: cint; norm: char; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; anorm: cdouble; rcond: ptr cdouble; work: ptr cdouble;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dgecon_work",
    dynlib: libName.}
proc dgeequ*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
            lda: LapackInt; r: ptr cdouble; c: ptr cdouble; rowcnd: ptr cdouble;
            colcnd: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgeequ", dynlib: libName.}
proc dgeequWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; r: ptr cdouble; c: ptr cdouble; rowcnd: ptr cdouble;
                colcnd: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgeequ_work", dynlib: libName.}
proc dgeequb*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
             lda: LapackInt; r: ptr cdouble; c: ptr cdouble; rowcnd: ptr cdouble;
             colcnd: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgeequb", dynlib: libName.}
proc dgeequbWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
                 lda: LapackInt; r: ptr cdouble; c: ptr cdouble; rowcnd: ptr cdouble;
                 colcnd: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgeequb_work", dynlib: libName.}
proc dgees*(matrixLayout: cint; jobvs: char; sort: char; select: Lapack_D_Select2;
           n: LapackInt; a: ptr cdouble; lda: LapackInt; sdim: ptr LapackInt;
           wr: ptr cdouble; wi: ptr cdouble; vs: ptr cdouble; ldvs: LapackInt): LapackInt {.
    importc: "LAPACKE_dgees", dynlib: libName.}
proc dgeesWork*(matrixLayout: cint; jobvs: char; sort: char; select: Lapack_D_Select2;
               n: LapackInt; a: ptr cdouble; lda: LapackInt; sdim: ptr LapackInt;
               wr: ptr cdouble; wi: ptr cdouble; vs: ptr cdouble; ldvs: LapackInt;
               work: ptr cdouble; lwork: LapackInt; bwork: ptr LapackLogical): LapackInt {.
    importc: "LAPACKE_dgees_work", dynlib: libName.}
proc dgeesx*(matrixLayout: cint; jobvs: char; sort: char; select: Lapack_D_Select2;
            sense: char; n: LapackInt; a: ptr cdouble; lda: LapackInt;
            sdim: ptr LapackInt; wr: ptr cdouble; wi: ptr cdouble; vs: ptr cdouble;
            ldvs: LapackInt; rconde: ptr cdouble; rcondv: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgeesx", dynlib: libName.}
proc dgeesxWork*(matrixLayout: cint; jobvs: char; sort: char; select: Lapack_D_Select2;
                sense: char; n: LapackInt; a: ptr cdouble; lda: LapackInt;
                sdim: ptr LapackInt; wr: ptr cdouble; wi: ptr cdouble; vs: ptr cdouble;
                ldvs: LapackInt; rconde: ptr cdouble; rcondv: ptr cdouble;
                work: ptr cdouble; lwork: LapackInt; iwork: ptr LapackInt;
                liwork: LapackInt; bwork: ptr LapackLogical): LapackInt {.
    importc: "LAPACKE_dgeesx_work", dynlib: libName.}
proc dgeev*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt; a: ptr cdouble;
           lda: LapackInt; wr: ptr cdouble; wi: ptr cdouble; vl: ptr cdouble;
           ldvl: LapackInt; vr: ptr cdouble; ldvr: LapackInt): LapackInt {.
    importc: "LAPACKE_dgeev", dynlib: libName.}
proc dgeevWork*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt;
               a: ptr cdouble; lda: LapackInt; wr: ptr cdouble; wi: ptr cdouble;
               vl: ptr cdouble; ldvl: LapackInt; vr: ptr cdouble; ldvr: LapackInt;
               work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dgeev_work", dynlib: libName.}
proc dgeevx*(matrixLayout: cint; balanc: char; jobvl: char; jobvr: char; sense: char;
            n: LapackInt; a: ptr cdouble; lda: LapackInt; wr: ptr cdouble;
            wi: ptr cdouble; vl: ptr cdouble; ldvl: LapackInt; vr: ptr cdouble;
            ldvr: LapackInt; ilo: ptr LapackInt; ihi: ptr LapackInt; scale: ptr cdouble;
            abnrm: ptr cdouble; rconde: ptr cdouble; rcondv: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgeevx", dynlib: libName.}
proc dgeevxWork*(matrixLayout: cint; balanc: char; jobvl: char; jobvr: char; sense: char;
                n: LapackInt; a: ptr cdouble; lda: LapackInt; wr: ptr cdouble;
                wi: ptr cdouble; vl: ptr cdouble; ldvl: LapackInt; vr: ptr cdouble;
                ldvr: LapackInt; ilo: ptr LapackInt; ihi: ptr LapackInt;
                scale: ptr cdouble; abnrm: ptr cdouble; rconde: ptr cdouble;
                rcondv: ptr cdouble; work: ptr cdouble; lwork: LapackInt;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dgeevx_work",
    dynlib: libName.}
proc dgehrd*(matrixLayout: cint; n: LapackInt; ilo: LapackInt; ihi: LapackInt;
            a: ptr cdouble; lda: LapackInt; tau: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgehrd", dynlib: libName.}
proc dgehrdWork*(matrixLayout: cint; n: LapackInt; ilo: LapackInt; ihi: LapackInt;
                a: ptr cdouble; lda: LapackInt; tau: ptr cdouble; work: ptr cdouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_dgehrd_work",
    dynlib: libName.}
proc dgejsv*(matrixLayout: cint; joba: char; jobu: char; jobv: char; jobr: char;
            jobt: char; jobp: char; m: LapackInt; n: LapackInt; a: ptr cdouble;
            lda: LapackInt; sva: ptr cdouble; u: ptr cdouble; ldu: LapackInt;
            v: ptr cdouble; ldv: LapackInt; stat: ptr cdouble; istat: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgejsv", dynlib: libName.}
proc dgejsvWork*(matrixLayout: cint; joba: char; jobu: char; jobv: char; jobr: char;
                jobt: char; jobp: char; m: LapackInt; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; sva: ptr cdouble; u: ptr cdouble; ldu: LapackInt;
                v: ptr cdouble; ldv: LapackInt; work: ptr cdouble; lwork: LapackInt;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dgejsv_work",
    dynlib: libName.}
proc dgelq2*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
            lda: LapackInt; tau: ptr cdouble): LapackInt {.importc: "LAPACKE_dgelq2",
    dynlib: libName.}
proc dgelq2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; tau: ptr cdouble; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgelq2_work", dynlib: libName.}
proc dgelqf*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
            lda: LapackInt; tau: ptr cdouble): LapackInt {.importc: "LAPACKE_dgelqf",
    dynlib: libName.}
proc dgelqfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; tau: ptr cdouble; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dgelqf_work", dynlib: libName.}
proc dgels*(matrixLayout: cint; trans: char; m: LapackInt; n: LapackInt; nrhs: LapackInt;
           a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dgels", dynlib: libName.}
proc dgelsWork*(matrixLayout: cint; trans: char; m: LapackInt; n: LapackInt;
               nrhs: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble;
               ldb: LapackInt; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dgels_work", dynlib: libName.}
proc dgelsd*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
            a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
            s: ptr cdouble; rcond: cdouble; rank: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgelsd", dynlib: libName.}
proc dgelsdWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
                a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
                s: ptr cdouble; rcond: cdouble; rank: ptr LapackInt; work: ptr cdouble;
                lwork: LapackInt; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgelsd_work", dynlib: libName.}
proc dgelss*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
            a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
            s: ptr cdouble; rcond: cdouble; rank: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgelss", dynlib: libName.}
proc dgelssWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
                a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
                s: ptr cdouble; rcond: cdouble; rank: ptr LapackInt; work: ptr cdouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_dgelss_work",
    dynlib: libName.}
proc dgelsy*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
            a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
            jpvt: ptr LapackInt; rcond: cdouble; rank: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgelsy", dynlib: libName.}
proc dgelsyWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
                a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
                jpvt: ptr LapackInt; rcond: cdouble; rank: ptr LapackInt;
                work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dgelsy_work", dynlib: libName.}
proc dgemqrt*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
             k: LapackInt; nb: LapackInt; v: ptr cdouble; ldv: LapackInt; t: ptr cdouble;
             ldt: LapackInt; c: ptr cdouble; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_dgemqrt", dynlib: libName.}
proc dgemqrtWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt;
                 n: LapackInt; k: LapackInt; nb: LapackInt; v: ptr cdouble;
                 ldv: LapackInt; t: ptr cdouble; ldt: LapackInt; c: ptr cdouble;
                 ldc: LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgemqrt_work", dynlib: libName.}
proc dgeqlf*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
            lda: LapackInt; tau: ptr cdouble): LapackInt {.importc: "LAPACKE_dgeqlf",
    dynlib: libName.}
proc dgeqlfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; tau: ptr cdouble; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dgeqlf_work", dynlib: libName.}
proc dgeqp3*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
            lda: LapackInt; jpvt: ptr LapackInt; tau: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgeqp3", dynlib: libName.}
proc dgeqp3Work*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; jpvt: ptr LapackInt; tau: ptr cdouble;
                work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dgeqp3_work", dynlib: libName.}
proc dgeqpf*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
            lda: LapackInt; jpvt: ptr LapackInt; tau: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgeqpf", dynlib: libName.}
proc dgeqpfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; jpvt: ptr LapackInt; tau: ptr cdouble; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgeqpf_work", dynlib: libName.}
proc dgeqr2*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
            lda: LapackInt; tau: ptr cdouble): LapackInt {.importc: "LAPACKE_dgeqr2",
    dynlib: libName.}
proc dgeqr2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; tau: ptr cdouble; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgeqr2_work", dynlib: libName.}
proc dgeqrf*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
            lda: LapackInt; tau: ptr cdouble): LapackInt {.importc: "LAPACKE_dgeqrf",
    dynlib: libName.}
proc dgeqrfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; tau: ptr cdouble; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dgeqrf_work", dynlib: libName.}
proc dgeqrfp*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
             lda: LapackInt; tau: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgeqrfp", dynlib: libName.}
proc dgeqrfpWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
                 lda: LapackInt; tau: ptr cdouble; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dgeqrfp_work", dynlib: libName.}
proc dgeqrt*(matrixLayout: cint; m: LapackInt; n: LapackInt; nb: LapackInt;
            a: ptr cdouble; lda: LapackInt; t: ptr cdouble; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_dgeqrt", dynlib: libName.}
proc dgeqrt2*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
             lda: LapackInt; t: ptr cdouble; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_dgeqrt2", dynlib: libName.}
proc dgeqrt2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
                 lda: LapackInt; t: ptr cdouble; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_dgeqrt2_work", dynlib: libName.}
proc dgeqrt3*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
             lda: LapackInt; t: ptr cdouble; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_dgeqrt3", dynlib: libName.}
proc dgeqrt3Work*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
                 lda: LapackInt; t: ptr cdouble; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_dgeqrt3_work", dynlib: libName.}
proc dgeqrtWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; nb: LapackInt;
                a: ptr cdouble; lda: LapackInt; t: ptr cdouble; ldt: LapackInt;
                work: ptr cdouble): LapackInt {.importc: "LAPACKE_dgeqrt_work",
    dynlib: libName.}
proc dgerfs*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
            a: ptr cdouble; lda: LapackInt; af: ptr cdouble; ldaf: LapackInt;
            ipiv: ptr LapackInt; b: ptr cdouble; ldb: LapackInt; x: ptr cdouble;
            ldx: LapackInt; ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgerfs", dynlib: libName.}
proc dgerfsWork*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
                a: ptr cdouble; lda: LapackInt; af: ptr cdouble; ldaf: LapackInt;
                ipiv: ptr LapackInt; b: ptr cdouble; ldb: LapackInt; x: ptr cdouble;
                ldx: LapackInt; ferr: ptr cdouble; berr: ptr cdouble; work: ptr cdouble;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dgerfs_work",
    dynlib: libName.}
proc dgerfsx*(matrixLayout: cint; trans: char; equed: char; n: LapackInt;
             nrhs: LapackInt; a: ptr cdouble; lda: LapackInt; af: ptr cdouble;
             ldaf: LapackInt; ipiv: ptr LapackInt; r: ptr cdouble; c: ptr cdouble;
             b: ptr cdouble; ldb: LapackInt; x: ptr cdouble; ldx: LapackInt;
             rcond: ptr cdouble; berr: ptr cdouble; nErrBnds: LapackInt;
             errBndsNorm: ptr cdouble; errBndsComp: ptr cdouble; nparams: LapackInt;
             params: ptr cdouble): LapackInt {.importc: "LAPACKE_dgerfsx",
    dynlib: libName.}
proc dgerfsxWork*(matrixLayout: cint; trans: char; equed: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr cdouble; lda: LapackInt; af: ptr cdouble;
                 ldaf: LapackInt; ipiv: ptr LapackInt; r: ptr cdouble; c: ptr cdouble;
                 b: ptr cdouble; ldb: LapackInt; x: ptr cdouble; ldx: LapackInt;
                 rcond: ptr cdouble; berr: ptr cdouble; nErrBnds: LapackInt;
                 errBndsNorm: ptr cdouble; errBndsComp: ptr cdouble;
                 nparams: LapackInt; params: ptr cdouble; work: ptr cdouble;
                 iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dgerfsx_work",
    dynlib: libName.}
proc dgerqf*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
            lda: LapackInt; tau: ptr cdouble): LapackInt {.importc: "LAPACKE_dgerqf",
    dynlib: libName.}
proc dgerqfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; tau: ptr cdouble; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dgerqf_work", dynlib: libName.}
proc dgesdd*(matrixLayout: cint; jobz: char; m: LapackInt; n: LapackInt; a: ptr cdouble;
            lda: LapackInt; s: ptr cdouble; u: ptr cdouble; ldu: LapackInt;
            vt: ptr cdouble; ldvt: LapackInt): LapackInt {.importc: "LAPACKE_dgesdd",
    dynlib: libName.}
proc dgesddWork*(matrixLayout: cint; jobz: char; m: LapackInt; n: LapackInt;
                a: ptr cdouble; lda: LapackInt; s: ptr cdouble; u: ptr cdouble;
                ldu: LapackInt; vt: ptr cdouble; ldvt: LapackInt; work: ptr cdouble;
                lwork: LapackInt; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgesdd_work", dynlib: libName.}
proc dgesv*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; a: ptr cdouble;
           lda: LapackInt; ipiv: ptr LapackInt; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dgesv", dynlib: libName.}
proc dgesvWork*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; a: ptr cdouble;
               lda: LapackInt; ipiv: ptr LapackInt; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dgesv_work", dynlib: libName.}
proc dgesvd*(matrixLayout: cint; jobu: char; jobvt: char; m: LapackInt; n: LapackInt;
            a: ptr cdouble; lda: LapackInt; s: ptr cdouble; u: ptr cdouble; ldu: LapackInt;
            vt: ptr cdouble; ldvt: LapackInt; superb: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgesvd", dynlib: libName.}
proc dgesvdWork*(matrixLayout: cint; jobu: char; jobvt: char; m: LapackInt; n: LapackInt;
                a: ptr cdouble; lda: LapackInt; s: ptr cdouble; u: ptr cdouble;
                ldu: LapackInt; vt: ptr cdouble; ldvt: LapackInt; work: ptr cdouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_dgesvd_work",
    dynlib: libName.}
proc dgesvdx*(matrixLayout: cint; jobu: char; jobvt: char; range: char; m: LapackInt;
             n: LapackInt; a: ptr cdouble; lda: LapackInt; vl: cdouble; vu: cdouble;
             il: LapackInt; iu: LapackInt; ns: ptr LapackInt; s: ptr cdouble;
             u: ptr cdouble; ldu: LapackInt; vt: ptr cdouble; ldvt: LapackInt;
             superb: ptr LapackInt): LapackInt {.importc: "LAPACKE_dgesvdx",
    dynlib: libName.}
proc dgesvdxWork*(matrixLayout: cint; jobu: char; jobvt: char; range: char; m: LapackInt;
                 n: LapackInt; a: ptr cdouble; lda: LapackInt; vl: cdouble; vu: cdouble;
                 il: LapackInt; iu: LapackInt; ns: ptr LapackInt; s: ptr cdouble;
                 u: ptr cdouble; ldu: LapackInt; vt: ptr cdouble; ldvt: LapackInt;
                 work: ptr cdouble; lwork: LapackInt; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgesvdx_work", dynlib: libName.}
proc dgesvj*(matrixLayout: cint; joba: char; jobu: char; jobv: char; m: LapackInt;
            n: LapackInt; a: ptr cdouble; lda: LapackInt; sva: ptr cdouble; mv: LapackInt;
            v: ptr cdouble; ldv: LapackInt; stat: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgesvj", dynlib: libName.}
proc dgesvjWork*(matrixLayout: cint; joba: char; jobu: char; jobv: char; m: LapackInt;
                n: LapackInt; a: ptr cdouble; lda: LapackInt; sva: ptr cdouble;
                mv: LapackInt; v: ptr cdouble; ldv: LapackInt; work: ptr cdouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_dgesvj_work",
    dynlib: libName.}
proc dgesvx*(matrixLayout: cint; fact: char; trans: char; n: LapackInt; nrhs: LapackInt;
            a: ptr cdouble; lda: LapackInt; af: ptr cdouble; ldaf: LapackInt;
            ipiv: ptr LapackInt; equed: cstring; r: ptr cdouble; c: ptr cdouble;
            b: ptr cdouble; ldb: LapackInt; x: ptr cdouble; ldx: LapackInt;
            rcond: ptr cdouble; ferr: ptr cdouble; berr: ptr cdouble; rpivot: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgesvx", dynlib: libName.}
proc dgesvxWork*(matrixLayout: cint; fact: char; trans: char; n: LapackInt;
                nrhs: LapackInt; a: ptr cdouble; lda: LapackInt; af: ptr cdouble;
                ldaf: LapackInt; ipiv: ptr LapackInt; equed: cstring; r: ptr cdouble;
                c: ptr cdouble; b: ptr cdouble; ldb: LapackInt; x: ptr cdouble;
                ldx: LapackInt; rcond: ptr cdouble; ferr: ptr cdouble;
                berr: ptr cdouble; work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgesvx_work", dynlib: libName.}
proc dgesvxx*(matrixLayout: cint; fact: char; trans: char; n: LapackInt; nrhs: LapackInt;
             a: ptr cdouble; lda: LapackInt; af: ptr cdouble; ldaf: LapackInt;
             ipiv: ptr LapackInt; equed: cstring; r: ptr cdouble; c: ptr cdouble;
             b: ptr cdouble; ldb: LapackInt; x: ptr cdouble; ldx: LapackInt;
             rcond: ptr cdouble; rpvgrw: ptr cdouble; berr: ptr cdouble;
             nErrBnds: LapackInt; errBndsNorm: ptr cdouble; errBndsComp: ptr cdouble;
             nparams: LapackInt; params: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgesvxx", dynlib: libName.}
proc dgesvxxWork*(matrixLayout: cint; fact: char; trans: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr cdouble; lda: LapackInt; af: ptr cdouble;
                 ldaf: LapackInt; ipiv: ptr LapackInt; equed: cstring; r: ptr cdouble;
                 c: ptr cdouble; b: ptr cdouble; ldb: LapackInt; x: ptr cdouble;
                 ldx: LapackInt; rcond: ptr cdouble; rpvgrw: ptr cdouble;
                 berr: ptr cdouble; nErrBnds: LapackInt; errBndsNorm: ptr cdouble;
                 errBndsComp: ptr cdouble; nparams: LapackInt; params: ptr cdouble;
                 work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgesvxx_work", dynlib: libName.}
proc dgetf2*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
            lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgetf2", dynlib: libName.}
proc dgetf2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgetf2_work", dynlib: libName.}
proc dgetrf*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
            lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgetrf", dynlib: libName.}
proc dgetrf2*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
             lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgetrf2", dynlib: libName.}
proc dgetrf2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
                 lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgetrf2_work", dynlib: libName.}
proc dgetrfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgetrf_work", dynlib: libName.}
proc dgetri*(matrixLayout: cint; n: LapackInt; a: ptr cdouble; lda: LapackInt;
            ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_dgetri",
    dynlib: libName.}
proc dgetriWork*(matrixLayout: cint; n: LapackInt; a: ptr cdouble; lda: LapackInt;
                ipiv: ptr LapackInt; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dgetri_work", dynlib: libName.}
proc dgetrs*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
            a: ptr cdouble; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cdouble;
            ldb: LapackInt): LapackInt {.importc: "LAPACKE_dgetrs", dynlib: libName.}
proc dgetrsWork*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
                a: ptr cdouble; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cdouble;
                ldb: LapackInt): LapackInt {.importc: "LAPACKE_dgetrs_work",
    dynlib: libName.}
proc dggbak*(matrixLayout: cint; job: char; side: char; n: LapackInt; ilo: LapackInt;
            ihi: LapackInt; lscale: ptr cdouble; rscale: ptr cdouble; m: LapackInt;
            v: ptr cdouble; ldv: LapackInt): LapackInt {.importc: "LAPACKE_dggbak",
    dynlib: libName.}
proc dggbakWork*(matrixLayout: cint; job: char; side: char; n: LapackInt; ilo: LapackInt;
                ihi: LapackInt; lscale: ptr cdouble; rscale: ptr cdouble; m: LapackInt;
                v: ptr cdouble; ldv: LapackInt): LapackInt {.
    importc: "LAPACKE_dggbak_work", dynlib: libName.}
proc dggbal*(matrixLayout: cint; job: char; n: LapackInt; a: ptr cdouble; lda: LapackInt;
            b: ptr cdouble; ldb: LapackInt; ilo: ptr LapackInt; ihi: ptr LapackInt;
            lscale: ptr cdouble; rscale: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dggbal", dynlib: libName.}
proc dggbalWork*(matrixLayout: cint; job: char; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; b: ptr cdouble; ldb: LapackInt; ilo: ptr LapackInt;
                ihi: ptr LapackInt; lscale: ptr cdouble; rscale: ptr cdouble;
                work: ptr cdouble): LapackInt {.importc: "LAPACKE_dggbal_work",
    dynlib: libName.}
proc dgges*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
           selctg: Lapack_D_Select3; n: LapackInt; a: ptr cdouble; lda: LapackInt;
           b: ptr cdouble; ldb: LapackInt; sdim: ptr LapackInt; alphar: ptr cdouble;
           alphai: ptr cdouble; beta: ptr cdouble; vsl: ptr cdouble; ldvsl: LapackInt;
           vsr: ptr cdouble; ldvsr: LapackInt): LapackInt {.importc: "LAPACKE_dgges",
    dynlib: libName.}
proc dgges3*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
            selctg: Lapack_D_Select3; n: LapackInt; a: ptr cdouble; lda: LapackInt;
            b: ptr cdouble; ldb: LapackInt; sdim: ptr LapackInt; alphar: ptr cdouble;
            alphai: ptr cdouble; beta: ptr cdouble; vsl: ptr cdouble; ldvsl: LapackInt;
            vsr: ptr cdouble; ldvsr: LapackInt): LapackInt {.
    importc: "LAPACKE_dgges3", dynlib: libName.}
proc dgges3Work*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
                selctg: Lapack_D_Select3; n: LapackInt; a: ptr cdouble; lda: LapackInt;
                b: ptr cdouble; ldb: LapackInt; sdim: ptr LapackInt;
                alphar: ptr cdouble; alphai: ptr cdouble; beta: ptr cdouble;
                vsl: ptr cdouble; ldvsl: LapackInt; vsr: ptr cdouble; ldvsr: LapackInt;
                work: ptr cdouble; lwork: LapackInt; bwork: ptr LapackLogical): LapackInt {.
    importc: "LAPACKE_dgges3_work", dynlib: libName.}
proc dggesWork*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
               selctg: Lapack_D_Select3; n: LapackInt; a: ptr cdouble; lda: LapackInt;
               b: ptr cdouble; ldb: LapackInt; sdim: ptr LapackInt; alphar: ptr cdouble;
               alphai: ptr cdouble; beta: ptr cdouble; vsl: ptr cdouble;
               ldvsl: LapackInt; vsr: ptr cdouble; ldvsr: LapackInt; work: ptr cdouble;
               lwork: LapackInt; bwork: ptr LapackLogical): LapackInt {.
    importc: "LAPACKE_dgges_work", dynlib: libName.}
proc dggesx*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
            selctg: Lapack_D_Select3; sense: char; n: LapackInt; a: ptr cdouble;
            lda: LapackInt; b: ptr cdouble; ldb: LapackInt; sdim: ptr LapackInt;
            alphar: ptr cdouble; alphai: ptr cdouble; beta: ptr cdouble;
            vsl: ptr cdouble; ldvsl: LapackInt; vsr: ptr cdouble; ldvsr: LapackInt;
            rconde: ptr cdouble; rcondv: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dggesx", dynlib: libName.}
proc dggesxWork*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
                selctg: Lapack_D_Select3; sense: char; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; b: ptr cdouble; ldb: LapackInt; sdim: ptr LapackInt;
                alphar: ptr cdouble; alphai: ptr cdouble; beta: ptr cdouble;
                vsl: ptr cdouble; ldvsl: LapackInt; vsr: ptr cdouble; ldvsr: LapackInt;
                rconde: ptr cdouble; rcondv: ptr cdouble; work: ptr cdouble;
                lwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt;
                bwork: ptr LapackLogical): LapackInt {.
    importc: "LAPACKE_dggesx_work", dynlib: libName.}
proc dggev*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt; a: ptr cdouble;
           lda: LapackInt; b: ptr cdouble; ldb: LapackInt; alphar: ptr cdouble;
           alphai: ptr cdouble; beta: ptr cdouble; vl: ptr cdouble; ldvl: LapackInt;
           vr: ptr cdouble; ldvr: LapackInt): LapackInt {.importc: "LAPACKE_dggev",
    dynlib: libName.}
proc dggev3*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt; a: ptr cdouble;
            lda: LapackInt; b: ptr cdouble; ldb: LapackInt; alphar: ptr cdouble;
            alphai: ptr cdouble; beta: ptr cdouble; vl: ptr cdouble; ldvl: LapackInt;
            vr: ptr cdouble; ldvr: LapackInt): LapackInt {.importc: "LAPACKE_dggev3",
    dynlib: libName.}
proc dggev3Work*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt;
                a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
                alphar: ptr cdouble; alphai: ptr cdouble; beta: ptr cdouble;
                vl: ptr cdouble; ldvl: LapackInt; vr: ptr cdouble; ldvr: LapackInt;
                work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dggev3_work", dynlib: libName.}
proc dggevWork*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt;
               a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
               alphar: ptr cdouble; alphai: ptr cdouble; beta: ptr cdouble;
               vl: ptr cdouble; ldvl: LapackInt; vr: ptr cdouble; ldvr: LapackInt;
               work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dggev_work", dynlib: libName.}
proc dggevx*(matrixLayout: cint; balanc: char; jobvl: char; jobvr: char; sense: char;
            n: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
            alphar: ptr cdouble; alphai: ptr cdouble; beta: ptr cdouble; vl: ptr cdouble;
            ldvl: LapackInt; vr: ptr cdouble; ldvr: LapackInt; ilo: ptr LapackInt;
            ihi: ptr LapackInt; lscale: ptr cdouble; rscale: ptr cdouble;
            abnrm: ptr cdouble; bbnrm: ptr cdouble; rconde: ptr cdouble;
            rcondv: ptr cdouble): LapackInt {.importc: "LAPACKE_dggevx",
    dynlib: libName.}
proc dggevxWork*(matrixLayout: cint; balanc: char; jobvl: char; jobvr: char; sense: char;
                n: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble;
                ldb: LapackInt; alphar: ptr cdouble; alphai: ptr cdouble;
                beta: ptr cdouble; vl: ptr cdouble; ldvl: LapackInt; vr: ptr cdouble;
                ldvr: LapackInt; ilo: ptr LapackInt; ihi: ptr LapackInt;
                lscale: ptr cdouble; rscale: ptr cdouble; abnrm: ptr cdouble;
                bbnrm: ptr cdouble; rconde: ptr cdouble; rcondv: ptr cdouble;
                work: ptr cdouble; lwork: LapackInt; iwork: ptr LapackInt;
                bwork: ptr LapackLogical): LapackInt {.
    importc: "LAPACKE_dggevx_work", dynlib: libName.}
proc dggglm*(matrixLayout: cint; n: LapackInt; m: LapackInt; p: LapackInt;
            a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
            d: ptr cdouble; x: ptr cdouble; y: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dggglm", dynlib: libName.}
proc dggglmWork*(matrixLayout: cint; n: LapackInt; m: LapackInt; p: LapackInt;
                a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
                d: ptr cdouble; x: ptr cdouble; y: ptr cdouble; work: ptr cdouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_dggglm_work",
    dynlib: libName.}
proc dgghd3*(matrixLayout: cint; compq: char; compz: char; n: LapackInt; ilo: LapackInt;
            ihi: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble;
            ldb: LapackInt; q: ptr cdouble; ldq: LapackInt; z: ptr cdouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_dgghd3", dynlib: libName.}
proc dgghd3Work*(matrixLayout: cint; compq: char; compz: char; n: LapackInt;
                ilo: LapackInt; ihi: LapackInt; a: ptr cdouble; lda: LapackInt;
                b: ptr cdouble; ldb: LapackInt; q: ptr cdouble; ldq: LapackInt;
                z: ptr cdouble; ldz: LapackInt; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dgghd3_work", dynlib: libName.}
proc dgghrd*(matrixLayout: cint; compq: char; compz: char; n: LapackInt; ilo: LapackInt;
            ihi: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble;
            ldb: LapackInt; q: ptr cdouble; ldq: LapackInt; z: ptr cdouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_dgghrd", dynlib: libName.}
proc dgghrdWork*(matrixLayout: cint; compq: char; compz: char; n: LapackInt;
                ilo: LapackInt; ihi: LapackInt; a: ptr cdouble; lda: LapackInt;
                b: ptr cdouble; ldb: LapackInt; q: ptr cdouble; ldq: LapackInt;
                z: ptr cdouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_dgghrd_work", dynlib: libName.}
proc dgglse*(matrixLayout: cint; m: LapackInt; n: LapackInt; p: LapackInt;
            a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
            c: ptr cdouble; d: ptr cdouble; x: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgglse", dynlib: libName.}
proc dgglseWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; p: LapackInt;
                a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
                c: ptr cdouble; d: ptr cdouble; x: ptr cdouble; work: ptr cdouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_dgglse_work",
    dynlib: libName.}
proc dggqrf*(matrixLayout: cint; n: LapackInt; m: LapackInt; p: LapackInt;
            a: ptr cdouble; lda: LapackInt; taua: ptr cdouble; b: ptr cdouble;
            ldb: LapackInt; taub: ptr cdouble): LapackInt {.importc: "LAPACKE_dggqrf",
    dynlib: libName.}
proc dggqrfWork*(matrixLayout: cint; n: LapackInt; m: LapackInt; p: LapackInt;
                a: ptr cdouble; lda: LapackInt; taua: ptr cdouble; b: ptr cdouble;
                ldb: LapackInt; taub: ptr cdouble; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dggqrf_work", dynlib: libName.}
proc dggrqf*(matrixLayout: cint; m: LapackInt; p: LapackInt; n: LapackInt;
            a: ptr cdouble; lda: LapackInt; taua: ptr cdouble; b: ptr cdouble;
            ldb: LapackInt; taub: ptr cdouble): LapackInt {.importc: "LAPACKE_dggrqf",
    dynlib: libName.}
proc dggrqfWork*(matrixLayout: cint; m: LapackInt; p: LapackInt; n: LapackInt;
                a: ptr cdouble; lda: LapackInt; taua: ptr cdouble; b: ptr cdouble;
                ldb: LapackInt; taub: ptr cdouble; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dggrqf_work", dynlib: libName.}
proc dggsvd*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
            n: LapackInt; p: LapackInt; k: ptr LapackInt; l: ptr LapackInt;
            a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
            alpha: ptr cdouble; beta: ptr cdouble; u: ptr cdouble; ldu: LapackInt;
            v: ptr cdouble; ldv: LapackInt; q: ptr cdouble; ldq: LapackInt;
            iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dggsvd",
    dynlib: libName.}
proc dggsvd3*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
             n: LapackInt; p: LapackInt; k: ptr LapackInt; l: ptr LapackInt;
             a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
             alpha: ptr cdouble; beta: ptr cdouble; u: ptr cdouble; ldu: LapackInt;
             v: ptr cdouble; ldv: LapackInt; q: ptr cdouble; ldq: LapackInt;
             iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dggsvd3",
    dynlib: libName.}
proc dggsvd3Work*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
                 n: LapackInt; p: LapackInt; k: ptr LapackInt; l: ptr LapackInt;
                 a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
                 alpha: ptr cdouble; beta: ptr cdouble; u: ptr cdouble; ldu: LapackInt;
                 v: ptr cdouble; ldv: LapackInt; q: ptr cdouble; ldq: LapackInt;
                 work: ptr cdouble; lwork: LapackInt; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dggsvd3_work", dynlib: libName.}
proc dggsvdWork*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
                n: LapackInt; p: LapackInt; k: ptr LapackInt; l: ptr LapackInt;
                a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
                alpha: ptr cdouble; beta: ptr cdouble; u: ptr cdouble; ldu: LapackInt;
                v: ptr cdouble; ldv: LapackInt; q: ptr cdouble; ldq: LapackInt;
                work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dggsvd_work", dynlib: libName.}
proc dggsvp*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
            p: LapackInt; n: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble;
            ldb: LapackInt; tola: cdouble; tolb: cdouble; k: ptr LapackInt;
            l: ptr LapackInt; u: ptr cdouble; ldu: LapackInt; v: ptr cdouble;
            ldv: LapackInt; q: ptr cdouble; ldq: LapackInt): LapackInt {.
    importc: "LAPACKE_dggsvp", dynlib: libName.}
proc dggsvp3*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
             p: LapackInt; n: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble;
             ldb: LapackInt; tola: cdouble; tolb: cdouble; k: ptr LapackInt;
             l: ptr LapackInt; u: ptr cdouble; ldu: LapackInt; v: ptr cdouble;
             ldv: LapackInt; q: ptr cdouble; ldq: LapackInt): LapackInt {.
    importc: "LAPACKE_dggsvp3", dynlib: libName.}
proc dggsvp3Work*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
                 p: LapackInt; n: LapackInt; a: ptr cdouble; lda: LapackInt;
                 b: ptr cdouble; ldb: LapackInt; tola: cdouble; tolb: cdouble;
                 k: ptr LapackInt; l: ptr LapackInt; u: ptr cdouble; ldu: LapackInt;
                 v: ptr cdouble; ldv: LapackInt; q: ptr cdouble; ldq: LapackInt;
                 iwork: ptr LapackInt; tau: ptr cdouble; work: ptr cdouble;
                 lwork: LapackInt): LapackInt {.importc: "LAPACKE_dggsvp3_work",
    dynlib: libName.}
proc dggsvpWork*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
                p: LapackInt; n: LapackInt; a: ptr cdouble; lda: LapackInt;
                b: ptr cdouble; ldb: LapackInt; tola: cdouble; tolb: cdouble;
                k: ptr LapackInt; l: ptr LapackInt; u: ptr cdouble; ldu: LapackInt;
                v: ptr cdouble; ldv: LapackInt; q: ptr cdouble; ldq: LapackInt;
                iwork: ptr LapackInt; tau: ptr cdouble; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dggsvp_work", dynlib: libName.}
proc dgtcon*(norm: char; n: LapackInt; dl: ptr cdouble; d: ptr cdouble; du: ptr cdouble;
            du2: ptr cdouble; ipiv: ptr LapackInt; anorm: cdouble; rcond: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgtcon", dynlib: libName.}
proc dgtconWork*(norm: char; n: LapackInt; dl: ptr cdouble; d: ptr cdouble;
                du: ptr cdouble; du2: ptr cdouble; ipiv: ptr LapackInt; anorm: cdouble;
                rcond: ptr cdouble; work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgtcon_work", dynlib: libName.}
proc dgtrfs*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
            dl: ptr cdouble; d: ptr cdouble; du: ptr cdouble; dlf: ptr cdouble;
            df: ptr cdouble; duf: ptr cdouble; du2: ptr cdouble; ipiv: ptr LapackInt;
            b: ptr cdouble; ldb: LapackInt; x: ptr cdouble; ldx: LapackInt;
            ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgtrfs", dynlib: libName.}
proc dgtrfsWork*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
                dl: ptr cdouble; d: ptr cdouble; du: ptr cdouble; dlf: ptr cdouble;
                df: ptr cdouble; duf: ptr cdouble; du2: ptr cdouble; ipiv: ptr LapackInt;
                b: ptr cdouble; ldb: LapackInt; x: ptr cdouble; ldx: LapackInt;
                ferr: ptr cdouble; berr: ptr cdouble; work: ptr cdouble;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dgtrfs_work",
    dynlib: libName.}
proc dgtsv*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; dl: ptr cdouble;
           d: ptr cdouble; du: ptr cdouble; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dgtsv", dynlib: libName.}
proc dgtsvWork*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; dl: ptr cdouble;
               d: ptr cdouble; du: ptr cdouble; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dgtsv_work", dynlib: libName.}
proc dgtsvx*(matrixLayout: cint; fact: char; trans: char; n: LapackInt; nrhs: LapackInt;
            dl: ptr cdouble; d: ptr cdouble; du: ptr cdouble; dlf: ptr cdouble;
            df: ptr cdouble; duf: ptr cdouble; du2: ptr cdouble; ipiv: ptr LapackInt;
            b: ptr cdouble; ldb: LapackInt; x: ptr cdouble; ldx: LapackInt;
            rcond: ptr cdouble; ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dgtsvx", dynlib: libName.}
proc dgtsvxWork*(matrixLayout: cint; fact: char; trans: char; n: LapackInt;
                nrhs: LapackInt; dl: ptr cdouble; d: ptr cdouble; du: ptr cdouble;
                dlf: ptr cdouble; df: ptr cdouble; duf: ptr cdouble; du2: ptr cdouble;
                ipiv: ptr LapackInt; b: ptr cdouble; ldb: LapackInt; x: ptr cdouble;
                ldx: LapackInt; rcond: ptr cdouble; ferr: ptr cdouble;
                berr: ptr cdouble; work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgtsvx_work", dynlib: libName.}
proc dgttrf*(n: LapackInt; dl: ptr cdouble; d: ptr cdouble; du: ptr cdouble;
            du2: ptr cdouble; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgttrf", dynlib: libName.}
proc dgttrfWork*(n: LapackInt; dl: ptr cdouble; d: ptr cdouble; du: ptr cdouble;
                du2: ptr cdouble; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dgttrf_work", dynlib: libName.}
proc dgttrs*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
            dl: ptr cdouble; d: ptr cdouble; du: ptr cdouble; du2: ptr cdouble;
            ipiv: ptr LapackInt; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dgttrs", dynlib: libName.}
proc dgttrsWork*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
                dl: ptr cdouble; d: ptr cdouble; du: ptr cdouble; du2: ptr cdouble;
                ipiv: ptr LapackInt; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dgttrs_work", dynlib: libName.}
proc dhgeqz*(matrixLayout: cint; job: char; compq: char; compz: char; n: LapackInt;
            ilo: LapackInt; ihi: LapackInt; h: ptr cdouble; ldh: LapackInt;
            t: ptr cdouble; ldt: LapackInt; alphar: ptr cdouble; alphai: ptr cdouble;
            beta: ptr cdouble; q: ptr cdouble; ldq: LapackInt; z: ptr cdouble;
            ldz: LapackInt): LapackInt {.importc: "LAPACKE_dhgeqz", dynlib: libName.}
proc dhgeqzWork*(matrixLayout: cint; job: char; compq: char; compz: char; n: LapackInt;
                ilo: LapackInt; ihi: LapackInt; h: ptr cdouble; ldh: LapackInt;
                t: ptr cdouble; ldt: LapackInt; alphar: ptr cdouble;
                alphai: ptr cdouble; beta: ptr cdouble; q: ptr cdouble; ldq: LapackInt;
                z: ptr cdouble; ldz: LapackInt; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dhgeqz_work", dynlib: libName.}
proc dhsein*(matrixLayout: cint; job: char; eigsrc: char; initv: char;
            select: ptr LapackLogical; n: LapackInt; h: ptr cdouble; ldh: LapackInt;
            wr: ptr cdouble; wi: ptr cdouble; vl: ptr cdouble; ldvl: LapackInt;
            vr: ptr cdouble; ldvr: LapackInt; mm: LapackInt; m: ptr LapackInt;
            ifaill: ptr LapackInt; ifailr: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dhsein", dynlib: libName.}
proc dhseinWork*(matrixLayout: cint; job: char; eigsrc: char; initv: char;
                select: ptr LapackLogical; n: LapackInt; h: ptr cdouble; ldh: LapackInt;
                wr: ptr cdouble; wi: ptr cdouble; vl: ptr cdouble; ldvl: LapackInt;
                vr: ptr cdouble; ldvr: LapackInt; mm: LapackInt; m: ptr LapackInt;
                work: ptr cdouble; ifaill: ptr LapackInt; ifailr: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dhsein_work", dynlib: libName.}
proc dhseqr*(matrixLayout: cint; job: char; compz: char; n: LapackInt; ilo: LapackInt;
            ihi: LapackInt; h: ptr cdouble; ldh: LapackInt; wr: ptr cdouble;
            wi: ptr cdouble; z: ptr cdouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_dhseqr", dynlib: libName.}
proc dhseqrWork*(matrixLayout: cint; job: char; compz: char; n: LapackInt;
                ilo: LapackInt; ihi: LapackInt; h: ptr cdouble; ldh: LapackInt;
                wr: ptr cdouble; wi: ptr cdouble; z: ptr cdouble; ldz: LapackInt;
                work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dhseqr_work", dynlib: libName.}
proc dlacn2*(n: LapackInt; v: ptr cdouble; x: ptr cdouble; isgn: ptr LapackInt;
            est: ptr cdouble; kase: ptr LapackInt; isave: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dlacn2", dynlib: libName.}
proc dlacn2Work*(n: LapackInt; v: ptr cdouble; x: ptr cdouble; isgn: ptr LapackInt;
                est: ptr cdouble; kase: ptr LapackInt; isave: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dlacn2_work", dynlib: libName.}
proc dlacpy*(matrixLayout: cint; uplo: char; m: LapackInt; n: LapackInt; a: ptr cdouble;
            lda: LapackInt; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dlacpy", dynlib: libName.}
proc dlacpyWork*(matrixLayout: cint; uplo: char; m: LapackInt; n: LapackInt;
                a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dlacpy_work", dynlib: libName.}
proc dlag2s*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
            lda: LapackInt; sa: ptr cfloat; ldsa: LapackInt): LapackInt {.
    importc: "LAPACKE_dlag2s", dynlib: libName.}
proc dlag2sWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; sa: ptr cfloat; ldsa: LapackInt): LapackInt {.
    importc: "LAPACKE_dlag2s_work", dynlib: libName.}
proc dlagge*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
            ku: LapackInt; d: ptr cdouble; a: ptr cdouble; lda: LapackInt;
            iseed: ptr LapackInt): LapackInt {.importc: "LAPACKE_dlagge",
    dynlib: libName.}
proc dlaggeWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
                ku: LapackInt; d: ptr cdouble; a: ptr cdouble; lda: LapackInt;
                iseed: ptr LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dlagge_work", dynlib: libName.}
proc dlagsy*(matrixLayout: cint; n: LapackInt; k: LapackInt; d: ptr cdouble;
            a: ptr cdouble; lda: LapackInt; iseed: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dlagsy", dynlib: libName.}
proc dlagsyWork*(matrixLayout: cint; n: LapackInt; k: LapackInt; d: ptr cdouble;
                a: ptr cdouble; lda: LapackInt; iseed: ptr LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dlagsy_work", dynlib: libName.}
proc dlapmr*(matrixLayout: cint; forwrd: LapackLogical; m: LapackInt; n: LapackInt;
            x: ptr cdouble; ldx: LapackInt; k: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dlapmr", dynlib: libName.}
proc dlapmrWork*(matrixLayout: cint; forwrd: LapackLogical; m: LapackInt; n: LapackInt;
                x: ptr cdouble; ldx: LapackInt; k: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dlapmr_work", dynlib: libName.}
proc dlapmt*(matrixLayout: cint; forwrd: LapackLogical; m: LapackInt; n: LapackInt;
            x: ptr cdouble; ldx: LapackInt; k: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dlapmt", dynlib: libName.}
proc dlapmtWork*(matrixLayout: cint; forwrd: LapackLogical; m: LapackInt; n: LapackInt;
                x: ptr cdouble; ldx: LapackInt; k: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dlapmt_work", dynlib: libName.}
proc dlarfb*(matrixLayout: cint; side: char; trans: char; direct: char; storev: char;
            m: LapackInt; n: LapackInt; k: LapackInt; v: ptr cdouble; ldv: LapackInt;
            t: ptr cdouble; ldt: LapackInt; c: ptr cdouble; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_dlarfb", dynlib: libName.}
proc dlarfbWork*(matrixLayout: cint; side: char; trans: char; direct: char; storev: char;
                m: LapackInt; n: LapackInt; k: LapackInt; v: ptr cdouble; ldv: LapackInt;
                t: ptr cdouble; ldt: LapackInt; c: ptr cdouble; ldc: LapackInt;
                work: ptr cdouble; ldwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dlarfb_work", dynlib: libName.}
proc dlarfg*(n: LapackInt; alpha: ptr cdouble; x: ptr cdouble; incx: LapackInt;
            tau: ptr cdouble): LapackInt {.importc: "LAPACKE_dlarfg", dynlib: libName.}
proc dlarfgWork*(n: LapackInt; alpha: ptr cdouble; x: ptr cdouble; incx: LapackInt;
                tau: ptr cdouble): LapackInt {.importc: "LAPACKE_dlarfg_work",
    dynlib: libName.}
proc dlarft*(matrixLayout: cint; direct: char; storev: char; n: LapackInt; k: LapackInt;
            v: ptr cdouble; ldv: LapackInt; tau: ptr cdouble; t: ptr cdouble;
            ldt: LapackInt): LapackInt {.importc: "LAPACKE_dlarft", dynlib: libName.}
proc dlarftWork*(matrixLayout: cint; direct: char; storev: char; n: LapackInt;
                k: LapackInt; v: ptr cdouble; ldv: LapackInt; tau: ptr cdouble;
                t: ptr cdouble; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_dlarft_work", dynlib: libName.}
proc dlarfx*(matrixLayout: cint; side: char; m: LapackInt; n: LapackInt; v: ptr cdouble;
            tau: cdouble; c: ptr cdouble; ldc: LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dlarfx", dynlib: libName.}
proc dlarfxWork*(matrixLayout: cint; side: char; m: LapackInt; n: LapackInt;
                v: ptr cdouble; tau: cdouble; c: ptr cdouble; ldc: LapackInt;
                work: ptr cdouble): LapackInt {.importc: "LAPACKE_dlarfx_work",
    dynlib: libName.}
proc dlarnv*(idist: LapackInt; iseed: ptr LapackInt; n: LapackInt; x: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dlarnv", dynlib: libName.}
proc dlarnvWork*(idist: LapackInt; iseed: ptr LapackInt; n: LapackInt; x: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dlarnv_work", dynlib: libName.}
proc dlartgp*(f: cdouble; g: cdouble; cs: ptr cdouble; sn: ptr cdouble; r: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dlartgp", dynlib: libName.}
proc dlartgpWork*(f: cdouble; g: cdouble; cs: ptr cdouble; sn: ptr cdouble; r: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dlartgp_work", dynlib: libName.}
proc dlartgs*(x: cdouble; y: cdouble; sigma: cdouble; cs: ptr cdouble; sn: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dlartgs", dynlib: libName.}
proc dlartgsWork*(x: cdouble; y: cdouble; sigma: cdouble; cs: ptr cdouble; sn: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dlartgs_work", dynlib: libName.}
proc dlascl*(matrixLayout: cint; `type`: char; kl: LapackInt; ku: LapackInt;
            cfrom: cdouble; cto: cdouble; m: LapackInt; n: LapackInt; a: ptr cdouble;
            lda: LapackInt): LapackInt {.importc: "LAPACKE_dlascl", dynlib: libName.}
proc dlasclWork*(matrixLayout: cint; `type`: char; kl: LapackInt; ku: LapackInt;
                cfrom: cdouble; cto: cdouble; m: LapackInt; n: LapackInt;
                a: ptr cdouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_dlascl_work", dynlib: libName.}
proc dlaset*(matrixLayout: cint; uplo: char; m: LapackInt; n: LapackInt; alpha: cdouble;
            beta: cdouble; a: ptr cdouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_dlaset", dynlib: libName.}
proc dlasetWork*(matrixLayout: cint; uplo: char; m: LapackInt; n: LapackInt;
                alpha: cdouble; beta: cdouble; a: ptr cdouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_dlaset_work", dynlib: libName.}
proc dlasrt*(id: char; n: LapackInt; d: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dlasrt", dynlib: libName.}
proc dlasrtWork*(id: char; n: LapackInt; d: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dlasrt_work", dynlib: libName.}
proc dlaswp*(matrixLayout: cint; n: LapackInt; a: ptr cdouble; lda: LapackInt;
            k1: LapackInt; k2: LapackInt; ipiv: ptr LapackInt; incx: LapackInt): LapackInt {.
    importc: "LAPACKE_dlaswp", dynlib: libName.}
proc dlaswpWork*(matrixLayout: cint; n: LapackInt; a: ptr cdouble; lda: LapackInt;
                k1: LapackInt; k2: LapackInt; ipiv: ptr LapackInt; incx: LapackInt): LapackInt {.
    importc: "LAPACKE_dlaswp_work", dynlib: libName.}
proc dlatms*(matrixLayout: cint; m: LapackInt; n: LapackInt; dist: char;
            iseed: ptr LapackInt; sym: char; d: ptr cdouble; mode: LapackInt;
            cond: cdouble; dmax: cdouble; kl: LapackInt; ku: LapackInt; pack: char;
            a: ptr cdouble; lda: LapackInt): LapackInt {.importc: "LAPACKE_dlatms",
    dynlib: libName.}
proc dlatmsWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; dist: char;
                iseed: ptr LapackInt; sym: char; d: ptr cdouble; mode: LapackInt;
                cond: cdouble; dmax: cdouble; kl: LapackInt; ku: LapackInt; pack: char;
                a: ptr cdouble; lda: LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dlatms_work", dynlib: libName.}
proc dlauum*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_dlauum", dynlib: libName.}
proc dlauumWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                lda: LapackInt): LapackInt {.importc: "LAPACKE_dlauum_work",
    dynlib: libName.}
proc dopgtr*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cdouble;
            tau: ptr cdouble; q: ptr cdouble; ldq: LapackInt): LapackInt {.
    importc: "LAPACKE_dopgtr", dynlib: libName.}
proc dopgtrWork*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cdouble;
                tau: ptr cdouble; q: ptr cdouble; ldq: LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dopgtr_work", dynlib: libName.}
proc dopmtr*(matrixLayout: cint; side: char; uplo: char; trans: char; m: LapackInt;
            n: LapackInt; ap: ptr cdouble; tau: ptr cdouble; c: ptr cdouble; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_dopmtr", dynlib: libName.}
proc dopmtrWork*(matrixLayout: cint; side: char; uplo: char; trans: char; m: LapackInt;
                n: LapackInt; ap: ptr cdouble; tau: ptr cdouble; c: ptr cdouble;
                ldc: LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dopmtr_work", dynlib: libName.}
proc dorbdb*(matrixLayout: cint; trans: char; signs: char; m: LapackInt; p: LapackInt;
            q: LapackInt; x11: ptr cdouble; ldx11: LapackInt; x12: ptr cdouble;
            ldx12: LapackInt; x21: ptr cdouble; ldx21: LapackInt; x22: ptr cdouble;
            ldx22: LapackInt; theta: ptr cdouble; phi: ptr cdouble; taup1: ptr cdouble;
            taup2: ptr cdouble; tauq1: ptr cdouble; tauq2: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dorbdb", dynlib: libName.}
proc dorbdbWork*(matrixLayout: cint; trans: char; signs: char; m: LapackInt;
                p: LapackInt; q: LapackInt; x11: ptr cdouble; ldx11: LapackInt;
                x12: ptr cdouble; ldx12: LapackInt; x21: ptr cdouble; ldx21: LapackInt;
                x22: ptr cdouble; ldx22: LapackInt; theta: ptr cdouble;
                phi: ptr cdouble; taup1: ptr cdouble; taup2: ptr cdouble;
                tauq1: ptr cdouble; tauq2: ptr cdouble; work: ptr cdouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_dorbdb_work",
    dynlib: libName.}
proc dorcsd*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char; jobv2t: char;
            trans: char; signs: char; m: LapackInt; p: LapackInt; q: LapackInt;
            x11: ptr cdouble; ldx11: LapackInt; x12: ptr cdouble; ldx12: LapackInt;
            x21: ptr cdouble; ldx21: LapackInt; x22: ptr cdouble; ldx22: LapackInt;
            theta: ptr cdouble; u1: ptr cdouble; ldu1: LapackInt; u2: ptr cdouble;
            ldu2: LapackInt; v1t: ptr cdouble; ldv1t: LapackInt; v2t: ptr cdouble;
            ldv2t: LapackInt): LapackInt {.importc: "LAPACKE_dorcsd", dynlib: libName.}
proc dorcsd2by1*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char;
                m: LapackInt; p: LapackInt; q: LapackInt; x11: ptr cdouble;
                ldx11: LapackInt; x21: ptr cdouble; ldx21: LapackInt;
                theta: ptr cdouble; u1: ptr cdouble; ldu1: LapackInt; u2: ptr cdouble;
                ldu2: LapackInt; v1t: ptr cdouble; ldv1t: LapackInt): LapackInt {.
    importc: "LAPACKE_dorcsd2by1", dynlib: libName.}
proc dorcsd2by1Work*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char;
                    m: LapackInt; p: LapackInt; q: LapackInt; x11: ptr cdouble;
                    ldx11: LapackInt; x21: ptr cdouble; ldx21: LapackInt;
                    theta: ptr cdouble; u1: ptr cdouble; ldu1: LapackInt;
                    u2: ptr cdouble; ldu2: LapackInt; v1t: ptr cdouble;
                    ldv1t: LapackInt; work: ptr cdouble; lwork: LapackInt;
                    iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dorcsd2by1_work", dynlib: libName.}
proc dorcsdWork*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char;
                jobv2t: char; trans: char; signs: char; m: LapackInt; p: LapackInt;
                q: LapackInt; x11: ptr cdouble; ldx11: LapackInt; x12: ptr cdouble;
                ldx12: LapackInt; x21: ptr cdouble; ldx21: LapackInt; x22: ptr cdouble;
                ldx22: LapackInt; theta: ptr cdouble; u1: ptr cdouble; ldu1: LapackInt;
                u2: ptr cdouble; ldu2: LapackInt; v1t: ptr cdouble; ldv1t: LapackInt;
                v2t: ptr cdouble; ldv2t: LapackInt; work: ptr cdouble; lwork: LapackInt;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dorcsd_work",
    dynlib: libName.}
proc dorgbr*(matrixLayout: cint; vect: char; m: LapackInt; n: LapackInt; k: LapackInt;
            a: ptr cdouble; lda: LapackInt; tau: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dorgbr", dynlib: libName.}
proc dorgbrWork*(matrixLayout: cint; vect: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr cdouble; lda: LapackInt; tau: ptr cdouble;
                work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dorgbr_work", dynlib: libName.}
proc dorghr*(matrixLayout: cint; n: LapackInt; ilo: LapackInt; ihi: LapackInt;
            a: ptr cdouble; lda: LapackInt; tau: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dorghr", dynlib: libName.}
proc dorghrWork*(matrixLayout: cint; n: LapackInt; ilo: LapackInt; ihi: LapackInt;
                a: ptr cdouble; lda: LapackInt; tau: ptr cdouble; work: ptr cdouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_dorghr_work",
    dynlib: libName.}
proc dorglq*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
            a: ptr cdouble; lda: LapackInt; tau: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dorglq", dynlib: libName.}
proc dorglqWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
                a: ptr cdouble; lda: LapackInt; tau: ptr cdouble; work: ptr cdouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_dorglq_work",
    dynlib: libName.}
proc dorgql*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
            a: ptr cdouble; lda: LapackInt; tau: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dorgql", dynlib: libName.}
proc dorgqlWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
                a: ptr cdouble; lda: LapackInt; tau: ptr cdouble; work: ptr cdouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_dorgql_work",
    dynlib: libName.}
proc dorgqr*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
            a: ptr cdouble; lda: LapackInt; tau: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dorgqr", dynlib: libName.}
proc dorgqrWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
                a: ptr cdouble; lda: LapackInt; tau: ptr cdouble; work: ptr cdouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_dorgqr_work",
    dynlib: libName.}
proc dorgrq*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
            a: ptr cdouble; lda: LapackInt; tau: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dorgrq", dynlib: libName.}
proc dorgrqWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
                a: ptr cdouble; lda: LapackInt; tau: ptr cdouble; work: ptr cdouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_dorgrq_work",
    dynlib: libName.}
proc dorgtr*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble; lda: LapackInt;
            tau: ptr cdouble): LapackInt {.importc: "LAPACKE_dorgtr", dynlib: libName.}
proc dorgtrWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; tau: ptr cdouble; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dorgtr_work", dynlib: libName.}
proc dormbr*(matrixLayout: cint; vect: char; side: char; trans: char; m: LapackInt;
            n: LapackInt; k: LapackInt; a: ptr cdouble; lda: LapackInt; tau: ptr cdouble;
            c: ptr cdouble; ldc: LapackInt): LapackInt {.importc: "LAPACKE_dormbr",
    dynlib: libName.}
proc dormbrWork*(matrixLayout: cint; vect: char; side: char; trans: char; m: LapackInt;
                n: LapackInt; k: LapackInt; a: ptr cdouble; lda: LapackInt;
                tau: ptr cdouble; c: ptr cdouble; ldc: LapackInt; work: ptr cdouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_dormbr_work",
    dynlib: libName.}
proc dormhr*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            ilo: LapackInt; ihi: LapackInt; a: ptr cdouble; lda: LapackInt;
            tau: ptr cdouble; c: ptr cdouble; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_dormhr", dynlib: libName.}
proc dormhrWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                ilo: LapackInt; ihi: LapackInt; a: ptr cdouble; lda: LapackInt;
                tau: ptr cdouble; c: ptr cdouble; ldc: LapackInt; work: ptr cdouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_dormhr_work",
    dynlib: libName.}
proc dormlq*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr cdouble; lda: LapackInt; tau: ptr cdouble;
            c: ptr cdouble; ldc: LapackInt): LapackInt {.importc: "LAPACKE_dormlq",
    dynlib: libName.}
proc dormlqWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr cdouble; lda: LapackInt; tau: ptr cdouble;
                c: ptr cdouble; ldc: LapackInt; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dormlq_work", dynlib: libName.}
proc dormql*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr cdouble; lda: LapackInt; tau: ptr cdouble;
            c: ptr cdouble; ldc: LapackInt): LapackInt {.importc: "LAPACKE_dormql",
    dynlib: libName.}
proc dormqlWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr cdouble; lda: LapackInt; tau: ptr cdouble;
                c: ptr cdouble; ldc: LapackInt; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dormql_work", dynlib: libName.}
proc dormqr*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr cdouble; lda: LapackInt; tau: ptr cdouble;
            c: ptr cdouble; ldc: LapackInt): LapackInt {.importc: "LAPACKE_dormqr",
    dynlib: libName.}
proc dormqrWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr cdouble; lda: LapackInt; tau: ptr cdouble;
                c: ptr cdouble; ldc: LapackInt; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dormqr_work", dynlib: libName.}
proc dormrq*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr cdouble; lda: LapackInt; tau: ptr cdouble;
            c: ptr cdouble; ldc: LapackInt): LapackInt {.importc: "LAPACKE_dormrq",
    dynlib: libName.}
proc dormrqWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr cdouble; lda: LapackInt; tau: ptr cdouble;
                c: ptr cdouble; ldc: LapackInt; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dormrq_work", dynlib: libName.}
proc dormrz*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; l: LapackInt; a: ptr cdouble; lda: LapackInt; tau: ptr cdouble;
            c: ptr cdouble; ldc: LapackInt): LapackInt {.importc: "LAPACKE_dormrz",
    dynlib: libName.}
proc dormrzWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; l: LapackInt; a: ptr cdouble; lda: LapackInt;
                tau: ptr cdouble; c: ptr cdouble; ldc: LapackInt; work: ptr cdouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_dormrz_work",
    dynlib: libName.}
proc dormtr*(matrixLayout: cint; side: char; uplo: char; trans: char; m: LapackInt;
            n: LapackInt; a: ptr cdouble; lda: LapackInt; tau: ptr cdouble;
            c: ptr cdouble; ldc: LapackInt): LapackInt {.importc: "LAPACKE_dormtr",
    dynlib: libName.}
proc dormtrWork*(matrixLayout: cint; side: char; uplo: char; trans: char; m: LapackInt;
                n: LapackInt; a: ptr cdouble; lda: LapackInt; tau: ptr cdouble;
                c: ptr cdouble; ldc: LapackInt; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dormtr_work", dynlib: libName.}
proc dpbcon*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
            ab: ptr cdouble; ldab: LapackInt; anorm: cdouble; rcond: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dpbcon", dynlib: libName.}
proc dpbconWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
                ab: ptr cdouble; ldab: LapackInt; anorm: cdouble; rcond: ptr cdouble;
                work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dpbcon_work", dynlib: libName.}
proc dpbequ*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
            ab: ptr cdouble; ldab: LapackInt; s: ptr cdouble; scond: ptr cdouble;
            amax: ptr cdouble): LapackInt {.importc: "LAPACKE_dpbequ", dynlib: libName.}
proc dpbequWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
                ab: ptr cdouble; ldab: LapackInt; s: ptr cdouble; scond: ptr cdouble;
                amax: ptr cdouble): LapackInt {.importc: "LAPACKE_dpbequ_work",
    dynlib: libName.}
proc dpbrfs*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
            nrhs: LapackInt; ab: ptr cdouble; ldab: LapackInt; afb: ptr cdouble;
            ldafb: LapackInt; b: ptr cdouble; ldb: LapackInt; x: ptr cdouble;
            ldx: LapackInt; ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dpbrfs", dynlib: libName.}
proc dpbrfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
                nrhs: LapackInt; ab: ptr cdouble; ldab: LapackInt; afb: ptr cdouble;
                ldafb: LapackInt; b: ptr cdouble; ldb: LapackInt; x: ptr cdouble;
                ldx: LapackInt; ferr: ptr cdouble; berr: ptr cdouble; work: ptr cdouble;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dpbrfs_work",
    dynlib: libName.}
proc dpbstf*(matrixLayout: cint; uplo: char; n: LapackInt; kb: LapackInt;
            bb: ptr cdouble; ldbb: LapackInt): LapackInt {.importc: "LAPACKE_dpbstf",
    dynlib: libName.}
proc dpbstfWork*(matrixLayout: cint; uplo: char; n: LapackInt; kb: LapackInt;
                bb: ptr cdouble; ldbb: LapackInt): LapackInt {.
    importc: "LAPACKE_dpbstf_work", dynlib: libName.}
proc dpbsv*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt; nrhs: LapackInt;
           ab: ptr cdouble; ldab: LapackInt; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dpbsv", dynlib: libName.}
proc dpbsvWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
               nrhs: LapackInt; ab: ptr cdouble; ldab: LapackInt; b: ptr cdouble;
               ldb: LapackInt): LapackInt {.importc: "LAPACKE_dpbsv_work",
    dynlib: libName.}
proc dpbsvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; kd: LapackInt;
            nrhs: LapackInt; ab: ptr cdouble; ldab: LapackInt; afb: ptr cdouble;
            ldafb: LapackInt; equed: cstring; s: ptr cdouble; b: ptr cdouble;
            ldb: LapackInt; x: ptr cdouble; ldx: LapackInt; rcond: ptr cdouble;
            ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dpbsvx", dynlib: libName.}
proc dpbsvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; kd: LapackInt;
                nrhs: LapackInt; ab: ptr cdouble; ldab: LapackInt; afb: ptr cdouble;
                ldafb: LapackInt; equed: cstring; s: ptr cdouble; b: ptr cdouble;
                ldb: LapackInt; x: ptr cdouble; ldx: LapackInt; rcond: ptr cdouble;
                ferr: ptr cdouble; berr: ptr cdouble; work: ptr cdouble;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dpbsvx_work",
    dynlib: libName.}
proc dpbtrf*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
            ab: ptr cdouble; ldab: LapackInt): LapackInt {.importc: "LAPACKE_dpbtrf",
    dynlib: libName.}
proc dpbtrfWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
                ab: ptr cdouble; ldab: LapackInt): LapackInt {.
    importc: "LAPACKE_dpbtrf_work", dynlib: libName.}
proc dpbtrs*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
            nrhs: LapackInt; ab: ptr cdouble; ldab: LapackInt; b: ptr cdouble;
            ldb: LapackInt): LapackInt {.importc: "LAPACKE_dpbtrs", dynlib: libName.}
proc dpbtrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
                nrhs: LapackInt; ab: ptr cdouble; ldab: LapackInt; b: ptr cdouble;
                ldb: LapackInt): LapackInt {.importc: "LAPACKE_dpbtrs_work",
    dynlib: libName.}
proc dpftrf*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt; a: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dpftrf", dynlib: libName.}
proc dpftrfWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                a: ptr cdouble): LapackInt {.importc: "LAPACKE_dpftrf_work",
    dynlib: libName.}
proc dpftri*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt; a: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dpftri", dynlib: libName.}
proc dpftriWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                a: ptr cdouble): LapackInt {.importc: "LAPACKE_dpftri_work",
    dynlib: libName.}
proc dpftrs*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr cdouble; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dpftrs", dynlib: libName.}
proc dpftrsWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; a: ptr cdouble; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dpftrs_work", dynlib: libName.}
proc dpocon*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble; lda: LapackInt;
            anorm: cdouble; rcond: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dpocon", dynlib: libName.}
proc dpoconWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; anorm: cdouble; rcond: ptr cdouble; work: ptr cdouble;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dpocon_work",
    dynlib: libName.}
proc dpoequ*(matrixLayout: cint; n: LapackInt; a: ptr cdouble; lda: LapackInt;
            s: ptr cdouble; scond: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dpoequ", dynlib: libName.}
proc dpoequWork*(matrixLayout: cint; n: LapackInt; a: ptr cdouble; lda: LapackInt;
                s: ptr cdouble; scond: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dpoequ_work", dynlib: libName.}
proc dpoequb*(matrixLayout: cint; n: LapackInt; a: ptr cdouble; lda: LapackInt;
             s: ptr cdouble; scond: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dpoequb", dynlib: libName.}
proc dpoequbWork*(matrixLayout: cint; n: LapackInt; a: ptr cdouble; lda: LapackInt;
                 s: ptr cdouble; scond: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dpoequb_work", dynlib: libName.}
proc dporfs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr cdouble; lda: LapackInt; af: ptr cdouble; ldaf: LapackInt;
            b: ptr cdouble; ldb: LapackInt; x: ptr cdouble; ldx: LapackInt;
            ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dporfs", dynlib: libName.}
proc dporfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr cdouble; lda: LapackInt; af: ptr cdouble; ldaf: LapackInt;
                b: ptr cdouble; ldb: LapackInt; x: ptr cdouble; ldx: LapackInt;
                ferr: ptr cdouble; berr: ptr cdouble; work: ptr cdouble;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dporfs_work",
    dynlib: libName.}
proc dporfsx*(matrixLayout: cint; uplo: char; equed: char; n: LapackInt; nrhs: LapackInt;
             a: ptr cdouble; lda: LapackInt; af: ptr cdouble; ldaf: LapackInt;
             s: ptr cdouble; b: ptr cdouble; ldb: LapackInt; x: ptr cdouble;
             ldx: LapackInt; rcond: ptr cdouble; berr: ptr cdouble; nErrBnds: LapackInt;
             errBndsNorm: ptr cdouble; errBndsComp: ptr cdouble; nparams: LapackInt;
             params: ptr cdouble): LapackInt {.importc: "LAPACKE_dporfsx",
    dynlib: libName.}
proc dporfsxWork*(matrixLayout: cint; uplo: char; equed: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr cdouble; lda: LapackInt; af: ptr cdouble;
                 ldaf: LapackInt; s: ptr cdouble; b: ptr cdouble; ldb: LapackInt;
                 x: ptr cdouble; ldx: LapackInt; rcond: ptr cdouble; berr: ptr cdouble;
                 nErrBnds: LapackInt; errBndsNorm: ptr cdouble;
                 errBndsComp: ptr cdouble; nparams: LapackInt; params: ptr cdouble;
                 work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dporfsx_work", dynlib: libName.}
proc dposv*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
           a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dposv", dynlib: libName.}
proc dposvWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dposv_work", dynlib: libName.}
proc dposvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr cdouble; lda: LapackInt; af: ptr cdouble; ldaf: LapackInt;
            equed: cstring; s: ptr cdouble; b: ptr cdouble; ldb: LapackInt;
            x: ptr cdouble; ldx: LapackInt; rcond: ptr cdouble; ferr: ptr cdouble;
            berr: ptr cdouble): LapackInt {.importc: "LAPACKE_dposvx", dynlib: libName.}
proc dposvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; a: ptr cdouble; lda: LapackInt; af: ptr cdouble;
                ldaf: LapackInt; equed: cstring; s: ptr cdouble; b: ptr cdouble;
                ldb: LapackInt; x: ptr cdouble; ldx: LapackInt; rcond: ptr cdouble;
                ferr: ptr cdouble; berr: ptr cdouble; work: ptr cdouble;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dposvx_work",
    dynlib: libName.}
proc dposvxx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr cdouble; lda: LapackInt; af: ptr cdouble; ldaf: LapackInt;
             equed: cstring; s: ptr cdouble; b: ptr cdouble; ldb: LapackInt;
             x: ptr cdouble; ldx: LapackInt; rcond: ptr cdouble; rpvgrw: ptr cdouble;
             berr: ptr cdouble; nErrBnds: LapackInt; errBndsNorm: ptr cdouble;
             errBndsComp: ptr cdouble; nparams: LapackInt; params: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dposvxx", dynlib: libName.}
proc dposvxxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr cdouble; lda: LapackInt; af: ptr cdouble;
                 ldaf: LapackInt; equed: cstring; s: ptr cdouble; b: ptr cdouble;
                 ldb: LapackInt; x: ptr cdouble; ldx: LapackInt; rcond: ptr cdouble;
                 rpvgrw: ptr cdouble; berr: ptr cdouble; nErrBnds: LapackInt;
                 errBndsNorm: ptr cdouble; errBndsComp: ptr cdouble;
                 nparams: LapackInt; params: ptr cdouble; work: ptr cdouble;
                 iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dposvxx_work",
    dynlib: libName.}
proc dpotrf*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_dpotrf", dynlib: libName.}
proc dpotrf2*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_dpotrf2", dynlib: libName.}
proc dpotrf2Work*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                 lda: LapackInt): LapackInt {.importc: "LAPACKE_dpotrf2_work",
    dynlib: libName.}
proc dpotrfWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                lda: LapackInt): LapackInt {.importc: "LAPACKE_dpotrf_work",
    dynlib: libName.}
proc dpotri*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_dpotri", dynlib: libName.}
proc dpotriWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                lda: LapackInt): LapackInt {.importc: "LAPACKE_dpotri_work",
    dynlib: libName.}
proc dpotrs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dpotrs", dynlib: libName.}
proc dpotrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dpotrs_work", dynlib: libName.}
proc dppcon*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cdouble;
            anorm: cdouble; rcond: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dppcon", dynlib: libName.}
proc dppconWork*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cdouble;
                anorm: cdouble; rcond: ptr cdouble; work: ptr cdouble;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dppcon_work",
    dynlib: libName.}
proc dppequ*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cdouble;
            s: ptr cdouble; scond: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dppequ", dynlib: libName.}
proc dppequWork*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cdouble;
                s: ptr cdouble; scond: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dppequ_work", dynlib: libName.}
proc dpprfs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr cdouble; afp: ptr cdouble; b: ptr cdouble; ldb: LapackInt;
            x: ptr cdouble; ldx: LapackInt; ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dpprfs", dynlib: libName.}
proc dpprfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                ap: ptr cdouble; afp: ptr cdouble; b: ptr cdouble; ldb: LapackInt;
                x: ptr cdouble; ldx: LapackInt; ferr: ptr cdouble; berr: ptr cdouble;
                work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dpprfs_work", dynlib: libName.}
proc dppsv*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
           ap: ptr cdouble; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dppsv", dynlib: libName.}
proc dppsvWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               ap: ptr cdouble; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dppsv_work", dynlib: libName.}
proc dppsvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr cdouble; afp: ptr cdouble; equed: cstring; s: ptr cdouble;
            b: ptr cdouble; ldb: LapackInt; x: ptr cdouble; ldx: LapackInt;
            rcond: ptr cdouble; ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dppsvx", dynlib: libName.}
proc dppsvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; ap: ptr cdouble; afp: ptr cdouble; equed: cstring;
                s: ptr cdouble; b: ptr cdouble; ldb: LapackInt; x: ptr cdouble;
                ldx: LapackInt; rcond: ptr cdouble; ferr: ptr cdouble;
                berr: ptr cdouble; work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dppsvx_work", dynlib: libName.}
proc dpptrf*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dpptrf", dynlib: libName.}
proc dpptrfWork*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dpptrf_work", dynlib: libName.}
proc dpptri*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dpptri", dynlib: libName.}
proc dpptriWork*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dpptri_work", dynlib: libName.}
proc dpptrs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr cdouble; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dpptrs", dynlib: libName.}
proc dpptrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                ap: ptr cdouble; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dpptrs_work", dynlib: libName.}
proc dpstrf*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble; lda: LapackInt;
            piv: ptr LapackInt; rank: ptr LapackInt; tol: cdouble): LapackInt {.
    importc: "LAPACKE_dpstrf", dynlib: libName.}
proc dpstrfWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; piv: ptr LapackInt; rank: ptr LapackInt; tol: cdouble;
                work: ptr cdouble): LapackInt {.importc: "LAPACKE_dpstrf_work",
    dynlib: libName.}
proc dptcon*(n: LapackInt; d: ptr cdouble; e: ptr cdouble; anorm: cdouble;
            rcond: ptr cdouble): LapackInt {.importc: "LAPACKE_dptcon",
    dynlib: libName.}
proc dptconWork*(n: LapackInt; d: ptr cdouble; e: ptr cdouble; anorm: cdouble;
                rcond: ptr cdouble; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dptcon_work", dynlib: libName.}
proc dpteqr*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cdouble;
            e: ptr cdouble; z: ptr cdouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_dpteqr", dynlib: libName.}
proc dpteqrWork*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cdouble;
                e: ptr cdouble; z: ptr cdouble; ldz: LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dpteqr_work", dynlib: libName.}
proc dptrfs*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; d: ptr cdouble;
            e: ptr cdouble; df: ptr cdouble; ef: ptr cdouble; b: ptr cdouble;
            ldb: LapackInt; x: ptr cdouble; ldx: LapackInt; ferr: ptr cdouble;
            berr: ptr cdouble): LapackInt {.importc: "LAPACKE_dptrfs", dynlib: libName.}
proc dptrfsWork*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; d: ptr cdouble;
                e: ptr cdouble; df: ptr cdouble; ef: ptr cdouble; b: ptr cdouble;
                ldb: LapackInt; x: ptr cdouble; ldx: LapackInt; ferr: ptr cdouble;
                berr: ptr cdouble; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dptrfs_work", dynlib: libName.}
proc dptsv*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; d: ptr cdouble;
           e: ptr cdouble; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dptsv", dynlib: libName.}
proc dptsvWork*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; d: ptr cdouble;
               e: ptr cdouble; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dptsv_work", dynlib: libName.}
proc dptsvx*(matrixLayout: cint; fact: char; n: LapackInt; nrhs: LapackInt;
            d: ptr cdouble; e: ptr cdouble; df: ptr cdouble; ef: ptr cdouble;
            b: ptr cdouble; ldb: LapackInt; x: ptr cdouble; ldx: LapackInt;
            rcond: ptr cdouble; ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dptsvx", dynlib: libName.}
proc dptsvxWork*(matrixLayout: cint; fact: char; n: LapackInt; nrhs: LapackInt;
                d: ptr cdouble; e: ptr cdouble; df: ptr cdouble; ef: ptr cdouble;
                b: ptr cdouble; ldb: LapackInt; x: ptr cdouble; ldx: LapackInt;
                rcond: ptr cdouble; ferr: ptr cdouble; berr: ptr cdouble;
                work: ptr cdouble): LapackInt {.importc: "LAPACKE_dptsvx_work",
    dynlib: libName.}
proc dpttrf*(n: LapackInt; d: ptr cdouble; e: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dpttrf", dynlib: libName.}
proc dpttrfWork*(n: LapackInt; d: ptr cdouble; e: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dpttrf_work", dynlib: libName.}
proc dpttrs*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; d: ptr cdouble;
            e: ptr cdouble; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dpttrs", dynlib: libName.}
proc dpttrsWork*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; d: ptr cdouble;
                e: ptr cdouble; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dpttrs_work", dynlib: libName.}
proc dsbev*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; kd: LapackInt;
           ab: ptr cdouble; ldab: LapackInt; w: ptr cdouble; z: ptr cdouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_dsbev", dynlib: libName.}
proc dsbevWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; kd: LapackInt;
               ab: ptr cdouble; ldab: LapackInt; w: ptr cdouble; z: ptr cdouble;
               ldz: LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsbev_work", dynlib: libName.}
proc dsbevd*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; kd: LapackInt;
            ab: ptr cdouble; ldab: LapackInt; w: ptr cdouble; z: ptr cdouble;
            ldz: LapackInt): LapackInt {.importc: "LAPACKE_dsbevd", dynlib: libName.}
proc dsbevdWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; kd: LapackInt;
                ab: ptr cdouble; ldab: LapackInt; w: ptr cdouble; z: ptr cdouble;
                ldz: LapackInt; work: ptr cdouble; lwork: LapackInt;
                iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dsbevd_work", dynlib: libName.}
proc dsbevx*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
            kd: LapackInt; ab: ptr cdouble; ldab: LapackInt; q: ptr cdouble;
            ldq: LapackInt; vl: cdouble; vu: cdouble; il: LapackInt; iu: LapackInt;
            abstol: cdouble; m: ptr LapackInt; w: ptr cdouble; z: ptr cdouble;
            ldz: LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsbevx", dynlib: libName.}
proc dsbevxWork*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                kd: LapackInt; ab: ptr cdouble; ldab: LapackInt; q: ptr cdouble;
                ldq: LapackInt; vl: cdouble; vu: cdouble; il: LapackInt; iu: LapackInt;
                abstol: cdouble; m: ptr LapackInt; w: ptr cdouble; z: ptr cdouble;
                ldz: LapackInt; work: ptr cdouble; iwork: ptr LapackInt;
                ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_dsbevx_work",
    dynlib: libName.}
proc dsbgst*(matrixLayout: cint; vect: char; uplo: char; n: LapackInt; ka: LapackInt;
            kb: LapackInt; ab: ptr cdouble; ldab: LapackInt; bb: ptr cdouble;
            ldbb: LapackInt; x: ptr cdouble; ldx: LapackInt): LapackInt {.
    importc: "LAPACKE_dsbgst", dynlib: libName.}
proc dsbgstWork*(matrixLayout: cint; vect: char; uplo: char; n: LapackInt; ka: LapackInt;
                kb: LapackInt; ab: ptr cdouble; ldab: LapackInt; bb: ptr cdouble;
                ldbb: LapackInt; x: ptr cdouble; ldx: LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsbgst_work", dynlib: libName.}
proc dsbgv*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; ka: LapackInt;
           kb: LapackInt; ab: ptr cdouble; ldab: LapackInt; bb: ptr cdouble;
           ldbb: LapackInt; w: ptr cdouble; z: ptr cdouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_dsbgv", dynlib: libName.}
proc dsbgvWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; ka: LapackInt;
               kb: LapackInt; ab: ptr cdouble; ldab: LapackInt; bb: ptr cdouble;
               ldbb: LapackInt; w: ptr cdouble; z: ptr cdouble; ldz: LapackInt;
               work: ptr cdouble): LapackInt {.importc: "LAPACKE_dsbgv_work",
    dynlib: libName.}
proc dsbgvd*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; ka: LapackInt;
            kb: LapackInt; ab: ptr cdouble; ldab: LapackInt; bb: ptr cdouble;
            ldbb: LapackInt; w: ptr cdouble; z: ptr cdouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_dsbgvd", dynlib: libName.}
proc dsbgvdWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; ka: LapackInt;
                kb: LapackInt; ab: ptr cdouble; ldab: LapackInt; bb: ptr cdouble;
                ldbb: LapackInt; w: ptr cdouble; z: ptr cdouble; ldz: LapackInt;
                work: ptr cdouble; lwork: LapackInt; iwork: ptr LapackInt;
                liwork: LapackInt): LapackInt {.importc: "LAPACKE_dsbgvd_work",
    dynlib: libName.}
proc dsbgvx*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
            ka: LapackInt; kb: LapackInt; ab: ptr cdouble; ldab: LapackInt;
            bb: ptr cdouble; ldbb: LapackInt; q: ptr cdouble; ldq: LapackInt; vl: cdouble;
            vu: cdouble; il: LapackInt; iu: LapackInt; abstol: cdouble;
            m: ptr LapackInt; w: ptr cdouble; z: ptr cdouble; ldz: LapackInt;
            ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_dsbgvx",
    dynlib: libName.}
proc dsbgvxWork*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                ka: LapackInt; kb: LapackInt; ab: ptr cdouble; ldab: LapackInt;
                bb: ptr cdouble; ldbb: LapackInt; q: ptr cdouble; ldq: LapackInt;
                vl: cdouble; vu: cdouble; il: LapackInt; iu: LapackInt; abstol: cdouble;
                m: ptr LapackInt; w: ptr cdouble; z: ptr cdouble; ldz: LapackInt;
                work: ptr cdouble; iwork: ptr LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsbgvx_work", dynlib: libName.}
proc dsbtrd*(matrixLayout: cint; vect: char; uplo: char; n: LapackInt; kd: LapackInt;
            ab: ptr cdouble; ldab: LapackInt; d: ptr cdouble; e: ptr cdouble;
            q: ptr cdouble; ldq: LapackInt): LapackInt {.importc: "LAPACKE_dsbtrd",
    dynlib: libName.}
proc dsbtrdWork*(matrixLayout: cint; vect: char; uplo: char; n: LapackInt; kd: LapackInt;
                ab: ptr cdouble; ldab: LapackInt; d: ptr cdouble; e: ptr cdouble;
                q: ptr cdouble; ldq: LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsbtrd_work", dynlib: libName.}
proc dsfrk*(matrixLayout: cint; transr: char; uplo: char; trans: char; n: LapackInt;
           k: LapackInt; alpha: cdouble; a: ptr cdouble; lda: LapackInt; beta: cdouble;
           c: ptr cdouble): LapackInt {.importc: "LAPACKE_dsfrk", dynlib: libName.}
proc dsfrkWork*(matrixLayout: cint; transr: char; uplo: char; trans: char; n: LapackInt;
               k: LapackInt; alpha: cdouble; a: ptr cdouble; lda: LapackInt;
               beta: cdouble; c: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsfrk_work", dynlib: libName.}
proc dsgesv*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; a: ptr cdouble;
            lda: LapackInt; ipiv: ptr LapackInt; b: ptr cdouble; ldb: LapackInt;
            x: ptr cdouble; ldx: LapackInt; iter: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsgesv", dynlib: libName.}
proc dsgesvWork*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; a: ptr cdouble;
                lda: LapackInt; ipiv: ptr LapackInt; b: ptr cdouble; ldb: LapackInt;
                x: ptr cdouble; ldx: LapackInt; work: ptr cdouble; swork: ptr cfloat;
                iter: ptr LapackInt): LapackInt {.importc: "LAPACKE_dsgesv_work",
    dynlib: libName.}
proc dspcon*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cdouble;
            ipiv: ptr LapackInt; anorm: cdouble; rcond: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dspcon", dynlib: libName.}
proc dspconWork*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cdouble;
                ipiv: ptr LapackInt; anorm: cdouble; rcond: ptr cdouble;
                work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dspcon_work", dynlib: libName.}
proc dspev*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; ap: ptr cdouble;
           w: ptr cdouble; z: ptr cdouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_dspev", dynlib: libName.}
proc dspevWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
               ap: ptr cdouble; w: ptr cdouble; z: ptr cdouble; ldz: LapackInt;
               work: ptr cdouble): LapackInt {.importc: "LAPACKE_dspev_work",
    dynlib: libName.}
proc dspevd*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; ap: ptr cdouble;
            w: ptr cdouble; z: ptr cdouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_dspevd", dynlib: libName.}
proc dspevdWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                ap: ptr cdouble; w: ptr cdouble; z: ptr cdouble; ldz: LapackInt;
                work: ptr cdouble; lwork: LapackInt; iwork: ptr LapackInt;
                liwork: LapackInt): LapackInt {.importc: "LAPACKE_dspevd_work",
    dynlib: libName.}
proc dspevx*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
            ap: ptr cdouble; vl: cdouble; vu: cdouble; il: LapackInt; iu: LapackInt;
            abstol: cdouble; m: ptr LapackInt; w: ptr cdouble; z: ptr cdouble;
            ldz: LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dspevx", dynlib: libName.}
proc dspevxWork*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                ap: ptr cdouble; vl: cdouble; vu: cdouble; il: LapackInt; iu: LapackInt;
                abstol: cdouble; m: ptr LapackInt; w: ptr cdouble; z: ptr cdouble;
                ldz: LapackInt; work: ptr cdouble; iwork: ptr LapackInt;
                ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_dspevx_work",
    dynlib: libName.}
proc dspgst*(matrixLayout: cint; itype: LapackInt; uplo: char; n: LapackInt;
            ap: ptr cdouble; bp: ptr cdouble): LapackInt {.importc: "LAPACKE_dspgst",
    dynlib: libName.}
proc dspgstWork*(matrixLayout: cint; itype: LapackInt; uplo: char; n: LapackInt;
                ap: ptr cdouble; bp: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dspgst_work", dynlib: libName.}
proc dspgv*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char; n: LapackInt;
           ap: ptr cdouble; bp: ptr cdouble; w: ptr cdouble; z: ptr cdouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_dspgv", dynlib: libName.}
proc dspgvWork*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
               n: LapackInt; ap: ptr cdouble; bp: ptr cdouble; w: ptr cdouble;
               z: ptr cdouble; ldz: LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dspgv_work", dynlib: libName.}
proc dspgvd*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char; n: LapackInt;
            ap: ptr cdouble; bp: ptr cdouble; w: ptr cdouble; z: ptr cdouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_dspgvd", dynlib: libName.}
proc dspgvdWork*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
                n: LapackInt; ap: ptr cdouble; bp: ptr cdouble; w: ptr cdouble;
                z: ptr cdouble; ldz: LapackInt; work: ptr cdouble; lwork: LapackInt;
                iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dspgvd_work", dynlib: libName.}
proc dspgvx*(matrixLayout: cint; itype: LapackInt; jobz: char; range: char; uplo: char;
            n: LapackInt; ap: ptr cdouble; bp: ptr cdouble; vl: cdouble; vu: cdouble;
            il: LapackInt; iu: LapackInt; abstol: cdouble; m: ptr LapackInt;
            w: ptr cdouble; z: ptr cdouble; ldz: LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dspgvx", dynlib: libName.}
proc dspgvxWork*(matrixLayout: cint; itype: LapackInt; jobz: char; range: char;
                uplo: char; n: LapackInt; ap: ptr cdouble; bp: ptr cdouble; vl: cdouble;
                vu: cdouble; il: LapackInt; iu: LapackInt; abstol: cdouble;
                m: ptr LapackInt; w: ptr cdouble; z: ptr cdouble; ldz: LapackInt;
                work: ptr cdouble; iwork: ptr LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dspgvx_work", dynlib: libName.}
proc dsposv*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
            x: ptr cdouble; ldx: LapackInt; iter: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsposv", dynlib: libName.}
proc dsposvWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
                x: ptr cdouble; ldx: LapackInt; work: ptr cdouble; swork: ptr cfloat;
                iter: ptr LapackInt): LapackInt {.importc: "LAPACKE_dsposv_work",
    dynlib: libName.}
proc dsprfs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr cdouble; afp: ptr cdouble; ipiv: ptr LapackInt; b: ptr cdouble;
            ldb: LapackInt; x: ptr cdouble; ldx: LapackInt; ferr: ptr cdouble;
            berr: ptr cdouble): LapackInt {.importc: "LAPACKE_dsprfs", dynlib: libName.}
proc dsprfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                ap: ptr cdouble; afp: ptr cdouble; ipiv: ptr LapackInt; b: ptr cdouble;
                ldb: LapackInt; x: ptr cdouble; ldx: LapackInt; ferr: ptr cdouble;
                berr: ptr cdouble; work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsprfs_work", dynlib: libName.}
proc dspsv*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
           ap: ptr cdouble; ipiv: ptr LapackInt; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dspsv", dynlib: libName.}
proc dspsvWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               ap: ptr cdouble; ipiv: ptr LapackInt; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dspsv_work", dynlib: libName.}
proc dspsvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr cdouble; afp: ptr cdouble; ipiv: ptr LapackInt; b: ptr cdouble;
            ldb: LapackInt; x: ptr cdouble; ldx: LapackInt; rcond: ptr cdouble;
            ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dspsvx", dynlib: libName.}
proc dspsvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; ap: ptr cdouble; afp: ptr cdouble; ipiv: ptr LapackInt;
                b: ptr cdouble; ldb: LapackInt; x: ptr cdouble; ldx: LapackInt;
                rcond: ptr cdouble; ferr: ptr cdouble; berr: ptr cdouble;
                work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dspsvx_work", dynlib: libName.}
proc dsptrd*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cdouble;
            d: ptr cdouble; e: ptr cdouble; tau: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsptrd", dynlib: libName.}
proc dsptrdWork*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cdouble;
                d: ptr cdouble; e: ptr cdouble; tau: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsptrd_work", dynlib: libName.}
proc dsptrf*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cdouble;
            ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_dsptrf",
    dynlib: libName.}
proc dsptrfWork*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cdouble;
                ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_dsptrf_work",
    dynlib: libName.}
proc dsptri*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cdouble;
            ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_dsptri",
    dynlib: libName.}
proc dsptriWork*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cdouble;
                ipiv: ptr LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsptri_work", dynlib: libName.}
proc dsptrs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr cdouble; ipiv: ptr LapackInt; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dsptrs", dynlib: libName.}
proc dsptrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                ap: ptr cdouble; ipiv: ptr LapackInt; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dsptrs_work", dynlib: libName.}
proc dstebz*(range: char; order: char; n: LapackInt; vl: cdouble; vu: cdouble;
            il: LapackInt; iu: LapackInt; abstol: cdouble; d: ptr cdouble;
            e: ptr cdouble; m: ptr LapackInt; nsplit: ptr LapackInt; w: ptr cdouble;
            iblock: ptr LapackInt; isplit: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dstebz", dynlib: libName.}
proc dstebzWork*(range: char; order: char; n: LapackInt; vl: cdouble; vu: cdouble;
                il: LapackInt; iu: LapackInt; abstol: cdouble; d: ptr cdouble;
                e: ptr cdouble; m: ptr LapackInt; nsplit: ptr LapackInt; w: ptr cdouble;
                iblock: ptr LapackInt; isplit: ptr LapackInt; work: ptr cdouble;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dstebz_work",
    dynlib: libName.}
proc dstedc*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cdouble;
            e: ptr cdouble; z: ptr cdouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_dstedc", dynlib: libName.}
proc dstedcWork*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cdouble;
                e: ptr cdouble; z: ptr cdouble; ldz: LapackInt; work: ptr cdouble;
                lwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dstedc_work", dynlib: libName.}
proc dstegr*(matrixLayout: cint; jobz: char; range: char; n: LapackInt; d: ptr cdouble;
            e: ptr cdouble; vl: cdouble; vu: cdouble; il: LapackInt; iu: LapackInt;
            abstol: cdouble; m: ptr LapackInt; w: ptr cdouble; z: ptr cdouble;
            ldz: LapackInt; isuppz: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dstegr", dynlib: libName.}
proc dstegrWork*(matrixLayout: cint; jobz: char; range: char; n: LapackInt;
                d: ptr cdouble; e: ptr cdouble; vl: cdouble; vu: cdouble; il: LapackInt;
                iu: LapackInt; abstol: cdouble; m: ptr LapackInt; w: ptr cdouble;
                z: ptr cdouble; ldz: LapackInt; isuppz: ptr LapackInt;
                work: ptr cdouble; lwork: LapackInt; iwork: ptr LapackInt;
                liwork: LapackInt): LapackInt {.importc: "LAPACKE_dstegr_work",
    dynlib: libName.}
proc dstein*(matrixLayout: cint; n: LapackInt; d: ptr cdouble; e: ptr cdouble;
            m: LapackInt; w: ptr cdouble; iblock: ptr LapackInt; isplit: ptr LapackInt;
            z: ptr cdouble; ldz: LapackInt; ifailv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dstein", dynlib: libName.}
proc dsteinWork*(matrixLayout: cint; n: LapackInt; d: ptr cdouble; e: ptr cdouble;
                m: LapackInt; w: ptr cdouble; iblock: ptr LapackInt;
                isplit: ptr LapackInt; z: ptr cdouble; ldz: LapackInt;
                work: ptr cdouble; iwork: ptr LapackInt; ifailv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dstein_work", dynlib: libName.}
proc dstemr*(matrixLayout: cint; jobz: char; range: char; n: LapackInt; d: ptr cdouble;
            e: ptr cdouble; vl: cdouble; vu: cdouble; il: LapackInt; iu: LapackInt;
            m: ptr LapackInt; w: ptr cdouble; z: ptr cdouble; ldz: LapackInt;
            nzc: LapackInt; isuppz: ptr LapackInt; tryrac: ptr LapackLogical): LapackInt {.
    importc: "LAPACKE_dstemr", dynlib: libName.}
proc dstemrWork*(matrixLayout: cint; jobz: char; range: char; n: LapackInt;
                d: ptr cdouble; e: ptr cdouble; vl: cdouble; vu: cdouble; il: LapackInt;
                iu: LapackInt; m: ptr LapackInt; w: ptr cdouble; z: ptr cdouble;
                ldz: LapackInt; nzc: LapackInt; isuppz: ptr LapackInt;
                tryrac: ptr LapackLogical; work: ptr cdouble; lwork: LapackInt;
                iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dstemr_work", dynlib: libName.}
proc dsteqr*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cdouble;
            e: ptr cdouble; z: ptr cdouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_dsteqr", dynlib: libName.}
proc dsteqrWork*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cdouble;
                e: ptr cdouble; z: ptr cdouble; ldz: LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsteqr_work", dynlib: libName.}
proc dsterf*(n: LapackInt; d: ptr cdouble; e: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsterf", dynlib: libName.}
proc dsterfWork*(n: LapackInt; d: ptr cdouble; e: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsterf_work", dynlib: libName.}
proc dstev*(matrixLayout: cint; jobz: char; n: LapackInt; d: ptr cdouble; e: ptr cdouble;
           z: ptr cdouble; ldz: LapackInt): LapackInt {.importc: "LAPACKE_dstev",
    dynlib: libName.}
proc dstevWork*(matrixLayout: cint; jobz: char; n: LapackInt; d: ptr cdouble;
               e: ptr cdouble; z: ptr cdouble; ldz: LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dstev_work", dynlib: libName.}
proc dstevd*(matrixLayout: cint; jobz: char; n: LapackInt; d: ptr cdouble; e: ptr cdouble;
            z: ptr cdouble; ldz: LapackInt): LapackInt {.importc: "LAPACKE_dstevd",
    dynlib: libName.}
proc dstevdWork*(matrixLayout: cint; jobz: char; n: LapackInt; d: ptr cdouble;
                e: ptr cdouble; z: ptr cdouble; ldz: LapackInt; work: ptr cdouble;
                lwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dstevd_work", dynlib: libName.}
proc dstevr*(matrixLayout: cint; jobz: char; range: char; n: LapackInt; d: ptr cdouble;
            e: ptr cdouble; vl: cdouble; vu: cdouble; il: LapackInt; iu: LapackInt;
            abstol: cdouble; m: ptr LapackInt; w: ptr cdouble; z: ptr cdouble;
            ldz: LapackInt; isuppz: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dstevr", dynlib: libName.}
proc dstevrWork*(matrixLayout: cint; jobz: char; range: char; n: LapackInt;
                d: ptr cdouble; e: ptr cdouble; vl: cdouble; vu: cdouble; il: LapackInt;
                iu: LapackInt; abstol: cdouble; m: ptr LapackInt; w: ptr cdouble;
                z: ptr cdouble; ldz: LapackInt; isuppz: ptr LapackInt;
                work: ptr cdouble; lwork: LapackInt; iwork: ptr LapackInt;
                liwork: LapackInt): LapackInt {.importc: "LAPACKE_dstevr_work",
    dynlib: libName.}
proc dstevx*(matrixLayout: cint; jobz: char; range: char; n: LapackInt; d: ptr cdouble;
            e: ptr cdouble; vl: cdouble; vu: cdouble; il: LapackInt; iu: LapackInt;
            abstol: cdouble; m: ptr LapackInt; w: ptr cdouble; z: ptr cdouble;
            ldz: LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dstevx", dynlib: libName.}
proc dstevxWork*(matrixLayout: cint; jobz: char; range: char; n: LapackInt;
                d: ptr cdouble; e: ptr cdouble; vl: cdouble; vu: cdouble; il: LapackInt;
                iu: LapackInt; abstol: cdouble; m: ptr LapackInt; w: ptr cdouble;
                z: ptr cdouble; ldz: LapackInt; work: ptr cdouble; iwork: ptr LapackInt;
                ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_dstevx_work",
    dynlib: libName.}
proc dsycon*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble; lda: LapackInt;
            ipiv: ptr LapackInt; anorm: cdouble; rcond: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsycon", dynlib: libName.}
proc dsyconWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; ipiv: ptr LapackInt; anorm: cdouble;
                rcond: ptr cdouble; work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsycon_work", dynlib: libName.}
proc dsyconv*(matrixLayout: cint; uplo: char; way: char; n: LapackInt; a: ptr cdouble;
             lda: LapackInt; ipiv: ptr LapackInt; e: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsyconv", dynlib: libName.}
proc dsyconvWork*(matrixLayout: cint; uplo: char; way: char; n: LapackInt;
                 a: ptr cdouble; lda: LapackInt; ipiv: ptr LapackInt; e: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsyconv_work", dynlib: libName.}
proc dsyequb*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
             lda: LapackInt; s: ptr cdouble; scond: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsyequb", dynlib: libName.}
proc dsyequbWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                 lda: LapackInt; s: ptr cdouble; scond: ptr cdouble; amax: ptr cdouble;
                 work: ptr cdouble): LapackInt {.importc: "LAPACKE_dsyequb_work",
    dynlib: libName.}
proc dsyev*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; a: ptr cdouble;
           lda: LapackInt; w: ptr cdouble): LapackInt {.importc: "LAPACKE_dsyev",
    dynlib: libName.}
proc dsyevWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; a: ptr cdouble;
               lda: LapackInt; w: ptr cdouble; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dsyev_work", dynlib: libName.}
proc dsyevd*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; a: ptr cdouble;
            lda: LapackInt; w: ptr cdouble): LapackInt {.importc: "LAPACKE_dsyevd",
    dynlib: libName.}
proc dsyevdWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                a: ptr cdouble; lda: LapackInt; w: ptr cdouble; work: ptr cdouble;
                lwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dsyevd_work", dynlib: libName.}
proc dsyevr*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
            a: ptr cdouble; lda: LapackInt; vl: cdouble; vu: cdouble; il: LapackInt;
            iu: LapackInt; abstol: cdouble; m: ptr LapackInt; w: ptr cdouble;
            z: ptr cdouble; ldz: LapackInt; isuppz: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsyevr", dynlib: libName.}
proc dsyevrWork*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                a: ptr cdouble; lda: LapackInt; vl: cdouble; vu: cdouble; il: LapackInt;
                iu: LapackInt; abstol: cdouble; m: ptr LapackInt; w: ptr cdouble;
                z: ptr cdouble; ldz: LapackInt; isuppz: ptr LapackInt;
                work: ptr cdouble; lwork: LapackInt; iwork: ptr LapackInt;
                liwork: LapackInt): LapackInt {.importc: "LAPACKE_dsyevr_work",
    dynlib: libName.}
proc dsyevx*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
            a: ptr cdouble; lda: LapackInt; vl: cdouble; vu: cdouble; il: LapackInt;
            iu: LapackInt; abstol: cdouble; m: ptr LapackInt; w: ptr cdouble;
            z: ptr cdouble; ldz: LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsyevx", dynlib: libName.}
proc dsyevxWork*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                a: ptr cdouble; lda: LapackInt; vl: cdouble; vu: cdouble; il: LapackInt;
                iu: LapackInt; abstol: cdouble; m: ptr LapackInt; w: ptr cdouble;
                z: ptr cdouble; ldz: LapackInt; work: ptr cdouble; lwork: LapackInt;
                iwork: ptr LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsyevx_work", dynlib: libName.}
proc dsygst*(matrixLayout: cint; itype: LapackInt; uplo: char; n: LapackInt;
            a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dsygst", dynlib: libName.}
proc dsygstWork*(matrixLayout: cint; itype: LapackInt; uplo: char; n: LapackInt;
                a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dsygst_work", dynlib: libName.}
proc dsygv*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char; n: LapackInt;
           a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt; w: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsygv", dynlib: libName.}
proc dsygvWork*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
               n: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble;
               ldb: LapackInt; w: ptr cdouble; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dsygv_work", dynlib: libName.}
proc dsygvd*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char; n: LapackInt;
            a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt; w: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsygvd", dynlib: libName.}
proc dsygvdWork*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
                n: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble;
                ldb: LapackInt; w: ptr cdouble; work: ptr cdouble; lwork: LapackInt;
                iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dsygvd_work", dynlib: libName.}
proc dsygvx*(matrixLayout: cint; itype: LapackInt; jobz: char; range: char; uplo: char;
            n: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
            vl: cdouble; vu: cdouble; il: LapackInt; iu: LapackInt; abstol: cdouble;
            m: ptr LapackInt; w: ptr cdouble; z: ptr cdouble; ldz: LapackInt;
            ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_dsygvx",
    dynlib: libName.}
proc dsygvxWork*(matrixLayout: cint; itype: LapackInt; jobz: char; range: char;
                uplo: char; n: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble;
                ldb: LapackInt; vl: cdouble; vu: cdouble; il: LapackInt; iu: LapackInt;
                abstol: cdouble; m: ptr LapackInt; w: ptr cdouble; z: ptr cdouble;
                ldz: LapackInt; work: ptr cdouble; lwork: LapackInt;
                iwork: ptr LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsygvx_work", dynlib: libName.}
proc dsyrfs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr cdouble; lda: LapackInt; af: ptr cdouble; ldaf: LapackInt;
            ipiv: ptr LapackInt; b: ptr cdouble; ldb: LapackInt; x: ptr cdouble;
            ldx: LapackInt; ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsyrfs", dynlib: libName.}
proc dsyrfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr cdouble; lda: LapackInt; af: ptr cdouble; ldaf: LapackInt;
                ipiv: ptr LapackInt; b: ptr cdouble; ldb: LapackInt; x: ptr cdouble;
                ldx: LapackInt; ferr: ptr cdouble; berr: ptr cdouble; work: ptr cdouble;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dsyrfs_work",
    dynlib: libName.}
proc dsyrfsx*(matrixLayout: cint; uplo: char; equed: char; n: LapackInt; nrhs: LapackInt;
             a: ptr cdouble; lda: LapackInt; af: ptr cdouble; ldaf: LapackInt;
             ipiv: ptr LapackInt; s: ptr cdouble; b: ptr cdouble; ldb: LapackInt;
             x: ptr cdouble; ldx: LapackInt; rcond: ptr cdouble; berr: ptr cdouble;
             nErrBnds: LapackInt; errBndsNorm: ptr cdouble; errBndsComp: ptr cdouble;
             nparams: LapackInt; params: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsyrfsx", dynlib: libName.}
proc dsyrfsxWork*(matrixLayout: cint; uplo: char; equed: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr cdouble; lda: LapackInt; af: ptr cdouble;
                 ldaf: LapackInt; ipiv: ptr LapackInt; s: ptr cdouble; b: ptr cdouble;
                 ldb: LapackInt; x: ptr cdouble; ldx: LapackInt; rcond: ptr cdouble;
                 berr: ptr cdouble; nErrBnds: LapackInt; errBndsNorm: ptr cdouble;
                 errBndsComp: ptr cdouble; nparams: LapackInt; params: ptr cdouble;
                 work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsyrfsx_work", dynlib: libName.}
proc dsysv*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
           a: ptr cdouble; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cdouble;
           ldb: LapackInt): LapackInt {.importc: "LAPACKE_dsysv", dynlib: libName.}
proc dsysvRook*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               a: ptr cdouble; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cdouble;
               ldb: LapackInt): LapackInt {.importc: "LAPACKE_dsysv_rook",
    dynlib: libName.}
proc dsysvRookWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                   a: ptr cdouble; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cdouble;
                   ldb: LapackInt; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dsysv_rook_work", dynlib: libName.}
proc dsysvWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               a: ptr cdouble; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cdouble;
               ldb: LapackInt; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dsysv_work", dynlib: libName.}
proc dsysvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr cdouble; lda: LapackInt; af: ptr cdouble; ldaf: LapackInt;
            ipiv: ptr LapackInt; b: ptr cdouble; ldb: LapackInt; x: ptr cdouble;
            ldx: LapackInt; rcond: ptr cdouble; ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsysvx", dynlib: libName.}
proc dsysvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; a: ptr cdouble; lda: LapackInt; af: ptr cdouble;
                ldaf: LapackInt; ipiv: ptr LapackInt; b: ptr cdouble; ldb: LapackInt;
                x: ptr cdouble; ldx: LapackInt; rcond: ptr cdouble; ferr: ptr cdouble;
                berr: ptr cdouble; work: ptr cdouble; lwork: LapackInt;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dsysvx_work",
    dynlib: libName.}
proc dsysvxx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr cdouble; lda: LapackInt; af: ptr cdouble; ldaf: LapackInt;
             ipiv: ptr LapackInt; equed: cstring; s: ptr cdouble; b: ptr cdouble;
             ldb: LapackInt; x: ptr cdouble; ldx: LapackInt; rcond: ptr cdouble;
             rpvgrw: ptr cdouble; berr: ptr cdouble; nErrBnds: LapackInt;
             errBndsNorm: ptr cdouble; errBndsComp: ptr cdouble; nparams: LapackInt;
             params: ptr cdouble): LapackInt {.importc: "LAPACKE_dsysvxx",
    dynlib: libName.}
proc dsysvxxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr cdouble; lda: LapackInt; af: ptr cdouble;
                 ldaf: LapackInt; ipiv: ptr LapackInt; equed: cstring; s: ptr cdouble;
                 b: ptr cdouble; ldb: LapackInt; x: ptr cdouble; ldx: LapackInt;
                 rcond: ptr cdouble; rpvgrw: ptr cdouble; berr: ptr cdouble;
                 nErrBnds: LapackInt; errBndsNorm: ptr cdouble;
                 errBndsComp: ptr cdouble; nparams: LapackInt; params: ptr cdouble;
                 work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsysvxx_work", dynlib: libName.}
proc dsyswapr*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
              lda: LapackInt; i1: LapackInt; i2: LapackInt): LapackInt {.
    importc: "LAPACKE_dsyswapr", dynlib: libName.}
proc dsyswaprWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                  lda: LapackInt; i1: LapackInt; i2: LapackInt): LapackInt {.
    importc: "LAPACKE_dsyswapr_work", dynlib: libName.}
proc dsytrd*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble; lda: LapackInt;
            d: ptr cdouble; e: ptr cdouble; tau: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsytrd", dynlib: libName.}
proc dsytrdWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; d: ptr cdouble; e: ptr cdouble; tau: ptr cdouble;
                work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dsytrd_work", dynlib: libName.}
proc dsytrf*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble; lda: LapackInt;
            ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_dsytrf",
    dynlib: libName.}
proc dsytrfRook*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsytrf_rook", dynlib: libName.}
proc dsytrfRookWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                    lda: LapackInt; ipiv: ptr LapackInt; work: ptr cdouble;
                    lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dsytrf_rook_work", dynlib: libName.}
proc dsytrfWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; ipiv: ptr LapackInt; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dsytrf_work", dynlib: libName.}
proc dsytri*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble; lda: LapackInt;
            ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_dsytri",
    dynlib: libName.}
proc dsytri2*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
             lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsytri2", dynlib: libName.}
proc dsytri2Work*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                 lda: LapackInt; ipiv: ptr LapackInt; work: ptr cdouble;
                 lwork: LapackInt): LapackInt {.importc: "LAPACKE_dsytri2_work",
    dynlib: libName.}
proc dsytri2x*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
              lda: LapackInt; ipiv: ptr LapackInt; nb: LapackInt): LapackInt {.
    importc: "LAPACKE_dsytri2x", dynlib: libName.}
proc dsytri2xWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                  lda: LapackInt; ipiv: ptr LapackInt; work: ptr cdouble; nb: LapackInt): LapackInt {.
    importc: "LAPACKE_dsytri2x_work", dynlib: libName.}
proc dsytriWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; ipiv: ptr LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsytri_work", dynlib: libName.}
proc dsytrs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr cdouble; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cdouble;
            ldb: LapackInt): LapackInt {.importc: "LAPACKE_dsytrs", dynlib: libName.}
proc dsytrs2*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr cdouble; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cdouble;
             ldb: LapackInt): LapackInt {.importc: "LAPACKE_dsytrs2", dynlib: libName.}
proc dsytrs2Work*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr cdouble; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cdouble;
                 ldb: LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsytrs2_work", dynlib: libName.}
proc dsytrsRook*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr cdouble; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cdouble;
                ldb: LapackInt): LapackInt {.importc: "LAPACKE_dsytrs_rook",
    dynlib: libName.}
proc dsytrsRookWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                    a: ptr cdouble; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cdouble;
                    ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dsytrs_rook_work", dynlib: libName.}
proc dsytrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr cdouble; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cdouble;
                ldb: LapackInt): LapackInt {.importc: "LAPACKE_dsytrs_work",
    dynlib: libName.}
proc dtbcon*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
            kd: LapackInt; ab: ptr cdouble; ldab: LapackInt; rcond: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dtbcon", dynlib: libName.}
proc dtbconWork*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
                kd: LapackInt; ab: ptr cdouble; ldab: LapackInt; rcond: ptr cdouble;
                work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dtbcon_work", dynlib: libName.}
proc dtbrfs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            kd: LapackInt; nrhs: LapackInt; ab: ptr cdouble; ldab: LapackInt;
            b: ptr cdouble; ldb: LapackInt; x: ptr cdouble; ldx: LapackInt;
            ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dtbrfs", dynlib: libName.}
proc dtbrfsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                kd: LapackInt; nrhs: LapackInt; ab: ptr cdouble; ldab: LapackInt;
                b: ptr cdouble; ldb: LapackInt; x: ptr cdouble; ldx: LapackInt;
                ferr: ptr cdouble; berr: ptr cdouble; work: ptr cdouble;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dtbrfs_work",
    dynlib: libName.}
proc dtbtrs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            kd: LapackInt; nrhs: LapackInt; ab: ptr cdouble; ldab: LapackInt;
            b: ptr cdouble; ldb: LapackInt): LapackInt {.importc: "LAPACKE_dtbtrs",
    dynlib: libName.}
proc dtbtrsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                kd: LapackInt; nrhs: LapackInt; ab: ptr cdouble; ldab: LapackInt;
                b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dtbtrs_work", dynlib: libName.}
proc dtfsm*(matrixLayout: cint; transr: char; side: char; uplo: char; trans: char;
           diag: char; m: LapackInt; n: LapackInt; alpha: cdouble; a: ptr cdouble;
           b: ptr cdouble; ldb: LapackInt): LapackInt {.importc: "LAPACKE_dtfsm",
    dynlib: libName.}
proc dtfsmWork*(matrixLayout: cint; transr: char; side: char; uplo: char; trans: char;
               diag: char; m: LapackInt; n: LapackInt; alpha: cdouble; a: ptr cdouble;
               b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dtfsm_work", dynlib: libName.}
proc dtftri*(matrixLayout: cint; transr: char; uplo: char; diag: char; n: LapackInt;
            a: ptr cdouble): LapackInt {.importc: "LAPACKE_dtftri", dynlib: libName.}
proc dtftriWork*(matrixLayout: cint; transr: char; uplo: char; diag: char; n: LapackInt;
                a: ptr cdouble): LapackInt {.importc: "LAPACKE_dtftri_work",
    dynlib: libName.}
proc dtfttp*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
            arf: ptr cdouble; ap: ptr cdouble): LapackInt {.importc: "LAPACKE_dtfttp",
    dynlib: libName.}
proc dtfttpWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                arf: ptr cdouble; ap: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dtfttp_work", dynlib: libName.}
proc dtfttr*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
            arf: ptr cdouble; a: ptr cdouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_dtfttr", dynlib: libName.}
proc dtfttrWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                arf: ptr cdouble; a: ptr cdouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_dtfttr_work", dynlib: libName.}
proc dtgevc*(matrixLayout: cint; side: char; howmny: char; select: ptr LapackLogical;
            n: LapackInt; s: ptr cdouble; lds: LapackInt; p: ptr cdouble; ldp: LapackInt;
            vl: ptr cdouble; ldvl: LapackInt; vr: ptr cdouble; ldvr: LapackInt;
            mm: LapackInt; m: ptr LapackInt): LapackInt {.importc: "LAPACKE_dtgevc",
    dynlib: libName.}
proc dtgevcWork*(matrixLayout: cint; side: char; howmny: char;
                select: ptr LapackLogical; n: LapackInt; s: ptr cdouble; lds: LapackInt;
                p: ptr cdouble; ldp: LapackInt; vl: ptr cdouble; ldvl: LapackInt;
                vr: ptr cdouble; ldvr: LapackInt; mm: LapackInt; m: ptr LapackInt;
                work: ptr cdouble): LapackInt {.importc: "LAPACKE_dtgevc_work",
    dynlib: libName.}
proc dtgexc*(matrixLayout: cint; wantq: LapackLogical; wantz: LapackLogical;
            n: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
            q: ptr cdouble; ldq: LapackInt; z: ptr cdouble; ldz: LapackInt;
            ifst: ptr LapackInt; ilst: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dtgexc", dynlib: libName.}
proc dtgexcWork*(matrixLayout: cint; wantq: LapackLogical; wantz: LapackLogical;
                n: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble;
                ldb: LapackInt; q: ptr cdouble; ldq: LapackInt; z: ptr cdouble;
                ldz: LapackInt; ifst: ptr LapackInt; ilst: ptr LapackInt;
                work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dtgexc_work", dynlib: libName.}
proc dtgsen*(matrixLayout: cint; ijob: LapackInt; wantq: LapackLogical;
            wantz: LapackLogical; select: ptr LapackLogical; n: LapackInt;
            a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
            alphar: ptr cdouble; alphai: ptr cdouble; beta: ptr cdouble; q: ptr cdouble;
            ldq: LapackInt; z: ptr cdouble; ldz: LapackInt; m: ptr LapackInt;
            pl: ptr cdouble; pr: ptr cdouble; dif: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dtgsen", dynlib: libName.}
proc dtgsenWork*(matrixLayout: cint; ijob: LapackInt; wantq: LapackLogical;
                wantz: LapackLogical; select: ptr LapackLogical; n: LapackInt;
                a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
                alphar: ptr cdouble; alphai: ptr cdouble; beta: ptr cdouble;
                q: ptr cdouble; ldq: LapackInt; z: ptr cdouble; ldz: LapackInt;
                m: ptr LapackInt; pl: ptr cdouble; pr: ptr cdouble; dif: ptr cdouble;
                work: ptr cdouble; lwork: LapackInt; iwork: ptr LapackInt;
                liwork: LapackInt): LapackInt {.importc: "LAPACKE_dtgsen_work",
    dynlib: libName.}
proc dtgsja*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
            p: LapackInt; n: LapackInt; k: LapackInt; l: LapackInt; a: ptr cdouble;
            lda: LapackInt; b: ptr cdouble; ldb: LapackInt; tola: cdouble; tolb: cdouble;
            alpha: ptr cdouble; beta: ptr cdouble; u: ptr cdouble; ldu: LapackInt;
            v: ptr cdouble; ldv: LapackInt; q: ptr cdouble; ldq: LapackInt;
            ncycle: ptr LapackInt): LapackInt {.importc: "LAPACKE_dtgsja",
    dynlib: libName.}
proc dtgsjaWork*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
                p: LapackInt; n: LapackInt; k: LapackInt; l: LapackInt; a: ptr cdouble;
                lda: LapackInt; b: ptr cdouble; ldb: LapackInt; tola: cdouble;
                tolb: cdouble; alpha: ptr cdouble; beta: ptr cdouble; u: ptr cdouble;
                ldu: LapackInt; v: ptr cdouble; ldv: LapackInt; q: ptr cdouble;
                ldq: LapackInt; work: ptr cdouble; ncycle: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dtgsja_work", dynlib: libName.}
proc dtgsna*(matrixLayout: cint; job: char; howmny: char; select: ptr LapackLogical;
            n: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
            vl: ptr cdouble; ldvl: LapackInt; vr: ptr cdouble; ldvr: LapackInt;
            s: ptr cdouble; dif: ptr cdouble; mm: LapackInt; m: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dtgsna", dynlib: libName.}
proc dtgsnaWork*(matrixLayout: cint; job: char; howmny: char;
                select: ptr LapackLogical; n: LapackInt; a: ptr cdouble; lda: LapackInt;
                b: ptr cdouble; ldb: LapackInt; vl: ptr cdouble; ldvl: LapackInt;
                vr: ptr cdouble; ldvr: LapackInt; s: ptr cdouble; dif: ptr cdouble;
                mm: LapackInt; m: ptr LapackInt; work: ptr cdouble; lwork: LapackInt;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dtgsna_work",
    dynlib: libName.}
proc dtgsyl*(matrixLayout: cint; trans: char; ijob: LapackInt; m: LapackInt;
            n: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
            c: ptr cdouble; ldc: LapackInt; d: ptr cdouble; ldd: LapackInt;
            e: ptr cdouble; lde: LapackInt; f: ptr cdouble; ldf: LapackInt;
            scale: ptr cdouble; dif: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dtgsyl", dynlib: libName.}
proc dtgsylWork*(matrixLayout: cint; trans: char; ijob: LapackInt; m: LapackInt;
                n: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble;
                ldb: LapackInt; c: ptr cdouble; ldc: LapackInt; d: ptr cdouble;
                ldd: LapackInt; e: ptr cdouble; lde: LapackInt; f: ptr cdouble;
                ldf: LapackInt; scale: ptr cdouble; dif: ptr cdouble; work: ptr cdouble;
                lwork: LapackInt; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dtgsyl_work", dynlib: libName.}
proc dtpcon*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
            ap: ptr cdouble; rcond: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dtpcon", dynlib: libName.}
proc dtpconWork*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
                ap: ptr cdouble; rcond: ptr cdouble; work: ptr cdouble;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dtpcon_work",
    dynlib: libName.}
proc dtpmqrt*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
             k: LapackInt; l: LapackInt; nb: LapackInt; v: ptr cdouble; ldv: LapackInt;
             t: ptr cdouble; ldt: LapackInt; a: ptr cdouble; lda: LapackInt;
             b: ptr cdouble; ldb: LapackInt): LapackInt {.importc: "LAPACKE_dtpmqrt",
    dynlib: libName.}
proc dtpmqrtWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt;
                 n: LapackInt; k: LapackInt; l: LapackInt; nb: LapackInt;
                 v: ptr cdouble; ldv: LapackInt; t: ptr cdouble; ldt: LapackInt;
                 a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
                 work: ptr cdouble): LapackInt {.importc: "LAPACKE_dtpmqrt_work",
    dynlib: libName.}
proc dtpqrt*(matrixLayout: cint; m: LapackInt; n: LapackInt; l: LapackInt; nb: LapackInt;
            a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
            t: ptr cdouble; ldt: LapackInt): LapackInt {.importc: "LAPACKE_dtpqrt",
    dynlib: libName.}
proc dtpqrt2*(matrixLayout: cint; m: LapackInt; n: LapackInt; l: LapackInt;
             a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
             t: ptr cdouble; ldt: LapackInt): LapackInt {.importc: "LAPACKE_dtpqrt2",
    dynlib: libName.}
proc dtpqrt2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt; l: LapackInt;
                 a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
                 t: ptr cdouble; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_dtpqrt2_work", dynlib: libName.}
proc dtpqrtWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; l: LapackInt;
                nb: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble;
                ldb: LapackInt; t: ptr cdouble; ldt: LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dtpqrt_work", dynlib: libName.}
proc dtprfb*(matrixLayout: cint; side: char; trans: char; direct: char; storev: char;
            m: LapackInt; n: LapackInt; k: LapackInt; l: LapackInt; v: ptr cdouble;
            ldv: LapackInt; t: ptr cdouble; ldt: LapackInt; a: ptr cdouble;
            lda: LapackInt; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dtprfb", dynlib: libName.}
proc dtprfbWork*(matrixLayout: cint; side: char; trans: char; direct: char; storev: char;
                m: LapackInt; n: LapackInt; k: LapackInt; l: LapackInt; v: ptr cdouble;
                ldv: LapackInt; t: ptr cdouble; ldt: LapackInt; a: ptr cdouble;
                lda: LapackInt; b: ptr cdouble; ldb: LapackInt; work: ptr cdouble;
                ldwork: LapackInt): LapackInt {.importc: "LAPACKE_dtprfb_work",
    dynlib: libName.}
proc dtprfs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            nrhs: LapackInt; ap: ptr cdouble; b: ptr cdouble; ldb: LapackInt;
            x: ptr cdouble; ldx: LapackInt; ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dtprfs", dynlib: libName.}
proc dtprfsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                nrhs: LapackInt; ap: ptr cdouble; b: ptr cdouble; ldb: LapackInt;
                x: ptr cdouble; ldx: LapackInt; ferr: ptr cdouble; berr: ptr cdouble;
                work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dtprfs_work", dynlib: libName.}
proc dtptri*(matrixLayout: cint; uplo: char; diag: char; n: LapackInt; ap: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dtptri", dynlib: libName.}
proc dtptriWork*(matrixLayout: cint; uplo: char; diag: char; n: LapackInt;
                ap: ptr cdouble): LapackInt {.importc: "LAPACKE_dtptri_work",
    dynlib: libName.}
proc dtptrs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            nrhs: LapackInt; ap: ptr cdouble; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dtptrs", dynlib: libName.}
proc dtptrsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                nrhs: LapackInt; ap: ptr cdouble; b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dtptrs_work", dynlib: libName.}
proc dtpttf*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt; ap: ptr cdouble;
            arf: ptr cdouble): LapackInt {.importc: "LAPACKE_dtpttf", dynlib: libName.}
proc dtpttfWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                ap: ptr cdouble; arf: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dtpttf_work", dynlib: libName.}
proc dtpttr*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cdouble;
            a: ptr cdouble; lda: LapackInt): LapackInt {.importc: "LAPACKE_dtpttr",
    dynlib: libName.}
proc dtpttrWork*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cdouble;
                a: ptr cdouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_dtpttr_work", dynlib: libName.}
proc dtrcon*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
            a: ptr cdouble; lda: LapackInt; rcond: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dtrcon", dynlib: libName.}
proc dtrconWork*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
                a: ptr cdouble; lda: LapackInt; rcond: ptr cdouble; work: ptr cdouble;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_dtrcon_work",
    dynlib: libName.}
proc dtrevc*(matrixLayout: cint; side: char; howmny: char; select: ptr LapackLogical;
            n: LapackInt; t: ptr cdouble; ldt: LapackInt; vl: ptr cdouble;
            ldvl: LapackInt; vr: ptr cdouble; ldvr: LapackInt; mm: LapackInt;
            m: ptr LapackInt): LapackInt {.importc: "LAPACKE_dtrevc", dynlib: libName.}
proc dtrevcWork*(matrixLayout: cint; side: char; howmny: char;
                select: ptr LapackLogical; n: LapackInt; t: ptr cdouble; ldt: LapackInt;
                vl: ptr cdouble; ldvl: LapackInt; vr: ptr cdouble; ldvr: LapackInt;
                mm: LapackInt; m: ptr LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dtrevc_work", dynlib: libName.}
proc dtrexc*(matrixLayout: cint; compq: char; n: LapackInt; t: ptr cdouble;
            ldt: LapackInt; q: ptr cdouble; ldq: LapackInt; ifst: ptr LapackInt;
            ilst: ptr LapackInt): LapackInt {.importc: "LAPACKE_dtrexc",
    dynlib: libName.}
proc dtrexcWork*(matrixLayout: cint; compq: char; n: LapackInt; t: ptr cdouble;
                ldt: LapackInt; q: ptr cdouble; ldq: LapackInt; ifst: ptr LapackInt;
                ilst: ptr LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dtrexc_work", dynlib: libName.}
proc dtrrfs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            nrhs: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble;
            ldb: LapackInt; x: ptr cdouble; ldx: LapackInt; ferr: ptr cdouble;
            berr: ptr cdouble): LapackInt {.importc: "LAPACKE_dtrrfs", dynlib: libName.}
proc dtrrfsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                nrhs: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble;
                ldb: LapackInt; x: ptr cdouble; ldx: LapackInt; ferr: ptr cdouble;
                berr: ptr cdouble; work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dtrrfs_work", dynlib: libName.}
proc dtrsen*(matrixLayout: cint; job: char; compq: char; select: ptr LapackLogical;
            n: LapackInt; t: ptr cdouble; ldt: LapackInt; q: ptr cdouble; ldq: LapackInt;
            wr: ptr cdouble; wi: ptr cdouble; m: ptr LapackInt; s: ptr cdouble;
            sep: ptr cdouble): LapackInt {.importc: "LAPACKE_dtrsen", dynlib: libName.}
proc dtrsenWork*(matrixLayout: cint; job: char; compq: char; select: ptr LapackLogical;
                n: LapackInt; t: ptr cdouble; ldt: LapackInt; q: ptr cdouble;
                ldq: LapackInt; wr: ptr cdouble; wi: ptr cdouble; m: ptr LapackInt;
                s: ptr cdouble; sep: ptr cdouble; work: ptr cdouble; lwork: LapackInt;
                iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dtrsen_work", dynlib: libName.}
proc dtrsna*(matrixLayout: cint; job: char; howmny: char; select: ptr LapackLogical;
            n: LapackInt; t: ptr cdouble; ldt: LapackInt; vl: ptr cdouble;
            ldvl: LapackInt; vr: ptr cdouble; ldvr: LapackInt; s: ptr cdouble;
            sep: ptr cdouble; mm: LapackInt; m: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dtrsna", dynlib: libName.}
proc dtrsnaWork*(matrixLayout: cint; job: char; howmny: char;
                select: ptr LapackLogical; n: LapackInt; t: ptr cdouble; ldt: LapackInt;
                vl: ptr cdouble; ldvl: LapackInt; vr: ptr cdouble; ldvr: LapackInt;
                s: ptr cdouble; sep: ptr cdouble; mm: LapackInt; m: ptr LapackInt;
                work: ptr cdouble; ldwork: LapackInt; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dtrsna_work", dynlib: libName.}
proc dtrsyl*(matrixLayout: cint; trana: char; tranb: char; isgn: LapackInt; m: LapackInt;
            n: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble; ldb: LapackInt;
            c: ptr cdouble; ldc: LapackInt; scale: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dtrsyl", dynlib: libName.}
proc dtrsylWork*(matrixLayout: cint; trana: char; tranb: char; isgn: LapackInt;
                m: LapackInt; n: LapackInt; a: ptr cdouble; lda: LapackInt;
                b: ptr cdouble; ldb: LapackInt; c: ptr cdouble; ldc: LapackInt;
                scale: ptr cdouble): LapackInt {.importc: "LAPACKE_dtrsyl_work",
    dynlib: libName.}
proc dtrtri*(matrixLayout: cint; uplo: char; diag: char; n: LapackInt; a: ptr cdouble;
            lda: LapackInt): LapackInt {.importc: "LAPACKE_dtrtri", dynlib: libName.}
proc dtrtriWork*(matrixLayout: cint; uplo: char; diag: char; n: LapackInt;
                a: ptr cdouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_dtrtri_work", dynlib: libName.}
proc dtrtrs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            nrhs: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble;
            ldb: LapackInt): LapackInt {.importc: "LAPACKE_dtrtrs", dynlib: libName.}
proc dtrtrsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                nrhs: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble;
                ldb: LapackInt): LapackInt {.importc: "LAPACKE_dtrtrs_work",
    dynlib: libName.}
proc dtrttf*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt; a: ptr cdouble;
            lda: LapackInt; arf: ptr cdouble): LapackInt {.importc: "LAPACKE_dtrttf",
    dynlib: libName.}
proc dtrttfWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                a: ptr cdouble; lda: LapackInt; arf: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dtrttf_work", dynlib: libName.}
proc dtrttp*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble; lda: LapackInt;
            ap: ptr cdouble): LapackInt {.importc: "LAPACKE_dtrttp", dynlib: libName.}
proc dtrttpWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; ap: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dtrttp_work", dynlib: libName.}
proc dtzrzf*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
            lda: LapackInt; tau: ptr cdouble): LapackInt {.importc: "LAPACKE_dtzrzf",
    dynlib: libName.}
proc dtzrzfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
                lda: LapackInt; tau: ptr cdouble; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dtzrzf_work", dynlib: libName.}
proc mklCgetrfnpi*(matrixLayout: cint; m: LapackInt; n: LapackInt; nfact: LapackInt;
                  a: ptr LapackComplexFloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_cgetrfnpi", dynlib: libName.}
proc mklCgetrfnpiWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                      nfact: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_cgetrfnpi_work", dynlib: libName.}
proc mklCtppack*(matrixLayout: cint; uplo: char; trans: char; n: LapackInt;
                ap: ptr LapackComplexFloat; i: LapackInt; j: LapackInt;
                rows: LapackInt; cols: LapackInt; a: ptr LapackComplexFloat;
                lda: LapackInt): LapackInt {.importc: "LAPACKE_mkl_ctppack",
    dynlib: libName.}
proc mklCtppackWork*(matrixLayout: cint; uplo: char; trans: char; n: LapackInt;
                    ap: ptr LapackComplexFloat; i: LapackInt; j: LapackInt;
                    rows: LapackInt; cols: LapackInt; a: ptr LapackComplexFloat;
                    lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_ctppack_work", dynlib: libName.}
proc mklCtpunpack*(matrixLayout: cint; uplo: char; trans: char; n: LapackInt;
                  ap: ptr LapackComplexFloat; i: LapackInt; j: LapackInt;
                  rows: LapackInt; cols: LapackInt; a: ptr LapackComplexFloat;
                  lda: LapackInt): LapackInt {.importc: "LAPACKE_mkl_ctpunpack",
    dynlib: libName.}
proc mklCtpunpackWork*(matrixLayout: cint; uplo: char; trans: char; n: LapackInt;
                      ap: ptr LapackComplexFloat; i: LapackInt; j: LapackInt;
                      rows: LapackInt; cols: LapackInt; a: ptr LapackComplexFloat;
                      lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_ctpunpack_work", dynlib: libName.}
proc mklDgetrfnpi*(matrixLayout: cint; m: LapackInt; n: LapackInt; nfact: LapackInt;
                  a: ptr cdouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_dgetrfnpi", dynlib: libName.}
proc mklDgetrfnpiWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                      nfact: LapackInt; a: ptr cdouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_dgetrfnpi_work", dynlib: libName.}
proc mklDtppack*(matrixLayout: cint; uplo: char; trans: char; n: LapackInt;
                ap: ptr cdouble; i: LapackInt; j: LapackInt; rows: LapackInt;
                cols: LapackInt; a: ptr cdouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_dtppack", dynlib: libName.}
proc mklDtppackWork*(matrixLayout: cint; uplo: char; trans: char; n: LapackInt;
                    ap: ptr cdouble; i: LapackInt; j: LapackInt; rows: LapackInt;
                    cols: LapackInt; a: ptr cdouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_dtppack_work", dynlib: libName.}
proc mklDtpunpack*(matrixLayout: cint; uplo: char; trans: char; n: LapackInt;
                  ap: ptr cdouble; i: LapackInt; j: LapackInt; rows: LapackInt;
                  cols: LapackInt; a: ptr cdouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_dtpunpack", dynlib: libName.}
proc mklDtpunpackWork*(matrixLayout: cint; uplo: char; trans: char; n: LapackInt;
                      ap: ptr cdouble; i: LapackInt; j: LapackInt; rows: LapackInt;
                      cols: LapackInt; a: ptr cdouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_dtpunpack_work", dynlib: libName.}
proc mklSgetrfnpi*(matrixLayout: cint; m: LapackInt; n: LapackInt; nfact: LapackInt;
                  a: ptr cfloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_sgetrfnpi", dynlib: libName.}
proc mklSgetrfnpiWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                      nfact: LapackInt; a: ptr cfloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_sgetrfnpi_work", dynlib: libName.}
proc mklStppack*(matrixLayout: cint; uplo: char; trans: char; n: LapackInt;
                ap: ptr cfloat; i: LapackInt; j: LapackInt; rows: LapackInt;
                cols: LapackInt; a: ptr cfloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_stppack", dynlib: libName.}
proc mklStppackWork*(matrixLayout: cint; uplo: char; trans: char; n: LapackInt;
                    ap: ptr cfloat; i: LapackInt; j: LapackInt; rows: LapackInt;
                    cols: LapackInt; a: ptr cfloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_stppack_work", dynlib: libName.}
proc mklStpunpack*(matrixLayout: cint; uplo: char; trans: char; n: LapackInt;
                  ap: ptr cfloat; i: LapackInt; j: LapackInt; rows: LapackInt;
                  cols: LapackInt; a: ptr cfloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_stpunpack", dynlib: libName.}
proc mklStpunpackWork*(matrixLayout: cint; uplo: char; trans: char; n: LapackInt;
                      ap: ptr cfloat; i: LapackInt; j: LapackInt; rows: LapackInt;
                      cols: LapackInt; a: ptr cfloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_stpunpack_work", dynlib: libName.}
proc mklZgetrfnpi*(matrixLayout: cint; m: LapackInt; n: LapackInt; nfact: LapackInt;
                  a: ptr LapackComplexDouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_zgetrfnpi", dynlib: libName.}
proc mklZgetrfnpiWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                      nfact: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_zgetrfnpi_work", dynlib: libName.}
proc mklZtppack*(matrixLayout: cint; uplo: char; trans: char; n: LapackInt;
                ap: ptr LapackComplexDouble; i: LapackInt; j: LapackInt;
                rows: LapackInt; cols: LapackInt; a: ptr LapackComplexDouble;
                lda: LapackInt): LapackInt {.importc: "LAPACKE_mkl_ztppack",
    dynlib: libName.}
proc mklZtppackWork*(matrixLayout: cint; uplo: char; trans: char; n: LapackInt;
                    ap: ptr LapackComplexDouble; i: LapackInt; j: LapackInt;
                    rows: LapackInt; cols: LapackInt; a: ptr LapackComplexDouble;
                    lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_ztppack_work", dynlib: libName.}
proc mklZtpunpack*(matrixLayout: cint; uplo: char; trans: char; n: LapackInt;
                  ap: ptr LapackComplexDouble; i: LapackInt; j: LapackInt;
                  rows: LapackInt; cols: LapackInt; a: ptr LapackComplexDouble;
                  lda: LapackInt): LapackInt {.importc: "LAPACKE_mkl_ztpunpack",
    dynlib: libName.}
proc mklZtpunpackWork*(matrixLayout: cint; uplo: char; trans: char; n: LapackInt;
                      ap: ptr LapackComplexDouble; i: LapackInt; j: LapackInt;
                      rows: LapackInt; cols: LapackInt; a: ptr LapackComplexDouble;
                      lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_ztpunpack_work", dynlib: libName.}
proc sbbcsd*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char; jobv2t: char;
            trans: char; m: LapackInt; p: LapackInt; q: LapackInt; theta: ptr cfloat;
            phi: ptr cfloat; u1: ptr cfloat; ldu1: LapackInt; u2: ptr cfloat;
            ldu2: LapackInt; v1t: ptr cfloat; ldv1t: LapackInt; v2t: ptr cfloat;
            ldv2t: LapackInt; b11d: ptr cfloat; b11e: ptr cfloat; b12d: ptr cfloat;
            b12e: ptr cfloat; b21d: ptr cfloat; b21e: ptr cfloat; b22d: ptr cfloat;
            b22e: ptr cfloat): LapackInt {.importc: "LAPACKE_sbbcsd", dynlib: libName.}
proc sbbcsdWork*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char;
                jobv2t: char; trans: char; m: LapackInt; p: LapackInt; q: LapackInt;
                theta: ptr cfloat; phi: ptr cfloat; u1: ptr cfloat; ldu1: LapackInt;
                u2: ptr cfloat; ldu2: LapackInt; v1t: ptr cfloat; ldv1t: LapackInt;
                v2t: ptr cfloat; ldv2t: LapackInt; b11d: ptr cfloat; b11e: ptr cfloat;
                b12d: ptr cfloat; b12e: ptr cfloat; b21d: ptr cfloat; b21e: ptr cfloat;
                b22d: ptr cfloat; b22e: ptr cfloat; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sbbcsd_work", dynlib: libName.}
proc sbdsdc*(matrixLayout: cint; uplo: char; compq: char; n: LapackInt; d: ptr cfloat;
            e: ptr cfloat; u: ptr cfloat; ldu: LapackInt; vt: ptr cfloat; ldvt: LapackInt;
            q: ptr cfloat; iq: ptr LapackInt): LapackInt {.importc: "LAPACKE_sbdsdc",
    dynlib: libName.}
proc sbdsdcWork*(matrixLayout: cint; uplo: char; compq: char; n: LapackInt;
                d: ptr cfloat; e: ptr cfloat; u: ptr cfloat; ldu: LapackInt;
                vt: ptr cfloat; ldvt: LapackInt; q: ptr cfloat; iq: ptr LapackInt;
                work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sbdsdc_work", dynlib: libName.}
proc sbdsqr*(matrixLayout: cint; uplo: char; n: LapackInt; ncvt: LapackInt;
            nru: LapackInt; ncc: LapackInt; d: ptr cfloat; e: ptr cfloat; vt: ptr cfloat;
            ldvt: LapackInt; u: ptr cfloat; ldu: LapackInt; c: ptr cfloat; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_sbdsqr", dynlib: libName.}
proc sbdsqrWork*(matrixLayout: cint; uplo: char; n: LapackInt; ncvt: LapackInt;
                nru: LapackInt; ncc: LapackInt; d: ptr cfloat; e: ptr cfloat;
                vt: ptr cfloat; ldvt: LapackInt; u: ptr cfloat; ldu: LapackInt;
                c: ptr cfloat; ldc: LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sbdsqr_work", dynlib: libName.}
proc sbdsvdx*(matrixLayout: cint; uplo: char; jobz: char; range: char; n: LapackInt;
             d: ptr cfloat; e: ptr cfloat; vl: cfloat; vu: cfloat; il: LapackInt;
             iu: LapackInt; ns: ptr LapackInt; s: ptr cfloat; z: ptr cfloat;
             ldz: LapackInt; superb: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sbdsvdx", dynlib: libName.}
proc sbdsvdxWork*(matrixLayout: cint; uplo: char; jobz: char; range: char; n: LapackInt;
                 d: ptr cfloat; e: ptr cfloat; vl: cfloat; vu: cfloat; il: LapackInt;
                 iu: LapackInt; ns: ptr LapackInt; s: ptr cfloat; z: ptr cfloat;
                 ldz: LapackInt; work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sbdsvdx_work", dynlib: libName.}
proc sdisna*(job: char; m: LapackInt; n: LapackInt; d: ptr cfloat; sep: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sdisna", dynlib: libName.}
proc sdisnaWork*(job: char; m: LapackInt; n: LapackInt; d: ptr cfloat; sep: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sdisna_work", dynlib: libName.}
proc sgbbrd*(matrixLayout: cint; vect: char; m: LapackInt; n: LapackInt; ncc: LapackInt;
            kl: LapackInt; ku: LapackInt; ab: ptr cfloat; ldab: LapackInt; d: ptr cfloat;
            e: ptr cfloat; q: ptr cfloat; ldq: LapackInt; pt: ptr cfloat; ldpt: LapackInt;
            c: ptr cfloat; ldc: LapackInt): LapackInt {.importc: "LAPACKE_sgbbrd",
    dynlib: libName.}
proc sgbbrdWork*(matrixLayout: cint; vect: char; m: LapackInt; n: LapackInt;
                ncc: LapackInt; kl: LapackInt; ku: LapackInt; ab: ptr cfloat;
                ldab: LapackInt; d: ptr cfloat; e: ptr cfloat; q: ptr cfloat;
                ldq: LapackInt; pt: ptr cfloat; ldpt: LapackInt; c: ptr cfloat;
                ldc: LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgbbrd_work", dynlib: libName.}
proc sgbcon*(matrixLayout: cint; norm: char; n: LapackInt; kl: LapackInt; ku: LapackInt;
            ab: ptr cfloat; ldab: LapackInt; ipiv: ptr LapackInt; anorm: cfloat;
            rcond: ptr cfloat): LapackInt {.importc: "LAPACKE_sgbcon", dynlib: libName.}
proc sgbconWork*(matrixLayout: cint; norm: char; n: LapackInt; kl: LapackInt;
                ku: LapackInt; ab: ptr cfloat; ldab: LapackInt; ipiv: ptr LapackInt;
                anorm: cfloat; rcond: ptr cfloat; work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgbcon_work", dynlib: libName.}
proc sgbequ*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
            ku: LapackInt; ab: ptr cfloat; ldab: LapackInt; r: ptr cfloat; c: ptr cfloat;
            rowcnd: ptr cfloat; colcnd: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgbequ", dynlib: libName.}
proc sgbequWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
                ku: LapackInt; ab: ptr cfloat; ldab: LapackInt; r: ptr cfloat;
                c: ptr cfloat; rowcnd: ptr cfloat; colcnd: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgbequ_work", dynlib: libName.}
proc sgbequb*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
             ku: LapackInt; ab: ptr cfloat; ldab: LapackInt; r: ptr cfloat; c: ptr cfloat;
             rowcnd: ptr cfloat; colcnd: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgbequb", dynlib: libName.}
proc sgbequbWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
                 ku: LapackInt; ab: ptr cfloat; ldab: LapackInt; r: ptr cfloat;
                 c: ptr cfloat; rowcnd: ptr cfloat; colcnd: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgbequb_work", dynlib: libName.}
proc sgbrfs*(matrixLayout: cint; trans: char; n: LapackInt; kl: LapackInt; ku: LapackInt;
            nrhs: LapackInt; ab: ptr cfloat; ldab: LapackInt; afb: ptr cfloat;
            ldafb: LapackInt; ipiv: ptr LapackInt; b: ptr cfloat; ldb: LapackInt;
            x: ptr cfloat; ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgbrfs", dynlib: libName.}
proc sgbrfsWork*(matrixLayout: cint; trans: char; n: LapackInt; kl: LapackInt;
                ku: LapackInt; nrhs: LapackInt; ab: ptr cfloat; ldab: LapackInt;
                afb: ptr cfloat; ldafb: LapackInt; ipiv: ptr LapackInt; b: ptr cfloat;
                ldb: LapackInt; x: ptr cfloat; ldx: LapackInt; ferr: ptr cfloat;
                berr: ptr cfloat; work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgbrfs_work", dynlib: libName.}
proc sgbrfsx*(matrixLayout: cint; trans: char; equed: char; n: LapackInt; kl: LapackInt;
             ku: LapackInt; nrhs: LapackInt; ab: ptr cfloat; ldab: LapackInt;
             afb: ptr cfloat; ldafb: LapackInt; ipiv: ptr LapackInt; r: ptr cfloat;
             c: ptr cfloat; b: ptr cfloat; ldb: LapackInt; x: ptr cfloat; ldx: LapackInt;
             rcond: ptr cfloat; berr: ptr cfloat; nErrBnds: LapackInt;
             errBndsNorm: ptr cfloat; errBndsComp: ptr cfloat; nparams: LapackInt;
             params: ptr cfloat): LapackInt {.importc: "LAPACKE_sgbrfsx",
    dynlib: libName.}
proc sgbrfsxWork*(matrixLayout: cint; trans: char; equed: char; n: LapackInt;
                 kl: LapackInt; ku: LapackInt; nrhs: LapackInt; ab: ptr cfloat;
                 ldab: LapackInt; afb: ptr cfloat; ldafb: LapackInt;
                 ipiv: ptr LapackInt; r: ptr cfloat; c: ptr cfloat; b: ptr cfloat;
                 ldb: LapackInt; x: ptr cfloat; ldx: LapackInt; rcond: ptr cfloat;
                 berr: ptr cfloat; nErrBnds: LapackInt; errBndsNorm: ptr cfloat;
                 errBndsComp: ptr cfloat; nparams: LapackInt; params: ptr cfloat;
                 work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgbrfsx_work", dynlib: libName.}
proc sgbsv*(matrixLayout: cint; n: LapackInt; kl: LapackInt; ku: LapackInt;
           nrhs: LapackInt; ab: ptr cfloat; ldab: LapackInt; ipiv: ptr LapackInt;
           b: ptr cfloat; ldb: LapackInt): LapackInt {.importc: "LAPACKE_sgbsv",
    dynlib: libName.}
proc sgbsvWork*(matrixLayout: cint; n: LapackInt; kl: LapackInt; ku: LapackInt;
               nrhs: LapackInt; ab: ptr cfloat; ldab: LapackInt; ipiv: ptr LapackInt;
               b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_sgbsv_work", dynlib: libName.}
proc sgbsvx*(matrixLayout: cint; fact: char; trans: char; n: LapackInt; kl: LapackInt;
            ku: LapackInt; nrhs: LapackInt; ab: ptr cfloat; ldab: LapackInt;
            afb: ptr cfloat; ldafb: LapackInt; ipiv: ptr LapackInt; equed: cstring;
            r: ptr cfloat; c: ptr cfloat; b: ptr cfloat; ldb: LapackInt; x: ptr cfloat;
            ldx: LapackInt; rcond: ptr cfloat; ferr: ptr cfloat; berr: ptr cfloat;
            rpivot: ptr cfloat): LapackInt {.importc: "LAPACKE_sgbsvx",
    dynlib: libName.}
proc sgbsvxWork*(matrixLayout: cint; fact: char; trans: char; n: LapackInt;
                kl: LapackInt; ku: LapackInt; nrhs: LapackInt; ab: ptr cfloat;
                ldab: LapackInt; afb: ptr cfloat; ldafb: LapackInt;
                ipiv: ptr LapackInt; equed: cstring; r: ptr cfloat; c: ptr cfloat;
                b: ptr cfloat; ldb: LapackInt; x: ptr cfloat; ldx: LapackInt;
                rcond: ptr cfloat; ferr: ptr cfloat; berr: ptr cfloat; work: ptr cfloat;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_sgbsvx_work",
    dynlib: libName.}
proc sgbsvxx*(matrixLayout: cint; fact: char; trans: char; n: LapackInt; kl: LapackInt;
             ku: LapackInt; nrhs: LapackInt; ab: ptr cfloat; ldab: LapackInt;
             afb: ptr cfloat; ldafb: LapackInt; ipiv: ptr LapackInt; equed: cstring;
             r: ptr cfloat; c: ptr cfloat; b: ptr cfloat; ldb: LapackInt; x: ptr cfloat;
             ldx: LapackInt; rcond: ptr cfloat; rpvgrw: ptr cfloat; berr: ptr cfloat;
             nErrBnds: LapackInt; errBndsNorm: ptr cfloat; errBndsComp: ptr cfloat;
             nparams: LapackInt; params: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgbsvxx", dynlib: libName.}
proc sgbsvxxWork*(matrixLayout: cint; fact: char; trans: char; n: LapackInt;
                 kl: LapackInt; ku: LapackInt; nrhs: LapackInt; ab: ptr cfloat;
                 ldab: LapackInt; afb: ptr cfloat; ldafb: LapackInt;
                 ipiv: ptr LapackInt; equed: cstring; r: ptr cfloat; c: ptr cfloat;
                 b: ptr cfloat; ldb: LapackInt; x: ptr cfloat; ldx: LapackInt;
                 rcond: ptr cfloat; rpvgrw: ptr cfloat; berr: ptr cfloat;
                 nErrBnds: LapackInt; errBndsNorm: ptr cfloat;
                 errBndsComp: ptr cfloat; nparams: LapackInt; params: ptr cfloat;
                 work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgbsvxx_work", dynlib: libName.}
proc sgbtrf*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
            ku: LapackInt; ab: ptr cfloat; ldab: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgbtrf", dynlib: libName.}
proc sgbtrfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
                ku: LapackInt; ab: ptr cfloat; ldab: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgbtrf_work", dynlib: libName.}
proc sgbtrs*(matrixLayout: cint; trans: char; n: LapackInt; kl: LapackInt; ku: LapackInt;
            nrhs: LapackInt; ab: ptr cfloat; ldab: LapackInt; ipiv: ptr LapackInt;
            b: ptr cfloat; ldb: LapackInt): LapackInt {.importc: "LAPACKE_sgbtrs",
    dynlib: libName.}
proc sgbtrsWork*(matrixLayout: cint; trans: char; n: LapackInt; kl: LapackInt;
                ku: LapackInt; nrhs: LapackInt; ab: ptr cfloat; ldab: LapackInt;
                ipiv: ptr LapackInt; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_sgbtrs_work", dynlib: libName.}
proc sgebak*(matrixLayout: cint; job: char; side: char; n: LapackInt; ilo: LapackInt;
            ihi: LapackInt; scale: ptr cfloat; m: LapackInt; v: ptr cfloat; ldv: LapackInt): LapackInt {.
    importc: "LAPACKE_sgebak", dynlib: libName.}
proc sgebakWork*(matrixLayout: cint; job: char; side: char; n: LapackInt; ilo: LapackInt;
                ihi: LapackInt; scale: ptr cfloat; m: LapackInt; v: ptr cfloat;
                ldv: LapackInt): LapackInt {.importc: "LAPACKE_sgebak_work",
    dynlib: libName.}
proc sgebal*(matrixLayout: cint; job: char; n: LapackInt; a: ptr cfloat; lda: LapackInt;
            ilo: ptr LapackInt; ihi: ptr LapackInt; scale: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgebal", dynlib: libName.}
proc sgebalWork*(matrixLayout: cint; job: char; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; ilo: ptr LapackInt; ihi: ptr LapackInt;
                scale: ptr cfloat): LapackInt {.importc: "LAPACKE_sgebal_work",
    dynlib: libName.}
proc sgebrd*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
            lda: LapackInt; d: ptr cfloat; e: ptr cfloat; tauq: ptr cfloat;
            taup: ptr cfloat): LapackInt {.importc: "LAPACKE_sgebrd", dynlib: libName.}
proc sgebrdWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; d: ptr cfloat; e: ptr cfloat; tauq: ptr cfloat;
                taup: ptr cfloat; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sgebrd_work", dynlib: libName.}
proc sgecon*(matrixLayout: cint; norm: char; n: LapackInt; a: ptr cfloat; lda: LapackInt;
            anorm: cfloat; rcond: ptr cfloat): LapackInt {.importc: "LAPACKE_sgecon",
    dynlib: libName.}
proc sgeconWork*(matrixLayout: cint; norm: char; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; anorm: cfloat; rcond: ptr cfloat; work: ptr cfloat;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_sgecon_work",
    dynlib: libName.}
proc sgeequ*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
            lda: LapackInt; r: ptr cfloat; c: ptr cfloat; rowcnd: ptr cfloat;
            colcnd: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgeequ", dynlib: libName.}
proc sgeequWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; r: ptr cfloat; c: ptr cfloat; rowcnd: ptr cfloat;
                colcnd: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgeequ_work", dynlib: libName.}
proc sgeequb*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
             lda: LapackInt; r: ptr cfloat; c: ptr cfloat; rowcnd: ptr cfloat;
             colcnd: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgeequb", dynlib: libName.}
proc sgeequbWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
                 lda: LapackInt; r: ptr cfloat; c: ptr cfloat; rowcnd: ptr cfloat;
                 colcnd: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgeequb_work", dynlib: libName.}
proc sgees*(matrixLayout: cint; jobvs: char; sort: char; select: Lapack_S_Select2;
           n: LapackInt; a: ptr cfloat; lda: LapackInt; sdim: ptr LapackInt;
           wr: ptr cfloat; wi: ptr cfloat; vs: ptr cfloat; ldvs: LapackInt): LapackInt {.
    importc: "LAPACKE_sgees", dynlib: libName.}
proc sgeesWork*(matrixLayout: cint; jobvs: char; sort: char; select: Lapack_S_Select2;
               n: LapackInt; a: ptr cfloat; lda: LapackInt; sdim: ptr LapackInt;
               wr: ptr cfloat; wi: ptr cfloat; vs: ptr cfloat; ldvs: LapackInt;
               work: ptr cfloat; lwork: LapackInt; bwork: ptr LapackLogical): LapackInt {.
    importc: "LAPACKE_sgees_work", dynlib: libName.}
proc sgeesx*(matrixLayout: cint; jobvs: char; sort: char; select: Lapack_S_Select2;
            sense: char; n: LapackInt; a: ptr cfloat; lda: LapackInt;
            sdim: ptr LapackInt; wr: ptr cfloat; wi: ptr cfloat; vs: ptr cfloat;
            ldvs: LapackInt; rconde: ptr cfloat; rcondv: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgeesx", dynlib: libName.}
proc sgeesxWork*(matrixLayout: cint; jobvs: char; sort: char; select: Lapack_S_Select2;
                sense: char; n: LapackInt; a: ptr cfloat; lda: LapackInt;
                sdim: ptr LapackInt; wr: ptr cfloat; wi: ptr cfloat; vs: ptr cfloat;
                ldvs: LapackInt; rconde: ptr cfloat; rcondv: ptr cfloat;
                work: ptr cfloat; lwork: LapackInt; iwork: ptr LapackInt;
                liwork: LapackInt; bwork: ptr LapackLogical): LapackInt {.
    importc: "LAPACKE_sgeesx_work", dynlib: libName.}
proc sgeev*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt; a: ptr cfloat;
           lda: LapackInt; wr: ptr cfloat; wi: ptr cfloat; vl: ptr cfloat; ldvl: LapackInt;
           vr: ptr cfloat; ldvr: LapackInt): LapackInt {.importc: "LAPACKE_sgeev",
    dynlib: libName.}
proc sgeevWork*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt;
               a: ptr cfloat; lda: LapackInt; wr: ptr cfloat; wi: ptr cfloat;
               vl: ptr cfloat; ldvl: LapackInt; vr: ptr cfloat; ldvr: LapackInt;
               work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sgeev_work", dynlib: libName.}
proc sgeevx*(matrixLayout: cint; balanc: char; jobvl: char; jobvr: char; sense: char;
            n: LapackInt; a: ptr cfloat; lda: LapackInt; wr: ptr cfloat; wi: ptr cfloat;
            vl: ptr cfloat; ldvl: LapackInt; vr: ptr cfloat; ldvr: LapackInt;
            ilo: ptr LapackInt; ihi: ptr LapackInt; scale: ptr cfloat; abnrm: ptr cfloat;
            rconde: ptr cfloat; rcondv: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgeevx", dynlib: libName.}
proc sgeevxWork*(matrixLayout: cint; balanc: char; jobvl: char; jobvr: char; sense: char;
                n: LapackInt; a: ptr cfloat; lda: LapackInt; wr: ptr cfloat;
                wi: ptr cfloat; vl: ptr cfloat; ldvl: LapackInt; vr: ptr cfloat;
                ldvr: LapackInt; ilo: ptr LapackInt; ihi: ptr LapackInt;
                scale: ptr cfloat; abnrm: ptr cfloat; rconde: ptr cfloat;
                rcondv: ptr cfloat; work: ptr cfloat; lwork: LapackInt;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_sgeevx_work",
    dynlib: libName.}
proc sgehrd*(matrixLayout: cint; n: LapackInt; ilo: LapackInt; ihi: LapackInt;
            a: ptr cfloat; lda: LapackInt; tau: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgehrd", dynlib: libName.}
proc sgehrdWork*(matrixLayout: cint; n: LapackInt; ilo: LapackInt; ihi: LapackInt;
                a: ptr cfloat; lda: LapackInt; tau: ptr cfloat; work: ptr cfloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_sgehrd_work",
    dynlib: libName.}
proc sgejsv*(matrixLayout: cint; joba: char; jobu: char; jobv: char; jobr: char;
            jobt: char; jobp: char; m: LapackInt; n: LapackInt; a: ptr cfloat;
            lda: LapackInt; sva: ptr cfloat; u: ptr cfloat; ldu: LapackInt; v: ptr cfloat;
            ldv: LapackInt; stat: ptr cfloat; istat: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgejsv", dynlib: libName.}
proc sgejsvWork*(matrixLayout: cint; joba: char; jobu: char; jobv: char; jobr: char;
                jobt: char; jobp: char; m: LapackInt; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; sva: ptr cfloat; u: ptr cfloat; ldu: LapackInt;
                v: ptr cfloat; ldv: LapackInt; work: ptr cfloat; lwork: LapackInt;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_sgejsv_work",
    dynlib: libName.}
proc sgelq2*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
            lda: LapackInt; tau: ptr cfloat): LapackInt {.importc: "LAPACKE_sgelq2",
    dynlib: libName.}
proc sgelq2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; tau: ptr cfloat; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgelq2_work", dynlib: libName.}
proc sgelqf*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
            lda: LapackInt; tau: ptr cfloat): LapackInt {.importc: "LAPACKE_sgelqf",
    dynlib: libName.}
proc sgelqfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; tau: ptr cfloat; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sgelqf_work", dynlib: libName.}
proc sgels*(matrixLayout: cint; trans: char; m: LapackInt; n: LapackInt; nrhs: LapackInt;
           a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_sgels", dynlib: libName.}
proc sgelsWork*(matrixLayout: cint; trans: char; m: LapackInt; n: LapackInt;
               nrhs: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat;
               ldb: LapackInt; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sgels_work", dynlib: libName.}
proc sgelsd*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
            a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt; s: ptr cfloat;
            rcond: cfloat; rank: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgelsd", dynlib: libName.}
proc sgelsdWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
                a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
                s: ptr cfloat; rcond: cfloat; rank: ptr LapackInt; work: ptr cfloat;
                lwork: LapackInt; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgelsd_work", dynlib: libName.}
proc sgelss*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
            a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt; s: ptr cfloat;
            rcond: cfloat; rank: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgelss", dynlib: libName.}
proc sgelssWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
                a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
                s: ptr cfloat; rcond: cfloat; rank: ptr LapackInt; work: ptr cfloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_sgelss_work",
    dynlib: libName.}
proc sgelsy*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
            a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
            jpvt: ptr LapackInt; rcond: cfloat; rank: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgelsy", dynlib: libName.}
proc sgelsyWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
                a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
                jpvt: ptr LapackInt; rcond: cfloat; rank: ptr LapackInt;
                work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sgelsy_work", dynlib: libName.}
proc sgemqrt*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
             k: LapackInt; nb: LapackInt; v: ptr cfloat; ldv: LapackInt; t: ptr cfloat;
             ldt: LapackInt; c: ptr cfloat; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_sgemqrt", dynlib: libName.}
proc sgemqrtWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt;
                 n: LapackInt; k: LapackInt; nb: LapackInt; v: ptr cfloat;
                 ldv: LapackInt; t: ptr cfloat; ldt: LapackInt; c: ptr cfloat;
                 ldc: LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgemqrt_work", dynlib: libName.}
proc sgeqlf*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
            lda: LapackInt; tau: ptr cfloat): LapackInt {.importc: "LAPACKE_sgeqlf",
    dynlib: libName.}
proc sgeqlfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; tau: ptr cfloat; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sgeqlf_work", dynlib: libName.}
proc sgeqp3*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
            lda: LapackInt; jpvt: ptr LapackInt; tau: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgeqp3", dynlib: libName.}
proc sgeqp3Work*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; jpvt: ptr LapackInt; tau: ptr cfloat; work: ptr cfloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_sgeqp3_work",
    dynlib: libName.}
proc sgeqpf*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
            lda: LapackInt; jpvt: ptr LapackInt; tau: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgeqpf", dynlib: libName.}
proc sgeqpfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; jpvt: ptr LapackInt; tau: ptr cfloat; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgeqpf_work", dynlib: libName.}
proc sgeqr2*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
            lda: LapackInt; tau: ptr cfloat): LapackInt {.importc: "LAPACKE_sgeqr2",
    dynlib: libName.}
proc sgeqr2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; tau: ptr cfloat; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgeqr2_work", dynlib: libName.}
proc sgeqrf*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
            lda: LapackInt; tau: ptr cfloat): LapackInt {.importc: "LAPACKE_sgeqrf",
    dynlib: libName.}
proc sgeqrfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; tau: ptr cfloat; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sgeqrf_work", dynlib: libName.}
proc sgeqrfp*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
             lda: LapackInt; tau: ptr cfloat): LapackInt {.importc: "LAPACKE_sgeqrfp",
    dynlib: libName.}
proc sgeqrfpWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
                 lda: LapackInt; tau: ptr cfloat; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sgeqrfp_work", dynlib: libName.}
proc sgeqrt*(matrixLayout: cint; m: LapackInt; n: LapackInt; nb: LapackInt;
            a: ptr cfloat; lda: LapackInt; t: ptr cfloat; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_sgeqrt", dynlib: libName.}
proc sgeqrt2*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
             lda: LapackInt; t: ptr cfloat; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_sgeqrt2", dynlib: libName.}
proc sgeqrt2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
                 lda: LapackInt; t: ptr cfloat; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_sgeqrt2_work", dynlib: libName.}
proc sgeqrt3*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
             lda: LapackInt; t: ptr cfloat; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_sgeqrt3", dynlib: libName.}
proc sgeqrt3Work*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
                 lda: LapackInt; t: ptr cfloat; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_sgeqrt3_work", dynlib: libName.}
proc sgeqrtWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; nb: LapackInt;
                a: ptr cfloat; lda: LapackInt; t: ptr cfloat; ldt: LapackInt;
                work: ptr cfloat): LapackInt {.importc: "LAPACKE_sgeqrt_work",
    dynlib: libName.}
proc sgerfs*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
            a: ptr cfloat; lda: LapackInt; af: ptr cfloat; ldaf: LapackInt;
            ipiv: ptr LapackInt; b: ptr cfloat; ldb: LapackInt; x: ptr cfloat;
            ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgerfs", dynlib: libName.}
proc sgerfsWork*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
                a: ptr cfloat; lda: LapackInt; af: ptr cfloat; ldaf: LapackInt;
                ipiv: ptr LapackInt; b: ptr cfloat; ldb: LapackInt; x: ptr cfloat;
                ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat; work: ptr cfloat;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_sgerfs_work",
    dynlib: libName.}
proc sgerfsx*(matrixLayout: cint; trans: char; equed: char; n: LapackInt;
             nrhs: LapackInt; a: ptr cfloat; lda: LapackInt; af: ptr cfloat;
             ldaf: LapackInt; ipiv: ptr LapackInt; r: ptr cfloat; c: ptr cfloat;
             b: ptr cfloat; ldb: LapackInt; x: ptr cfloat; ldx: LapackInt;
             rcond: ptr cfloat; berr: ptr cfloat; nErrBnds: LapackInt;
             errBndsNorm: ptr cfloat; errBndsComp: ptr cfloat; nparams: LapackInt;
             params: ptr cfloat): LapackInt {.importc: "LAPACKE_sgerfsx",
    dynlib: libName.}
proc sgerfsxWork*(matrixLayout: cint; trans: char; equed: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr cfloat; lda: LapackInt; af: ptr cfloat;
                 ldaf: LapackInt; ipiv: ptr LapackInt; r: ptr cfloat; c: ptr cfloat;
                 b: ptr cfloat; ldb: LapackInt; x: ptr cfloat; ldx: LapackInt;
                 rcond: ptr cfloat; berr: ptr cfloat; nErrBnds: LapackInt;
                 errBndsNorm: ptr cfloat; errBndsComp: ptr cfloat; nparams: LapackInt;
                 params: ptr cfloat; work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgerfsx_work", dynlib: libName.}
proc sgerqf*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
            lda: LapackInt; tau: ptr cfloat): LapackInt {.importc: "LAPACKE_sgerqf",
    dynlib: libName.}
proc sgerqfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; tau: ptr cfloat; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sgerqf_work", dynlib: libName.}
proc sgesdd*(matrixLayout: cint; jobz: char; m: LapackInt; n: LapackInt; a: ptr cfloat;
            lda: LapackInt; s: ptr cfloat; u: ptr cfloat; ldu: LapackInt; vt: ptr cfloat;
            ldvt: LapackInt): LapackInt {.importc: "LAPACKE_sgesdd", dynlib: libName.}
proc sgesddWork*(matrixLayout: cint; jobz: char; m: LapackInt; n: LapackInt;
                a: ptr cfloat; lda: LapackInt; s: ptr cfloat; u: ptr cfloat;
                ldu: LapackInt; vt: ptr cfloat; ldvt: LapackInt; work: ptr cfloat;
                lwork: LapackInt; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgesdd_work", dynlib: libName.}
proc sgesv*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; a: ptr cfloat;
           lda: LapackInt; ipiv: ptr LapackInt; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_sgesv", dynlib: libName.}
proc sgesvWork*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; a: ptr cfloat;
               lda: LapackInt; ipiv: ptr LapackInt; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_sgesv_work", dynlib: libName.}
proc sgesvd*(matrixLayout: cint; jobu: char; jobvt: char; m: LapackInt; n: LapackInt;
            a: ptr cfloat; lda: LapackInt; s: ptr cfloat; u: ptr cfloat; ldu: LapackInt;
            vt: ptr cfloat; ldvt: LapackInt; superb: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgesvd", dynlib: libName.}
proc sgesvdWork*(matrixLayout: cint; jobu: char; jobvt: char; m: LapackInt; n: LapackInt;
                a: ptr cfloat; lda: LapackInt; s: ptr cfloat; u: ptr cfloat;
                ldu: LapackInt; vt: ptr cfloat; ldvt: LapackInt; work: ptr cfloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_sgesvd_work",
    dynlib: libName.}
proc sgesvdx*(matrixLayout: cint; jobu: char; jobvt: char; range: char; m: LapackInt;
             n: LapackInt; a: ptr cfloat; lda: LapackInt; vl: cfloat; vu: cfloat;
             il: LapackInt; iu: LapackInt; ns: ptr LapackInt; s: ptr cfloat;
             u: ptr cfloat; ldu: LapackInt; vt: ptr cfloat; ldvt: LapackInt;
             superb: ptr LapackInt): LapackInt {.importc: "LAPACKE_sgesvdx",
    dynlib: libName.}
proc sgesvdxWork*(matrixLayout: cint; jobu: char; jobvt: char; range: char; m: LapackInt;
                 n: LapackInt; a: ptr cfloat; lda: LapackInt; vl: cfloat; vu: cfloat;
                 il: LapackInt; iu: LapackInt; ns: ptr LapackInt; s: ptr cfloat;
                 u: ptr cfloat; ldu: LapackInt; vt: ptr cfloat; ldvt: LapackInt;
                 work: ptr cfloat; lwork: LapackInt; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgesvdx_work", dynlib: libName.}
proc sgesvj*(matrixLayout: cint; joba: char; jobu: char; jobv: char; m: LapackInt;
            n: LapackInt; a: ptr cfloat; lda: LapackInt; sva: ptr cfloat; mv: LapackInt;
            v: ptr cfloat; ldv: LapackInt; stat: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgesvj", dynlib: libName.}
proc sgesvjWork*(matrixLayout: cint; joba: char; jobu: char; jobv: char; m: LapackInt;
                n: LapackInt; a: ptr cfloat; lda: LapackInt; sva: ptr cfloat;
                mv: LapackInt; v: ptr cfloat; ldv: LapackInt; work: ptr cfloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_sgesvj_work",
    dynlib: libName.}
proc sgesvx*(matrixLayout: cint; fact: char; trans: char; n: LapackInt; nrhs: LapackInt;
            a: ptr cfloat; lda: LapackInt; af: ptr cfloat; ldaf: LapackInt;
            ipiv: ptr LapackInt; equed: cstring; r: ptr cfloat; c: ptr cfloat;
            b: ptr cfloat; ldb: LapackInt; x: ptr cfloat; ldx: LapackInt;
            rcond: ptr cfloat; ferr: ptr cfloat; berr: ptr cfloat; rpivot: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgesvx", dynlib: libName.}
proc sgesvxWork*(matrixLayout: cint; fact: char; trans: char; n: LapackInt;
                nrhs: LapackInt; a: ptr cfloat; lda: LapackInt; af: ptr cfloat;
                ldaf: LapackInt; ipiv: ptr LapackInt; equed: cstring; r: ptr cfloat;
                c: ptr cfloat; b: ptr cfloat; ldb: LapackInt; x: ptr cfloat;
                ldx: LapackInt; rcond: ptr cfloat; ferr: ptr cfloat; berr: ptr cfloat;
                work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgesvx_work", dynlib: libName.}
proc sgesvxx*(matrixLayout: cint; fact: char; trans: char; n: LapackInt; nrhs: LapackInt;
             a: ptr cfloat; lda: LapackInt; af: ptr cfloat; ldaf: LapackInt;
             ipiv: ptr LapackInt; equed: cstring; r: ptr cfloat; c: ptr cfloat;
             b: ptr cfloat; ldb: LapackInt; x: ptr cfloat; ldx: LapackInt;
             rcond: ptr cfloat; rpvgrw: ptr cfloat; berr: ptr cfloat;
             nErrBnds: LapackInt; errBndsNorm: ptr cfloat; errBndsComp: ptr cfloat;
             nparams: LapackInt; params: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgesvxx", dynlib: libName.}
proc sgesvxxWork*(matrixLayout: cint; fact: char; trans: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr cfloat; lda: LapackInt; af: ptr cfloat;
                 ldaf: LapackInt; ipiv: ptr LapackInt; equed: cstring; r: ptr cfloat;
                 c: ptr cfloat; b: ptr cfloat; ldb: LapackInt; x: ptr cfloat;
                 ldx: LapackInt; rcond: ptr cfloat; rpvgrw: ptr cfloat;
                 berr: ptr cfloat; nErrBnds: LapackInt; errBndsNorm: ptr cfloat;
                 errBndsComp: ptr cfloat; nparams: LapackInt; params: ptr cfloat;
                 work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgesvxx_work", dynlib: libName.}
proc sgetf2*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
            lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgetf2", dynlib: libName.}
proc sgetf2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgetf2_work", dynlib: libName.}
proc sgetrf*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
            lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgetrf", dynlib: libName.}
proc sgetrf2*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
             lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgetrf2", dynlib: libName.}
proc sgetrf2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
                 lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgetrf2_work", dynlib: libName.}
proc sgetrfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgetrf_work", dynlib: libName.}
proc sgetri*(matrixLayout: cint; n: LapackInt; a: ptr cfloat; lda: LapackInt;
            ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_sgetri",
    dynlib: libName.}
proc sgetriWork*(matrixLayout: cint; n: LapackInt; a: ptr cfloat; lda: LapackInt;
                ipiv: ptr LapackInt; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sgetri_work", dynlib: libName.}
proc sgetrs*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
            a: ptr cfloat; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cfloat;
            ldb: LapackInt): LapackInt {.importc: "LAPACKE_sgetrs", dynlib: libName.}
proc sgetrsWork*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
                a: ptr cfloat; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cfloat;
                ldb: LapackInt): LapackInt {.importc: "LAPACKE_sgetrs_work",
    dynlib: libName.}
proc sggbak*(matrixLayout: cint; job: char; side: char; n: LapackInt; ilo: LapackInt;
            ihi: LapackInt; lscale: ptr cfloat; rscale: ptr cfloat; m: LapackInt;
            v: ptr cfloat; ldv: LapackInt): LapackInt {.importc: "LAPACKE_sggbak",
    dynlib: libName.}
proc sggbakWork*(matrixLayout: cint; job: char; side: char; n: LapackInt; ilo: LapackInt;
                ihi: LapackInt; lscale: ptr cfloat; rscale: ptr cfloat; m: LapackInt;
                v: ptr cfloat; ldv: LapackInt): LapackInt {.
    importc: "LAPACKE_sggbak_work", dynlib: libName.}
proc sggbal*(matrixLayout: cint; job: char; n: LapackInt; a: ptr cfloat; lda: LapackInt;
            b: ptr cfloat; ldb: LapackInt; ilo: ptr LapackInt; ihi: ptr LapackInt;
            lscale: ptr cfloat; rscale: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sggbal", dynlib: libName.}
proc sggbalWork*(matrixLayout: cint; job: char; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; b: ptr cfloat; ldb: LapackInt; ilo: ptr LapackInt;
                ihi: ptr LapackInt; lscale: ptr cfloat; rscale: ptr cfloat;
                work: ptr cfloat): LapackInt {.importc: "LAPACKE_sggbal_work",
    dynlib: libName.}
proc sgges*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
           selctg: Lapack_S_Select3; n: LapackInt; a: ptr cfloat; lda: LapackInt;
           b: ptr cfloat; ldb: LapackInt; sdim: ptr LapackInt; alphar: ptr cfloat;
           alphai: ptr cfloat; beta: ptr cfloat; vsl: ptr cfloat; ldvsl: LapackInt;
           vsr: ptr cfloat; ldvsr: LapackInt): LapackInt {.importc: "LAPACKE_sgges",
    dynlib: libName.}
proc sgges3*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
            selctg: Lapack_S_Select3; n: LapackInt; a: ptr cfloat; lda: LapackInt;
            b: ptr cfloat; ldb: LapackInt; sdim: ptr LapackInt; alphar: ptr cfloat;
            alphai: ptr cfloat; beta: ptr cfloat; vsl: ptr cfloat; ldvsl: LapackInt;
            vsr: ptr cfloat; ldvsr: LapackInt): LapackInt {.importc: "LAPACKE_sgges3",
    dynlib: libName.}
proc sgges3Work*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
                selctg: Lapack_S_Select3; n: LapackInt; a: ptr cfloat; lda: LapackInt;
                b: ptr cfloat; ldb: LapackInt; sdim: ptr LapackInt; alphar: ptr cfloat;
                alphai: ptr cfloat; beta: ptr cfloat; vsl: ptr cfloat; ldvsl: LapackInt;
                vsr: ptr cfloat; ldvsr: LapackInt; work: ptr cfloat; lwork: LapackInt;
                bwork: ptr LapackLogical): LapackInt {.
    importc: "LAPACKE_sgges3_work", dynlib: libName.}
proc sggesWork*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
               selctg: Lapack_S_Select3; n: LapackInt; a: ptr cfloat; lda: LapackInt;
               b: ptr cfloat; ldb: LapackInt; sdim: ptr LapackInt; alphar: ptr cfloat;
               alphai: ptr cfloat; beta: ptr cfloat; vsl: ptr cfloat; ldvsl: LapackInt;
               vsr: ptr cfloat; ldvsr: LapackInt; work: ptr cfloat; lwork: LapackInt;
               bwork: ptr LapackLogical): LapackInt {.importc: "LAPACKE_sgges_work",
    dynlib: libName.}
proc sggesx*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
            selctg: Lapack_S_Select3; sense: char; n: LapackInt; a: ptr cfloat;
            lda: LapackInt; b: ptr cfloat; ldb: LapackInt; sdim: ptr LapackInt;
            alphar: ptr cfloat; alphai: ptr cfloat; beta: ptr cfloat; vsl: ptr cfloat;
            ldvsl: LapackInt; vsr: ptr cfloat; ldvsr: LapackInt; rconde: ptr cfloat;
            rcondv: ptr cfloat): LapackInt {.importc: "LAPACKE_sggesx",
    dynlib: libName.}
proc sggesxWork*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
                selctg: Lapack_S_Select3; sense: char; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; b: ptr cfloat; ldb: LapackInt; sdim: ptr LapackInt;
                alphar: ptr cfloat; alphai: ptr cfloat; beta: ptr cfloat;
                vsl: ptr cfloat; ldvsl: LapackInt; vsr: ptr cfloat; ldvsr: LapackInt;
                rconde: ptr cfloat; rcondv: ptr cfloat; work: ptr cfloat;
                lwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt;
                bwork: ptr LapackLogical): LapackInt {.
    importc: "LAPACKE_sggesx_work", dynlib: libName.}
proc sggev*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt; a: ptr cfloat;
           lda: LapackInt; b: ptr cfloat; ldb: LapackInt; alphar: ptr cfloat;
           alphai: ptr cfloat; beta: ptr cfloat; vl: ptr cfloat; ldvl: LapackInt;
           vr: ptr cfloat; ldvr: LapackInt): LapackInt {.importc: "LAPACKE_sggev",
    dynlib: libName.}
proc sggev3*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt; a: ptr cfloat;
            lda: LapackInt; b: ptr cfloat; ldb: LapackInt; alphar: ptr cfloat;
            alphai: ptr cfloat; beta: ptr cfloat; vl: ptr cfloat; ldvl: LapackInt;
            vr: ptr cfloat; ldvr: LapackInt): LapackInt {.importc: "LAPACKE_sggev3",
    dynlib: libName.}
proc sggev3Work*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt;
                a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
                alphar: ptr cfloat; alphai: ptr cfloat; beta: ptr cfloat; vl: ptr cfloat;
                ldvl: LapackInt; vr: ptr cfloat; ldvr: LapackInt; work: ptr cfloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_sggev3_work",
    dynlib: libName.}
proc sggevWork*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt;
               a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
               alphar: ptr cfloat; alphai: ptr cfloat; beta: ptr cfloat; vl: ptr cfloat;
               ldvl: LapackInt; vr: ptr cfloat; ldvr: LapackInt; work: ptr cfloat;
               lwork: LapackInt): LapackInt {.importc: "LAPACKE_sggev_work",
    dynlib: libName.}
proc sggevx*(matrixLayout: cint; balanc: char; jobvl: char; jobvr: char; sense: char;
            n: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
            alphar: ptr cfloat; alphai: ptr cfloat; beta: ptr cfloat; vl: ptr cfloat;
            ldvl: LapackInt; vr: ptr cfloat; ldvr: LapackInt; ilo: ptr LapackInt;
            ihi: ptr LapackInt; lscale: ptr cfloat; rscale: ptr cfloat;
            abnrm: ptr cfloat; bbnrm: ptr cfloat; rconde: ptr cfloat; rcondv: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sggevx", dynlib: libName.}
proc sggevxWork*(matrixLayout: cint; balanc: char; jobvl: char; jobvr: char; sense: char;
                n: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat;
                ldb: LapackInt; alphar: ptr cfloat; alphai: ptr cfloat;
                beta: ptr cfloat; vl: ptr cfloat; ldvl: LapackInt; vr: ptr cfloat;
                ldvr: LapackInt; ilo: ptr LapackInt; ihi: ptr LapackInt;
                lscale: ptr cfloat; rscale: ptr cfloat; abnrm: ptr cfloat;
                bbnrm: ptr cfloat; rconde: ptr cfloat; rcondv: ptr cfloat;
                work: ptr cfloat; lwork: LapackInt; iwork: ptr LapackInt;
                bwork: ptr LapackLogical): LapackInt {.
    importc: "LAPACKE_sggevx_work", dynlib: libName.}
proc sggglm*(matrixLayout: cint; n: LapackInt; m: LapackInt; p: LapackInt; a: ptr cfloat;
            lda: LapackInt; b: ptr cfloat; ldb: LapackInt; d: ptr cfloat; x: ptr cfloat;
            y: ptr cfloat): LapackInt {.importc: "LAPACKE_sggglm", dynlib: libName.}
proc sggglmWork*(matrixLayout: cint; n: LapackInt; m: LapackInt; p: LapackInt;
                a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
                d: ptr cfloat; x: ptr cfloat; y: ptr cfloat; work: ptr cfloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_sggglm_work",
    dynlib: libName.}
proc sgghd3*(matrixLayout: cint; compq: char; compz: char; n: LapackInt; ilo: LapackInt;
            ihi: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
            q: ptr cfloat; ldq: LapackInt; z: ptr cfloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_sgghd3", dynlib: libName.}
proc sgghd3Work*(matrixLayout: cint; compq: char; compz: char; n: LapackInt;
                ilo: LapackInt; ihi: LapackInt; a: ptr cfloat; lda: LapackInt;
                b: ptr cfloat; ldb: LapackInt; q: ptr cfloat; ldq: LapackInt;
                z: ptr cfloat; ldz: LapackInt; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sgghd3_work", dynlib: libName.}
proc sgghrd*(matrixLayout: cint; compq: char; compz: char; n: LapackInt; ilo: LapackInt;
            ihi: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
            q: ptr cfloat; ldq: LapackInt; z: ptr cfloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_sgghrd", dynlib: libName.}
proc sgghrdWork*(matrixLayout: cint; compq: char; compz: char; n: LapackInt;
                ilo: LapackInt; ihi: LapackInt; a: ptr cfloat; lda: LapackInt;
                b: ptr cfloat; ldb: LapackInt; q: ptr cfloat; ldq: LapackInt;
                z: ptr cfloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_sgghrd_work", dynlib: libName.}
proc sgglse*(matrixLayout: cint; m: LapackInt; n: LapackInt; p: LapackInt; a: ptr cfloat;
            lda: LapackInt; b: ptr cfloat; ldb: LapackInt; c: ptr cfloat; d: ptr cfloat;
            x: ptr cfloat): LapackInt {.importc: "LAPACKE_sgglse", dynlib: libName.}
proc sgglseWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; p: LapackInt;
                a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
                c: ptr cfloat; d: ptr cfloat; x: ptr cfloat; work: ptr cfloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_sgglse_work",
    dynlib: libName.}
proc sggqrf*(matrixLayout: cint; n: LapackInt; m: LapackInt; p: LapackInt; a: ptr cfloat;
            lda: LapackInt; taua: ptr cfloat; b: ptr cfloat; ldb: LapackInt;
            taub: ptr cfloat): LapackInt {.importc: "LAPACKE_sggqrf", dynlib: libName.}
proc sggqrfWork*(matrixLayout: cint; n: LapackInt; m: LapackInt; p: LapackInt;
                a: ptr cfloat; lda: LapackInt; taua: ptr cfloat; b: ptr cfloat;
                ldb: LapackInt; taub: ptr cfloat; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sggqrf_work", dynlib: libName.}
proc sggrqf*(matrixLayout: cint; m: LapackInt; p: LapackInt; n: LapackInt; a: ptr cfloat;
            lda: LapackInt; taua: ptr cfloat; b: ptr cfloat; ldb: LapackInt;
            taub: ptr cfloat): LapackInt {.importc: "LAPACKE_sggrqf", dynlib: libName.}
proc sggrqfWork*(matrixLayout: cint; m: LapackInt; p: LapackInt; n: LapackInt;
                a: ptr cfloat; lda: LapackInt; taua: ptr cfloat; b: ptr cfloat;
                ldb: LapackInt; taub: ptr cfloat; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sggrqf_work", dynlib: libName.}
proc sggsvd*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
            n: LapackInt; p: LapackInt; k: ptr LapackInt; l: ptr LapackInt; a: ptr cfloat;
            lda: LapackInt; b: ptr cfloat; ldb: LapackInt; alpha: ptr cfloat;
            beta: ptr cfloat; u: ptr cfloat; ldu: LapackInt; v: ptr cfloat; ldv: LapackInt;
            q: ptr cfloat; ldq: LapackInt; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sggsvd", dynlib: libName.}
proc sggsvd3*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
             n: LapackInt; p: LapackInt; k: ptr LapackInt; l: ptr LapackInt;
             a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
             alpha: ptr cfloat; beta: ptr cfloat; u: ptr cfloat; ldu: LapackInt;
             v: ptr cfloat; ldv: LapackInt; q: ptr cfloat; ldq: LapackInt;
             iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_sggsvd3",
    dynlib: libName.}
proc sggsvd3Work*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
                 n: LapackInt; p: LapackInt; k: ptr LapackInt; l: ptr LapackInt;
                 a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
                 alpha: ptr cfloat; beta: ptr cfloat; u: ptr cfloat; ldu: LapackInt;
                 v: ptr cfloat; ldv: LapackInt; q: ptr cfloat; ldq: LapackInt;
                 work: ptr cfloat; lwork: LapackInt; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sggsvd3_work", dynlib: libName.}
proc sggsvdWork*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
                n: LapackInt; p: LapackInt; k: ptr LapackInt; l: ptr LapackInt;
                a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
                alpha: ptr cfloat; beta: ptr cfloat; u: ptr cfloat; ldu: LapackInt;
                v: ptr cfloat; ldv: LapackInt; q: ptr cfloat; ldq: LapackInt;
                work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sggsvd_work", dynlib: libName.}
proc sggsvp*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
            p: LapackInt; n: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat;
            ldb: LapackInt; tola: cfloat; tolb: cfloat; k: ptr LapackInt;
            l: ptr LapackInt; u: ptr cfloat; ldu: LapackInt; v: ptr cfloat; ldv: LapackInt;
            q: ptr cfloat; ldq: LapackInt): LapackInt {.importc: "LAPACKE_sggsvp",
    dynlib: libName.}
proc sggsvp3*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
             p: LapackInt; n: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat;
             ldb: LapackInt; tola: cfloat; tolb: cfloat; k: ptr LapackInt;
             l: ptr LapackInt; u: ptr cfloat; ldu: LapackInt; v: ptr cfloat;
             ldv: LapackInt; q: ptr cfloat; ldq: LapackInt): LapackInt {.
    importc: "LAPACKE_sggsvp3", dynlib: libName.}
proc sggsvp3Work*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
                 p: LapackInt; n: LapackInt; a: ptr cfloat; lda: LapackInt;
                 b: ptr cfloat; ldb: LapackInt; tola: cfloat; tolb: cfloat;
                 k: ptr LapackInt; l: ptr LapackInt; u: ptr cfloat; ldu: LapackInt;
                 v: ptr cfloat; ldv: LapackInt; q: ptr cfloat; ldq: LapackInt;
                 iwork: ptr LapackInt; tau: ptr cfloat; work: ptr cfloat;
                 lwork: LapackInt): LapackInt {.importc: "LAPACKE_sggsvp3_work",
    dynlib: libName.}
proc sggsvpWork*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
                p: LapackInt; n: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat;
                ldb: LapackInt; tola: cfloat; tolb: cfloat; k: ptr LapackInt;
                l: ptr LapackInt; u: ptr cfloat; ldu: LapackInt; v: ptr cfloat;
                ldv: LapackInt; q: ptr cfloat; ldq: LapackInt; iwork: ptr LapackInt;
                tau: ptr cfloat; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sggsvp_work", dynlib: libName.}
proc sgtcon*(norm: char; n: LapackInt; dl: ptr cfloat; d: ptr cfloat; du: ptr cfloat;
            du2: ptr cfloat; ipiv: ptr LapackInt; anorm: cfloat; rcond: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sgtcon", dynlib: libName.}
proc sgtconWork*(norm: char; n: LapackInt; dl: ptr cfloat; d: ptr cfloat; du: ptr cfloat;
                du2: ptr cfloat; ipiv: ptr LapackInt; anorm: cfloat; rcond: ptr cfloat;
                work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgtcon_work", dynlib: libName.}
proc sgtrfs*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
            dl: ptr cfloat; d: ptr cfloat; du: ptr cfloat; dlf: ptr cfloat; df: ptr cfloat;
            duf: ptr cfloat; du2: ptr cfloat; ipiv: ptr LapackInt; b: ptr cfloat;
            ldb: LapackInt; x: ptr cfloat; ldx: LapackInt; ferr: ptr cfloat;
            berr: ptr cfloat): LapackInt {.importc: "LAPACKE_sgtrfs", dynlib: libName.}
proc sgtrfsWork*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
                dl: ptr cfloat; d: ptr cfloat; du: ptr cfloat; dlf: ptr cfloat;
                df: ptr cfloat; duf: ptr cfloat; du2: ptr cfloat; ipiv: ptr LapackInt;
                b: ptr cfloat; ldb: LapackInt; x: ptr cfloat; ldx: LapackInt;
                ferr: ptr cfloat; berr: ptr cfloat; work: ptr cfloat;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_sgtrfs_work",
    dynlib: libName.}
proc sgtsv*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; dl: ptr cfloat;
           d: ptr cfloat; du: ptr cfloat; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_sgtsv", dynlib: libName.}
proc sgtsvWork*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; dl: ptr cfloat;
               d: ptr cfloat; du: ptr cfloat; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_sgtsv_work", dynlib: libName.}
proc sgtsvx*(matrixLayout: cint; fact: char; trans: char; n: LapackInt; nrhs: LapackInt;
            dl: ptr cfloat; d: ptr cfloat; du: ptr cfloat; dlf: ptr cfloat; df: ptr cfloat;
            duf: ptr cfloat; du2: ptr cfloat; ipiv: ptr LapackInt; b: ptr cfloat;
            ldb: LapackInt; x: ptr cfloat; ldx: LapackInt; rcond: ptr cfloat;
            ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.importc: "LAPACKE_sgtsvx",
    dynlib: libName.}
proc sgtsvxWork*(matrixLayout: cint; fact: char; trans: char; n: LapackInt;
                nrhs: LapackInt; dl: ptr cfloat; d: ptr cfloat; du: ptr cfloat;
                dlf: ptr cfloat; df: ptr cfloat; duf: ptr cfloat; du2: ptr cfloat;
                ipiv: ptr LapackInt; b: ptr cfloat; ldb: LapackInt; x: ptr cfloat;
                ldx: LapackInt; rcond: ptr cfloat; ferr: ptr cfloat; berr: ptr cfloat;
                work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgtsvx_work", dynlib: libName.}
proc sgttrf*(n: LapackInt; dl: ptr cfloat; d: ptr cfloat; du: ptr cfloat; du2: ptr cfloat;
            ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_sgttrf",
    dynlib: libName.}
proc sgttrfWork*(n: LapackInt; dl: ptr cfloat; d: ptr cfloat; du: ptr cfloat;
                du2: ptr cfloat; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sgttrf_work", dynlib: libName.}
proc sgttrs*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
            dl: ptr cfloat; d: ptr cfloat; du: ptr cfloat; du2: ptr cfloat;
            ipiv: ptr LapackInt; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_sgttrs", dynlib: libName.}
proc sgttrsWork*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
                dl: ptr cfloat; d: ptr cfloat; du: ptr cfloat; du2: ptr cfloat;
                ipiv: ptr LapackInt; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_sgttrs_work", dynlib: libName.}
proc shgeqz*(matrixLayout: cint; job: char; compq: char; compz: char; n: LapackInt;
            ilo: LapackInt; ihi: LapackInt; h: ptr cfloat; ldh: LapackInt; t: ptr cfloat;
            ldt: LapackInt; alphar: ptr cfloat; alphai: ptr cfloat; beta: ptr cfloat;
            q: ptr cfloat; ldq: LapackInt; z: ptr cfloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_shgeqz", dynlib: libName.}
proc shgeqzWork*(matrixLayout: cint; job: char; compq: char; compz: char; n: LapackInt;
                ilo: LapackInt; ihi: LapackInt; h: ptr cfloat; ldh: LapackInt;
                t: ptr cfloat; ldt: LapackInt; alphar: ptr cfloat; alphai: ptr cfloat;
                beta: ptr cfloat; q: ptr cfloat; ldq: LapackInt; z: ptr cfloat;
                ldz: LapackInt; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_shgeqz_work", dynlib: libName.}
proc shsein*(matrixLayout: cint; job: char; eigsrc: char; initv: char;
            select: ptr LapackLogical; n: LapackInt; h: ptr cfloat; ldh: LapackInt;
            wr: ptr cfloat; wi: ptr cfloat; vl: ptr cfloat; ldvl: LapackInt;
            vr: ptr cfloat; ldvr: LapackInt; mm: LapackInt; m: ptr LapackInt;
            ifaill: ptr LapackInt; ifailr: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_shsein", dynlib: libName.}
proc shseinWork*(matrixLayout: cint; job: char; eigsrc: char; initv: char;
                select: ptr LapackLogical; n: LapackInt; h: ptr cfloat; ldh: LapackInt;
                wr: ptr cfloat; wi: ptr cfloat; vl: ptr cfloat; ldvl: LapackInt;
                vr: ptr cfloat; ldvr: LapackInt; mm: LapackInt; m: ptr LapackInt;
                work: ptr cfloat; ifaill: ptr LapackInt; ifailr: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_shsein_work", dynlib: libName.}
proc shseqr*(matrixLayout: cint; job: char; compz: char; n: LapackInt; ilo: LapackInt;
            ihi: LapackInt; h: ptr cfloat; ldh: LapackInt; wr: ptr cfloat; wi: ptr cfloat;
            z: ptr cfloat; ldz: LapackInt): LapackInt {.importc: "LAPACKE_shseqr",
    dynlib: libName.}
proc shseqrWork*(matrixLayout: cint; job: char; compz: char; n: LapackInt;
                ilo: LapackInt; ihi: LapackInt; h: ptr cfloat; ldh: LapackInt;
                wr: ptr cfloat; wi: ptr cfloat; z: ptr cfloat; ldz: LapackInt;
                work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_shseqr_work", dynlib: libName.}
proc slacn2*(n: LapackInt; v: ptr cfloat; x: ptr cfloat; isgn: ptr LapackInt;
            est: ptr cfloat; kase: ptr LapackInt; isave: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_slacn2", dynlib: libName.}
proc slacn2Work*(n: LapackInt; v: ptr cfloat; x: ptr cfloat; isgn: ptr LapackInt;
                est: ptr cfloat; kase: ptr LapackInt; isave: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_slacn2_work", dynlib: libName.}
proc slacpy*(matrixLayout: cint; uplo: char; m: LapackInt; n: LapackInt; a: ptr cfloat;
            lda: LapackInt; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_slacpy", dynlib: libName.}
proc slacpyWork*(matrixLayout: cint; uplo: char; m: LapackInt; n: LapackInt;
                a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_slacpy_work", dynlib: libName.}
proc slag2d*(matrixLayout: cint; m: LapackInt; n: LapackInt; sa: ptr cfloat;
            ldsa: LapackInt; a: ptr cdouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_slag2d", dynlib: libName.}
proc slag2dWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; sa: ptr cfloat;
                ldsa: LapackInt; a: ptr cdouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_slag2d_work", dynlib: libName.}
proc slagge*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
            ku: LapackInt; d: ptr cfloat; a: ptr cfloat; lda: LapackInt;
            iseed: ptr LapackInt): LapackInt {.importc: "LAPACKE_slagge",
    dynlib: libName.}
proc slaggeWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
                ku: LapackInt; d: ptr cfloat; a: ptr cfloat; lda: LapackInt;
                iseed: ptr LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_slagge_work", dynlib: libName.}
proc slagsy*(matrixLayout: cint; n: LapackInt; k: LapackInt; d: ptr cfloat; a: ptr cfloat;
            lda: LapackInt; iseed: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_slagsy", dynlib: libName.}
proc slagsyWork*(matrixLayout: cint; n: LapackInt; k: LapackInt; d: ptr cfloat;
                a: ptr cfloat; lda: LapackInt; iseed: ptr LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_slagsy_work", dynlib: libName.}
proc slapmr*(matrixLayout: cint; forwrd: LapackLogical; m: LapackInt; n: LapackInt;
            x: ptr cfloat; ldx: LapackInt; k: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_slapmr", dynlib: libName.}
proc slapmrWork*(matrixLayout: cint; forwrd: LapackLogical; m: LapackInt; n: LapackInt;
                x: ptr cfloat; ldx: LapackInt; k: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_slapmr_work", dynlib: libName.}
proc slapmt*(matrixLayout: cint; forwrd: LapackLogical; m: LapackInt; n: LapackInt;
            x: ptr cfloat; ldx: LapackInt; k: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_slapmt", dynlib: libName.}
proc slapmtWork*(matrixLayout: cint; forwrd: LapackLogical; m: LapackInt; n: LapackInt;
                x: ptr cfloat; ldx: LapackInt; k: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_slapmt_work", dynlib: libName.}
proc slarfb*(matrixLayout: cint; side: char; trans: char; direct: char; storev: char;
            m: LapackInt; n: LapackInt; k: LapackInt; v: ptr cfloat; ldv: LapackInt;
            t: ptr cfloat; ldt: LapackInt; c: ptr cfloat; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_slarfb", dynlib: libName.}
proc slarfbWork*(matrixLayout: cint; side: char; trans: char; direct: char; storev: char;
                m: LapackInt; n: LapackInt; k: LapackInt; v: ptr cfloat; ldv: LapackInt;
                t: ptr cfloat; ldt: LapackInt; c: ptr cfloat; ldc: LapackInt;
                work: ptr cfloat; ldwork: LapackInt): LapackInt {.
    importc: "LAPACKE_slarfb_work", dynlib: libName.}
proc slarfg*(n: LapackInt; alpha: ptr cfloat; x: ptr cfloat; incx: LapackInt;
            tau: ptr cfloat): LapackInt {.importc: "LAPACKE_slarfg", dynlib: libName.}
proc slarfgWork*(n: LapackInt; alpha: ptr cfloat; x: ptr cfloat; incx: LapackInt;
                tau: ptr cfloat): LapackInt {.importc: "LAPACKE_slarfg_work",
    dynlib: libName.}
proc slarft*(matrixLayout: cint; direct: char; storev: char; n: LapackInt; k: LapackInt;
            v: ptr cfloat; ldv: LapackInt; tau: ptr cfloat; t: ptr cfloat; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_slarft", dynlib: libName.}
proc slarftWork*(matrixLayout: cint; direct: char; storev: char; n: LapackInt;
                k: LapackInt; v: ptr cfloat; ldv: LapackInt; tau: ptr cfloat;
                t: ptr cfloat; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_slarft_work", dynlib: libName.}
proc slarfx*(matrixLayout: cint; side: char; m: LapackInt; n: LapackInt; v: ptr cfloat;
            tau: cfloat; c: ptr cfloat; ldc: LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_slarfx", dynlib: libName.}
proc slarfxWork*(matrixLayout: cint; side: char; m: LapackInt; n: LapackInt;
                v: ptr cfloat; tau: cfloat; c: ptr cfloat; ldc: LapackInt;
                work: ptr cfloat): LapackInt {.importc: "LAPACKE_slarfx_work",
    dynlib: libName.}
proc slarnv*(idist: LapackInt; iseed: ptr LapackInt; n: LapackInt; x: ptr cfloat): LapackInt {.
    importc: "LAPACKE_slarnv", dynlib: libName.}
proc slarnvWork*(idist: LapackInt; iseed: ptr LapackInt; n: LapackInt; x: ptr cfloat): LapackInt {.
    importc: "LAPACKE_slarnv_work", dynlib: libName.}
proc slartgp*(f: cfloat; g: cfloat; cs: ptr cfloat; sn: ptr cfloat; r: ptr cfloat): LapackInt {.
    importc: "LAPACKE_slartgp", dynlib: libName.}
proc slartgpWork*(f: cfloat; g: cfloat; cs: ptr cfloat; sn: ptr cfloat; r: ptr cfloat): LapackInt {.
    importc: "LAPACKE_slartgp_work", dynlib: libName.}
proc slartgs*(x: cfloat; y: cfloat; sigma: cfloat; cs: ptr cfloat; sn: ptr cfloat): LapackInt {.
    importc: "LAPACKE_slartgs", dynlib: libName.}
proc slartgsWork*(x: cfloat; y: cfloat; sigma: cfloat; cs: ptr cfloat; sn: ptr cfloat): LapackInt {.
    importc: "LAPACKE_slartgs_work", dynlib: libName.}
proc slascl*(matrixLayout: cint; `type`: char; kl: LapackInt; ku: LapackInt;
            cfrom: cfloat; cto: cfloat; m: LapackInt; n: LapackInt; a: ptr cfloat;
            lda: LapackInt): LapackInt {.importc: "LAPACKE_slascl", dynlib: libName.}
proc slasclWork*(matrixLayout: cint; `type`: char; kl: LapackInt; ku: LapackInt;
                cfrom: cfloat; cto: cfloat; m: LapackInt; n: LapackInt; a: ptr cfloat;
                lda: LapackInt): LapackInt {.importc: "LAPACKE_slascl_work",
    dynlib: libName.}
proc slaset*(matrixLayout: cint; uplo: char; m: LapackInt; n: LapackInt; alpha: cfloat;
            beta: cfloat; a: ptr cfloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_slaset", dynlib: libName.}
proc slasetWork*(matrixLayout: cint; uplo: char; m: LapackInt; n: LapackInt;
                alpha: cfloat; beta: cfloat; a: ptr cfloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_slaset_work", dynlib: libName.}
proc slasrt*(id: char; n: LapackInt; d: ptr cfloat): LapackInt {.
    importc: "LAPACKE_slasrt", dynlib: libName.}
proc slasrtWork*(id: char; n: LapackInt; d: ptr cfloat): LapackInt {.
    importc: "LAPACKE_slasrt_work", dynlib: libName.}
proc slaswp*(matrixLayout: cint; n: LapackInt; a: ptr cfloat; lda: LapackInt;
            k1: LapackInt; k2: LapackInt; ipiv: ptr LapackInt; incx: LapackInt): LapackInt {.
    importc: "LAPACKE_slaswp", dynlib: libName.}
proc slaswpWork*(matrixLayout: cint; n: LapackInt; a: ptr cfloat; lda: LapackInt;
                k1: LapackInt; k2: LapackInt; ipiv: ptr LapackInt; incx: LapackInt): LapackInt {.
    importc: "LAPACKE_slaswp_work", dynlib: libName.}
proc slatms*(matrixLayout: cint; m: LapackInt; n: LapackInt; dist: char;
            iseed: ptr LapackInt; sym: char; d: ptr cfloat; mode: LapackInt; cond: cfloat;
            dmax: cfloat; kl: LapackInt; ku: LapackInt; pack: char; a: ptr cfloat;
            lda: LapackInt): LapackInt {.importc: "LAPACKE_slatms", dynlib: libName.}
proc slatmsWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; dist: char;
                iseed: ptr LapackInt; sym: char; d: ptr cfloat; mode: LapackInt;
                cond: cfloat; dmax: cfloat; kl: LapackInt; ku: LapackInt; pack: char;
                a: ptr cfloat; lda: LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_slatms_work", dynlib: libName.}
proc slauum*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_slauum", dynlib: libName.}
proc slauumWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                lda: LapackInt): LapackInt {.importc: "LAPACKE_slauum_work",
    dynlib: libName.}
proc sopgtr*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cfloat;
            tau: ptr cfloat; q: ptr cfloat; ldq: LapackInt): LapackInt {.
    importc: "LAPACKE_sopgtr", dynlib: libName.}
proc sopgtrWork*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cfloat;
                tau: ptr cfloat; q: ptr cfloat; ldq: LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sopgtr_work", dynlib: libName.}
proc sopmtr*(matrixLayout: cint; side: char; uplo: char; trans: char; m: LapackInt;
            n: LapackInt; ap: ptr cfloat; tau: ptr cfloat; c: ptr cfloat; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_sopmtr", dynlib: libName.}
proc sopmtrWork*(matrixLayout: cint; side: char; uplo: char; trans: char; m: LapackInt;
                n: LapackInt; ap: ptr cfloat; tau: ptr cfloat; c: ptr cfloat;
                ldc: LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sopmtr_work", dynlib: libName.}
proc sorbdb*(matrixLayout: cint; trans: char; signs: char; m: LapackInt; p: LapackInt;
            q: LapackInt; x11: ptr cfloat; ldx11: LapackInt; x12: ptr cfloat;
            ldx12: LapackInt; x21: ptr cfloat; ldx21: LapackInt; x22: ptr cfloat;
            ldx22: LapackInt; theta: ptr cfloat; phi: ptr cfloat; taup1: ptr cfloat;
            taup2: ptr cfloat; tauq1: ptr cfloat; tauq2: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sorbdb", dynlib: libName.}
proc sorbdbWork*(matrixLayout: cint; trans: char; signs: char; m: LapackInt;
                p: LapackInt; q: LapackInt; x11: ptr cfloat; ldx11: LapackInt;
                x12: ptr cfloat; ldx12: LapackInt; x21: ptr cfloat; ldx21: LapackInt;
                x22: ptr cfloat; ldx22: LapackInt; theta: ptr cfloat; phi: ptr cfloat;
                taup1: ptr cfloat; taup2: ptr cfloat; tauq1: ptr cfloat;
                tauq2: ptr cfloat; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sorbdb_work", dynlib: libName.}
proc sorcsd*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char; jobv2t: char;
            trans: char; signs: char; m: LapackInt; p: LapackInt; q: LapackInt;
            x11: ptr cfloat; ldx11: LapackInt; x12: ptr cfloat; ldx12: LapackInt;
            x21: ptr cfloat; ldx21: LapackInt; x22: ptr cfloat; ldx22: LapackInt;
            theta: ptr cfloat; u1: ptr cfloat; ldu1: LapackInt; u2: ptr cfloat;
            ldu2: LapackInt; v1t: ptr cfloat; ldv1t: LapackInt; v2t: ptr cfloat;
            ldv2t: LapackInt): LapackInt {.importc: "LAPACKE_sorcsd", dynlib: libName.}
proc sorcsd2by1*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char;
                m: LapackInt; p: LapackInt; q: LapackInt; x11: ptr cfloat;
                ldx11: LapackInt; x21: ptr cfloat; ldx21: LapackInt; theta: ptr cfloat;
                u1: ptr cfloat; ldu1: LapackInt; u2: ptr cfloat; ldu2: LapackInt;
                v1t: ptr cfloat; ldv1t: LapackInt): LapackInt {.
    importc: "LAPACKE_sorcsd2by1", dynlib: libName.}
proc sorcsd2by1Work*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char;
                    m: LapackInt; p: LapackInt; q: LapackInt; x11: ptr cfloat;
                    ldx11: LapackInt; x21: ptr cfloat; ldx21: LapackInt;
                    theta: ptr cfloat; u1: ptr cfloat; ldu1: LapackInt; u2: ptr cfloat;
                    ldu2: LapackInt; v1t: ptr cfloat; ldv1t: LapackInt;
                    work: ptr cfloat; lwork: LapackInt; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sorcsd2by1_work", dynlib: libName.}
proc sorcsdWork*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char;
                jobv2t: char; trans: char; signs: char; m: LapackInt; p: LapackInt;
                q: LapackInt; x11: ptr cfloat; ldx11: LapackInt; x12: ptr cfloat;
                ldx12: LapackInt; x21: ptr cfloat; ldx21: LapackInt; x22: ptr cfloat;
                ldx22: LapackInt; theta: ptr cfloat; u1: ptr cfloat; ldu1: LapackInt;
                u2: ptr cfloat; ldu2: LapackInt; v1t: ptr cfloat; ldv1t: LapackInt;
                v2t: ptr cfloat; ldv2t: LapackInt; work: ptr cfloat; lwork: LapackInt;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_sorcsd_work",
    dynlib: libName.}
proc sorgbr*(matrixLayout: cint; vect: char; m: LapackInt; n: LapackInt; k: LapackInt;
            a: ptr cfloat; lda: LapackInt; tau: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sorgbr", dynlib: libName.}
proc sorgbrWork*(matrixLayout: cint; vect: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr cfloat; lda: LapackInt; tau: ptr cfloat;
                work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sorgbr_work", dynlib: libName.}
proc sorghr*(matrixLayout: cint; n: LapackInt; ilo: LapackInt; ihi: LapackInt;
            a: ptr cfloat; lda: LapackInt; tau: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sorghr", dynlib: libName.}
proc sorghrWork*(matrixLayout: cint; n: LapackInt; ilo: LapackInt; ihi: LapackInt;
                a: ptr cfloat; lda: LapackInt; tau: ptr cfloat; work: ptr cfloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_sorghr_work",
    dynlib: libName.}
proc sorglq*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt; a: ptr cfloat;
            lda: LapackInt; tau: ptr cfloat): LapackInt {.importc: "LAPACKE_sorglq",
    dynlib: libName.}
proc sorglqWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
                a: ptr cfloat; lda: LapackInt; tau: ptr cfloat; work: ptr cfloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_sorglq_work",
    dynlib: libName.}
proc sorgql*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt; a: ptr cfloat;
            lda: LapackInt; tau: ptr cfloat): LapackInt {.importc: "LAPACKE_sorgql",
    dynlib: libName.}
proc sorgqlWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
                a: ptr cfloat; lda: LapackInt; tau: ptr cfloat; work: ptr cfloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_sorgql_work",
    dynlib: libName.}
proc sorgqr*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt; a: ptr cfloat;
            lda: LapackInt; tau: ptr cfloat): LapackInt {.importc: "LAPACKE_sorgqr",
    dynlib: libName.}
proc sorgqrWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
                a: ptr cfloat; lda: LapackInt; tau: ptr cfloat; work: ptr cfloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_sorgqr_work",
    dynlib: libName.}
proc sorgrq*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt; a: ptr cfloat;
            lda: LapackInt; tau: ptr cfloat): LapackInt {.importc: "LAPACKE_sorgrq",
    dynlib: libName.}
proc sorgrqWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
                a: ptr cfloat; lda: LapackInt; tau: ptr cfloat; work: ptr cfloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_sorgrq_work",
    dynlib: libName.}
proc sorgtr*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat; lda: LapackInt;
            tau: ptr cfloat): LapackInt {.importc: "LAPACKE_sorgtr", dynlib: libName.}
proc sorgtrWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; tau: ptr cfloat; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sorgtr_work", dynlib: libName.}
proc sormbr*(matrixLayout: cint; vect: char; side: char; trans: char; m: LapackInt;
            n: LapackInt; k: LapackInt; a: ptr cfloat; lda: LapackInt; tau: ptr cfloat;
            c: ptr cfloat; ldc: LapackInt): LapackInt {.importc: "LAPACKE_sormbr",
    dynlib: libName.}
proc sormbrWork*(matrixLayout: cint; vect: char; side: char; trans: char; m: LapackInt;
                n: LapackInt; k: LapackInt; a: ptr cfloat; lda: LapackInt;
                tau: ptr cfloat; c: ptr cfloat; ldc: LapackInt; work: ptr cfloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_sormbr_work",
    dynlib: libName.}
proc sormhr*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            ilo: LapackInt; ihi: LapackInt; a: ptr cfloat; lda: LapackInt;
            tau: ptr cfloat; c: ptr cfloat; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_sormhr", dynlib: libName.}
proc sormhrWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                ilo: LapackInt; ihi: LapackInt; a: ptr cfloat; lda: LapackInt;
                tau: ptr cfloat; c: ptr cfloat; ldc: LapackInt; work: ptr cfloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_sormhr_work",
    dynlib: libName.}
proc sormlq*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr cfloat; lda: LapackInt; tau: ptr cfloat; c: ptr cfloat;
            ldc: LapackInt): LapackInt {.importc: "LAPACKE_sormlq", dynlib: libName.}
proc sormlqWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr cfloat; lda: LapackInt; tau: ptr cfloat;
                c: ptr cfloat; ldc: LapackInt; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sormlq_work", dynlib: libName.}
proc sormql*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr cfloat; lda: LapackInt; tau: ptr cfloat; c: ptr cfloat;
            ldc: LapackInt): LapackInt {.importc: "LAPACKE_sormql", dynlib: libName.}
proc sormqlWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr cfloat; lda: LapackInt; tau: ptr cfloat;
                c: ptr cfloat; ldc: LapackInt; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sormql_work", dynlib: libName.}
proc sormqr*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr cfloat; lda: LapackInt; tau: ptr cfloat; c: ptr cfloat;
            ldc: LapackInt): LapackInt {.importc: "LAPACKE_sormqr", dynlib: libName.}
proc sormqrWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr cfloat; lda: LapackInt; tau: ptr cfloat;
                c: ptr cfloat; ldc: LapackInt; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sormqr_work", dynlib: libName.}
proc sormrq*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr cfloat; lda: LapackInt; tau: ptr cfloat; c: ptr cfloat;
            ldc: LapackInt): LapackInt {.importc: "LAPACKE_sormrq", dynlib: libName.}
proc sormrqWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr cfloat; lda: LapackInt; tau: ptr cfloat;
                c: ptr cfloat; ldc: LapackInt; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sormrq_work", dynlib: libName.}
proc sormrz*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; l: LapackInt; a: ptr cfloat; lda: LapackInt; tau: ptr cfloat;
            c: ptr cfloat; ldc: LapackInt): LapackInt {.importc: "LAPACKE_sormrz",
    dynlib: libName.}
proc sormrzWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; l: LapackInt; a: ptr cfloat; lda: LapackInt;
                tau: ptr cfloat; c: ptr cfloat; ldc: LapackInt; work: ptr cfloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_sormrz_work",
    dynlib: libName.}
proc sormtr*(matrixLayout: cint; side: char; uplo: char; trans: char; m: LapackInt;
            n: LapackInt; a: ptr cfloat; lda: LapackInt; tau: ptr cfloat; c: ptr cfloat;
            ldc: LapackInt): LapackInt {.importc: "LAPACKE_sormtr", dynlib: libName.}
proc sormtrWork*(matrixLayout: cint; side: char; uplo: char; trans: char; m: LapackInt;
                n: LapackInt; a: ptr cfloat; lda: LapackInt; tau: ptr cfloat;
                c: ptr cfloat; ldc: LapackInt; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sormtr_work", dynlib: libName.}
proc spbcon*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt; ab: ptr cfloat;
            ldab: LapackInt; anorm: cfloat; rcond: ptr cfloat): LapackInt {.
    importc: "LAPACKE_spbcon", dynlib: libName.}
proc spbconWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
                ab: ptr cfloat; ldab: LapackInt; anorm: cfloat; rcond: ptr cfloat;
                work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_spbcon_work", dynlib: libName.}
proc spbequ*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt; ab: ptr cfloat;
            ldab: LapackInt; s: ptr cfloat; scond: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_spbequ", dynlib: libName.}
proc spbequWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
                ab: ptr cfloat; ldab: LapackInt; s: ptr cfloat; scond: ptr cfloat;
                amax: ptr cfloat): LapackInt {.importc: "LAPACKE_spbequ_work",
    dynlib: libName.}
proc spbrfs*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
            nrhs: LapackInt; ab: ptr cfloat; ldab: LapackInt; afb: ptr cfloat;
            ldafb: LapackInt; b: ptr cfloat; ldb: LapackInt; x: ptr cfloat;
            ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.
    importc: "LAPACKE_spbrfs", dynlib: libName.}
proc spbrfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
                nrhs: LapackInt; ab: ptr cfloat; ldab: LapackInt; afb: ptr cfloat;
                ldafb: LapackInt; b: ptr cfloat; ldb: LapackInt; x: ptr cfloat;
                ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat; work: ptr cfloat;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_spbrfs_work",
    dynlib: libName.}
proc spbstf*(matrixLayout: cint; uplo: char; n: LapackInt; kb: LapackInt; bb: ptr cfloat;
            ldbb: LapackInt): LapackInt {.importc: "LAPACKE_spbstf", dynlib: libName.}
proc spbstfWork*(matrixLayout: cint; uplo: char; n: LapackInt; kb: LapackInt;
                bb: ptr cfloat; ldbb: LapackInt): LapackInt {.
    importc: "LAPACKE_spbstf_work", dynlib: libName.}
proc spbsv*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt; nrhs: LapackInt;
           ab: ptr cfloat; ldab: LapackInt; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_spbsv", dynlib: libName.}
proc spbsvWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
               nrhs: LapackInt; ab: ptr cfloat; ldab: LapackInt; b: ptr cfloat;
               ldb: LapackInt): LapackInt {.importc: "LAPACKE_spbsv_work",
    dynlib: libName.}
proc spbsvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; kd: LapackInt;
            nrhs: LapackInt; ab: ptr cfloat; ldab: LapackInt; afb: ptr cfloat;
            ldafb: LapackInt; equed: cstring; s: ptr cfloat; b: ptr cfloat;
            ldb: LapackInt; x: ptr cfloat; ldx: LapackInt; rcond: ptr cfloat;
            ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.importc: "LAPACKE_spbsvx",
    dynlib: libName.}
proc spbsvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; kd: LapackInt;
                nrhs: LapackInt; ab: ptr cfloat; ldab: LapackInt; afb: ptr cfloat;
                ldafb: LapackInt; equed: cstring; s: ptr cfloat; b: ptr cfloat;
                ldb: LapackInt; x: ptr cfloat; ldx: LapackInt; rcond: ptr cfloat;
                ferr: ptr cfloat; berr: ptr cfloat; work: ptr cfloat;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_spbsvx_work",
    dynlib: libName.}
proc spbtrf*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt; ab: ptr cfloat;
            ldab: LapackInt): LapackInt {.importc: "LAPACKE_spbtrf", dynlib: libName.}
proc spbtrfWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
                ab: ptr cfloat; ldab: LapackInt): LapackInt {.
    importc: "LAPACKE_spbtrf_work", dynlib: libName.}
proc spbtrs*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
            nrhs: LapackInt; ab: ptr cfloat; ldab: LapackInt; b: ptr cfloat;
            ldb: LapackInt): LapackInt {.importc: "LAPACKE_spbtrs", dynlib: libName.}
proc spbtrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
                nrhs: LapackInt; ab: ptr cfloat; ldab: LapackInt; b: ptr cfloat;
                ldb: LapackInt): LapackInt {.importc: "LAPACKE_spbtrs_work",
    dynlib: libName.}
proc spftrf*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt; a: ptr cfloat): LapackInt {.
    importc: "LAPACKE_spftrf", dynlib: libName.}
proc spftrfWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                a: ptr cfloat): LapackInt {.importc: "LAPACKE_spftrf_work",
                                        dynlib: libName.}
proc spftri*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt; a: ptr cfloat): LapackInt {.
    importc: "LAPACKE_spftri", dynlib: libName.}
proc spftriWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                a: ptr cfloat): LapackInt {.importc: "LAPACKE_spftri_work",
                                        dynlib: libName.}
proc spftrs*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr cfloat; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_spftrs", dynlib: libName.}
proc spftrsWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; a: ptr cfloat; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_spftrs_work", dynlib: libName.}
proc spocon*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat; lda: LapackInt;
            anorm: cfloat; rcond: ptr cfloat): LapackInt {.importc: "LAPACKE_spocon",
    dynlib: libName.}
proc spoconWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; anorm: cfloat; rcond: ptr cfloat; work: ptr cfloat;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_spocon_work",
    dynlib: libName.}
proc spoequ*(matrixLayout: cint; n: LapackInt; a: ptr cfloat; lda: LapackInt;
            s: ptr cfloat; scond: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_spoequ", dynlib: libName.}
proc spoequWork*(matrixLayout: cint; n: LapackInt; a: ptr cfloat; lda: LapackInt;
                s: ptr cfloat; scond: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_spoequ_work", dynlib: libName.}
proc spoequb*(matrixLayout: cint; n: LapackInt; a: ptr cfloat; lda: LapackInt;
             s: ptr cfloat; scond: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_spoequb", dynlib: libName.}
proc spoequbWork*(matrixLayout: cint; n: LapackInt; a: ptr cfloat; lda: LapackInt;
                 s: ptr cfloat; scond: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_spoequb_work", dynlib: libName.}
proc sporfs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr cfloat; lda: LapackInt; af: ptr cfloat; ldaf: LapackInt; b: ptr cfloat;
            ldb: LapackInt; x: ptr cfloat; ldx: LapackInt; ferr: ptr cfloat;
            berr: ptr cfloat): LapackInt {.importc: "LAPACKE_sporfs", dynlib: libName.}
proc sporfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr cfloat; lda: LapackInt; af: ptr cfloat; ldaf: LapackInt;
                b: ptr cfloat; ldb: LapackInt; x: ptr cfloat; ldx: LapackInt;
                ferr: ptr cfloat; berr: ptr cfloat; work: ptr cfloat;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_sporfs_work",
    dynlib: libName.}
proc sporfsx*(matrixLayout: cint; uplo: char; equed: char; n: LapackInt; nrhs: LapackInt;
             a: ptr cfloat; lda: LapackInt; af: ptr cfloat; ldaf: LapackInt;
             s: ptr cfloat; b: ptr cfloat; ldb: LapackInt; x: ptr cfloat; ldx: LapackInt;
             rcond: ptr cfloat; berr: ptr cfloat; nErrBnds: LapackInt;
             errBndsNorm: ptr cfloat; errBndsComp: ptr cfloat; nparams: LapackInt;
             params: ptr cfloat): LapackInt {.importc: "LAPACKE_sporfsx",
    dynlib: libName.}
proc sporfsxWork*(matrixLayout: cint; uplo: char; equed: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr cfloat; lda: LapackInt; af: ptr cfloat;
                 ldaf: LapackInt; s: ptr cfloat; b: ptr cfloat; ldb: LapackInt;
                 x: ptr cfloat; ldx: LapackInt; rcond: ptr cfloat; berr: ptr cfloat;
                 nErrBnds: LapackInt; errBndsNorm: ptr cfloat;
                 errBndsComp: ptr cfloat; nparams: LapackInt; params: ptr cfloat;
                 work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sporfsx_work", dynlib: libName.}
proc sposv*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt; a: ptr cfloat;
           lda: LapackInt; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_sposv", dynlib: libName.}
proc sposvWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_sposv_work", dynlib: libName.}
proc sposvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr cfloat; lda: LapackInt; af: ptr cfloat; ldaf: LapackInt;
            equed: cstring; s: ptr cfloat; b: ptr cfloat; ldb: LapackInt; x: ptr cfloat;
            ldx: LapackInt; rcond: ptr cfloat; ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sposvx", dynlib: libName.}
proc sposvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; a: ptr cfloat; lda: LapackInt; af: ptr cfloat;
                ldaf: LapackInt; equed: cstring; s: ptr cfloat; b: ptr cfloat;
                ldb: LapackInt; x: ptr cfloat; ldx: LapackInt; rcond: ptr cfloat;
                ferr: ptr cfloat; berr: ptr cfloat; work: ptr cfloat;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_sposvx_work",
    dynlib: libName.}
proc sposvxx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr cfloat; lda: LapackInt; af: ptr cfloat; ldaf: LapackInt;
             equed: cstring; s: ptr cfloat; b: ptr cfloat; ldb: LapackInt; x: ptr cfloat;
             ldx: LapackInt; rcond: ptr cfloat; rpvgrw: ptr cfloat; berr: ptr cfloat;
             nErrBnds: LapackInt; errBndsNorm: ptr cfloat; errBndsComp: ptr cfloat;
             nparams: LapackInt; params: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sposvxx", dynlib: libName.}
proc sposvxxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr cfloat; lda: LapackInt; af: ptr cfloat;
                 ldaf: LapackInt; equed: cstring; s: ptr cfloat; b: ptr cfloat;
                 ldb: LapackInt; x: ptr cfloat; ldx: LapackInt; rcond: ptr cfloat;
                 rpvgrw: ptr cfloat; berr: ptr cfloat; nErrBnds: LapackInt;
                 errBndsNorm: ptr cfloat; errBndsComp: ptr cfloat; nparams: LapackInt;
                 params: ptr cfloat; work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sposvxx_work", dynlib: libName.}
proc spotrf*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_spotrf", dynlib: libName.}
proc spotrf2*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_spotrf2", dynlib: libName.}
proc spotrf2Work*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                 lda: LapackInt): LapackInt {.importc: "LAPACKE_spotrf2_work",
    dynlib: libName.}
proc spotrfWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                lda: LapackInt): LapackInt {.importc: "LAPACKE_spotrf_work",
    dynlib: libName.}
proc spotri*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_spotri", dynlib: libName.}
proc spotriWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                lda: LapackInt): LapackInt {.importc: "LAPACKE_spotri_work",
    dynlib: libName.}
proc spotrs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_spotrs", dynlib: libName.}
proc spotrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_spotrs_work", dynlib: libName.}
proc sppcon*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cfloat; anorm: cfloat;
            rcond: ptr cfloat): LapackInt {.importc: "LAPACKE_sppcon", dynlib: libName.}
proc sppconWork*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cfloat;
                anorm: cfloat; rcond: ptr cfloat; work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sppcon_work", dynlib: libName.}
proc sppequ*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cfloat; s: ptr cfloat;
            scond: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sppequ", dynlib: libName.}
proc sppequWork*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cfloat;
                s: ptr cfloat; scond: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sppequ_work", dynlib: libName.}
proc spprfs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr cfloat; afp: ptr cfloat; b: ptr cfloat; ldb: LapackInt; x: ptr cfloat;
            ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.
    importc: "LAPACKE_spprfs", dynlib: libName.}
proc spprfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                ap: ptr cfloat; afp: ptr cfloat; b: ptr cfloat; ldb: LapackInt;
                x: ptr cfloat; ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat;
                work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_spprfs_work", dynlib: libName.}
proc sppsv*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
           ap: ptr cfloat; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_sppsv", dynlib: libName.}
proc sppsvWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               ap: ptr cfloat; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_sppsv_work", dynlib: libName.}
proc sppsvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr cfloat; afp: ptr cfloat; equed: cstring; s: ptr cfloat; b: ptr cfloat;
            ldb: LapackInt; x: ptr cfloat; ldx: LapackInt; rcond: ptr cfloat;
            ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.importc: "LAPACKE_sppsvx",
    dynlib: libName.}
proc sppsvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; ap: ptr cfloat; afp: ptr cfloat; equed: cstring;
                s: ptr cfloat; b: ptr cfloat; ldb: LapackInt; x: ptr cfloat;
                ldx: LapackInt; rcond: ptr cfloat; ferr: ptr cfloat; berr: ptr cfloat;
                work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sppsvx_work", dynlib: libName.}
proc spptrf*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cfloat): LapackInt {.
    importc: "LAPACKE_spptrf", dynlib: libName.}
proc spptrfWork*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cfloat): LapackInt {.
    importc: "LAPACKE_spptrf_work", dynlib: libName.}
proc spptri*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cfloat): LapackInt {.
    importc: "LAPACKE_spptri", dynlib: libName.}
proc spptriWork*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cfloat): LapackInt {.
    importc: "LAPACKE_spptri_work", dynlib: libName.}
proc spptrs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr cfloat; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_spptrs", dynlib: libName.}
proc spptrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                ap: ptr cfloat; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_spptrs_work", dynlib: libName.}
proc spstrf*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat; lda: LapackInt;
            piv: ptr LapackInt; rank: ptr LapackInt; tol: cfloat): LapackInt {.
    importc: "LAPACKE_spstrf", dynlib: libName.}
proc spstrfWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; piv: ptr LapackInt; rank: ptr LapackInt; tol: cfloat;
                work: ptr cfloat): LapackInt {.importc: "LAPACKE_spstrf_work",
    dynlib: libName.}
proc sptcon*(n: LapackInt; d: ptr cfloat; e: ptr cfloat; anorm: cfloat; rcond: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sptcon", dynlib: libName.}
proc sptconWork*(n: LapackInt; d: ptr cfloat; e: ptr cfloat; anorm: cfloat;
                rcond: ptr cfloat; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sptcon_work", dynlib: libName.}
proc spteqr*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cfloat; e: ptr cfloat;
            z: ptr cfloat; ldz: LapackInt): LapackInt {.importc: "LAPACKE_spteqr",
    dynlib: libName.}
proc spteqrWork*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cfloat;
                e: ptr cfloat; z: ptr cfloat; ldz: LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_spteqr_work", dynlib: libName.}
proc sptrfs*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; d: ptr cfloat;
            e: ptr cfloat; df: ptr cfloat; ef: ptr cfloat; b: ptr cfloat; ldb: LapackInt;
            x: ptr cfloat; ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sptrfs", dynlib: libName.}
proc sptrfsWork*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; d: ptr cfloat;
                e: ptr cfloat; df: ptr cfloat; ef: ptr cfloat; b: ptr cfloat;
                ldb: LapackInt; x: ptr cfloat; ldx: LapackInt; ferr: ptr cfloat;
                berr: ptr cfloat; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sptrfs_work", dynlib: libName.}
proc sptsv*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; d: ptr cfloat;
           e: ptr cfloat; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_sptsv", dynlib: libName.}
proc sptsvWork*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; d: ptr cfloat;
               e: ptr cfloat; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_sptsv_work", dynlib: libName.}
proc sptsvx*(matrixLayout: cint; fact: char; n: LapackInt; nrhs: LapackInt;
            d: ptr cfloat; e: ptr cfloat; df: ptr cfloat; ef: ptr cfloat; b: ptr cfloat;
            ldb: LapackInt; x: ptr cfloat; ldx: LapackInt; rcond: ptr cfloat;
            ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.importc: "LAPACKE_sptsvx",
    dynlib: libName.}
proc sptsvxWork*(matrixLayout: cint; fact: char; n: LapackInt; nrhs: LapackInt;
                d: ptr cfloat; e: ptr cfloat; df: ptr cfloat; ef: ptr cfloat; b: ptr cfloat;
                ldb: LapackInt; x: ptr cfloat; ldx: LapackInt; rcond: ptr cfloat;
                ferr: ptr cfloat; berr: ptr cfloat; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sptsvx_work", dynlib: libName.}
proc spttrf*(n: LapackInt; d: ptr cfloat; e: ptr cfloat): LapackInt {.
    importc: "LAPACKE_spttrf", dynlib: libName.}
proc spttrfWork*(n: LapackInt; d: ptr cfloat; e: ptr cfloat): LapackInt {.
    importc: "LAPACKE_spttrf_work", dynlib: libName.}
proc spttrs*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; d: ptr cfloat;
            e: ptr cfloat; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_spttrs", dynlib: libName.}
proc spttrsWork*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; d: ptr cfloat;
                e: ptr cfloat; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_spttrs_work", dynlib: libName.}
proc ssbev*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; kd: LapackInt;
           ab: ptr cfloat; ldab: LapackInt; w: ptr cfloat; z: ptr cfloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_ssbev", dynlib: libName.}
proc ssbevWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; kd: LapackInt;
               ab: ptr cfloat; ldab: LapackInt; w: ptr cfloat; z: ptr cfloat;
               ldz: LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssbev_work", dynlib: libName.}
proc ssbevd*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; kd: LapackInt;
            ab: ptr cfloat; ldab: LapackInt; w: ptr cfloat; z: ptr cfloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_ssbevd", dynlib: libName.}
proc ssbevdWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; kd: LapackInt;
                ab: ptr cfloat; ldab: LapackInt; w: ptr cfloat; z: ptr cfloat;
                ldz: LapackInt; work: ptr cfloat; lwork: LapackInt;
                iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ssbevd_work", dynlib: libName.}
proc ssbevx*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
            kd: LapackInt; ab: ptr cfloat; ldab: LapackInt; q: ptr cfloat; ldq: LapackInt;
            vl: cfloat; vu: cfloat; il: LapackInt; iu: LapackInt; abstol: cfloat;
            m: ptr LapackInt; w: ptr cfloat; z: ptr cfloat; ldz: LapackInt;
            ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_ssbevx",
    dynlib: libName.}
proc ssbevxWork*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                kd: LapackInt; ab: ptr cfloat; ldab: LapackInt; q: ptr cfloat;
                ldq: LapackInt; vl: cfloat; vu: cfloat; il: LapackInt; iu: LapackInt;
                abstol: cfloat; m: ptr LapackInt; w: ptr cfloat; z: ptr cfloat;
                ldz: LapackInt; work: ptr cfloat; iwork: ptr LapackInt;
                ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_ssbevx_work",
    dynlib: libName.}
proc ssbgst*(matrixLayout: cint; vect: char; uplo: char; n: LapackInt; ka: LapackInt;
            kb: LapackInt; ab: ptr cfloat; ldab: LapackInt; bb: ptr cfloat;
            ldbb: LapackInt; x: ptr cfloat; ldx: LapackInt): LapackInt {.
    importc: "LAPACKE_ssbgst", dynlib: libName.}
proc ssbgstWork*(matrixLayout: cint; vect: char; uplo: char; n: LapackInt; ka: LapackInt;
                kb: LapackInt; ab: ptr cfloat; ldab: LapackInt; bb: ptr cfloat;
                ldbb: LapackInt; x: ptr cfloat; ldx: LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssbgst_work", dynlib: libName.}
proc ssbgv*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; ka: LapackInt;
           kb: LapackInt; ab: ptr cfloat; ldab: LapackInt; bb: ptr cfloat;
           ldbb: LapackInt; w: ptr cfloat; z: ptr cfloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_ssbgv", dynlib: libName.}
proc ssbgvWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; ka: LapackInt;
               kb: LapackInt; ab: ptr cfloat; ldab: LapackInt; bb: ptr cfloat;
               ldbb: LapackInt; w: ptr cfloat; z: ptr cfloat; ldz: LapackInt;
               work: ptr cfloat): LapackInt {.importc: "LAPACKE_ssbgv_work",
    dynlib: libName.}
proc ssbgvd*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; ka: LapackInt;
            kb: LapackInt; ab: ptr cfloat; ldab: LapackInt; bb: ptr cfloat;
            ldbb: LapackInt; w: ptr cfloat; z: ptr cfloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_ssbgvd", dynlib: libName.}
proc ssbgvdWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; ka: LapackInt;
                kb: LapackInt; ab: ptr cfloat; ldab: LapackInt; bb: ptr cfloat;
                ldbb: LapackInt; w: ptr cfloat; z: ptr cfloat; ldz: LapackInt;
                work: ptr cfloat; lwork: LapackInt; iwork: ptr LapackInt;
                liwork: LapackInt): LapackInt {.importc: "LAPACKE_ssbgvd_work",
    dynlib: libName.}
proc ssbgvx*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
            ka: LapackInt; kb: LapackInt; ab: ptr cfloat; ldab: LapackInt;
            bb: ptr cfloat; ldbb: LapackInt; q: ptr cfloat; ldq: LapackInt; vl: cfloat;
            vu: cfloat; il: LapackInt; iu: LapackInt; abstol: cfloat; m: ptr LapackInt;
            w: ptr cfloat; z: ptr cfloat; ldz: LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ssbgvx", dynlib: libName.}
proc ssbgvxWork*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                ka: LapackInt; kb: LapackInt; ab: ptr cfloat; ldab: LapackInt;
                bb: ptr cfloat; ldbb: LapackInt; q: ptr cfloat; ldq: LapackInt;
                vl: cfloat; vu: cfloat; il: LapackInt; iu: LapackInt; abstol: cfloat;
                m: ptr LapackInt; w: ptr cfloat; z: ptr cfloat; ldz: LapackInt;
                work: ptr cfloat; iwork: ptr LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ssbgvx_work", dynlib: libName.}
proc ssbtrd*(matrixLayout: cint; vect: char; uplo: char; n: LapackInt; kd: LapackInt;
            ab: ptr cfloat; ldab: LapackInt; d: ptr cfloat; e: ptr cfloat; q: ptr cfloat;
            ldq: LapackInt): LapackInt {.importc: "LAPACKE_ssbtrd", dynlib: libName.}
proc ssbtrdWork*(matrixLayout: cint; vect: char; uplo: char; n: LapackInt; kd: LapackInt;
                ab: ptr cfloat; ldab: LapackInt; d: ptr cfloat; e: ptr cfloat;
                q: ptr cfloat; ldq: LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssbtrd_work", dynlib: libName.}
proc ssfrk*(matrixLayout: cint; transr: char; uplo: char; trans: char; n: LapackInt;
           k: LapackInt; alpha: cfloat; a: ptr cfloat; lda: LapackInt; beta: cfloat;
           c: ptr cfloat): LapackInt {.importc: "LAPACKE_ssfrk", dynlib: libName.}
proc ssfrkWork*(matrixLayout: cint; transr: char; uplo: char; trans: char; n: LapackInt;
               k: LapackInt; alpha: cfloat; a: ptr cfloat; lda: LapackInt; beta: cfloat;
               c: ptr cfloat): LapackInt {.importc: "LAPACKE_ssfrk_work",
                                       dynlib: libName.}
proc sspcon*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cfloat;
            ipiv: ptr LapackInt; anorm: cfloat; rcond: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sspcon", dynlib: libName.}
proc sspconWork*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cfloat;
                ipiv: ptr LapackInt; anorm: cfloat; rcond: ptr cfloat; work: ptr cfloat;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_sspcon_work",
    dynlib: libName.}
proc sspev*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; ap: ptr cfloat;
           w: ptr cfloat; z: ptr cfloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_sspev", dynlib: libName.}
proc sspevWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; ap: ptr cfloat;
               w: ptr cfloat; z: ptr cfloat; ldz: LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sspev_work", dynlib: libName.}
proc sspevd*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; ap: ptr cfloat;
            w: ptr cfloat; z: ptr cfloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_sspevd", dynlib: libName.}
proc sspevdWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                ap: ptr cfloat; w: ptr cfloat; z: ptr cfloat; ldz: LapackInt;
                work: ptr cfloat; lwork: LapackInt; iwork: ptr LapackInt;
                liwork: LapackInt): LapackInt {.importc: "LAPACKE_sspevd_work",
    dynlib: libName.}
proc sspevx*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
            ap: ptr cfloat; vl: cfloat; vu: cfloat; il: LapackInt; iu: LapackInt;
            abstol: cfloat; m: ptr LapackInt; w: ptr cfloat; z: ptr cfloat; ldz: LapackInt;
            ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_sspevx",
    dynlib: libName.}
proc sspevxWork*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                ap: ptr cfloat; vl: cfloat; vu: cfloat; il: LapackInt; iu: LapackInt;
                abstol: cfloat; m: ptr LapackInt; w: ptr cfloat; z: ptr cfloat;
                ldz: LapackInt; work: ptr cfloat; iwork: ptr LapackInt;
                ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_sspevx_work",
    dynlib: libName.}
proc sspgst*(matrixLayout: cint; itype: LapackInt; uplo: char; n: LapackInt;
            ap: ptr cfloat; bp: ptr cfloat): LapackInt {.importc: "LAPACKE_sspgst",
    dynlib: libName.}
proc sspgstWork*(matrixLayout: cint; itype: LapackInt; uplo: char; n: LapackInt;
                ap: ptr cfloat; bp: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sspgst_work", dynlib: libName.}
proc sspgv*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char; n: LapackInt;
           ap: ptr cfloat; bp: ptr cfloat; w: ptr cfloat; z: ptr cfloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_sspgv", dynlib: libName.}
proc sspgvWork*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
               n: LapackInt; ap: ptr cfloat; bp: ptr cfloat; w: ptr cfloat; z: ptr cfloat;
               ldz: LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sspgv_work", dynlib: libName.}
proc sspgvd*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char; n: LapackInt;
            ap: ptr cfloat; bp: ptr cfloat; w: ptr cfloat; z: ptr cfloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_sspgvd", dynlib: libName.}
proc sspgvdWork*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
                n: LapackInt; ap: ptr cfloat; bp: ptr cfloat; w: ptr cfloat; z: ptr cfloat;
                ldz: LapackInt; work: ptr cfloat; lwork: LapackInt;
                iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sspgvd_work", dynlib: libName.}
proc sspgvx*(matrixLayout: cint; itype: LapackInt; jobz: char; range: char; uplo: char;
            n: LapackInt; ap: ptr cfloat; bp: ptr cfloat; vl: cfloat; vu: cfloat;
            il: LapackInt; iu: LapackInt; abstol: cfloat; m: ptr LapackInt;
            w: ptr cfloat; z: ptr cfloat; ldz: LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sspgvx", dynlib: libName.}
proc sspgvxWork*(matrixLayout: cint; itype: LapackInt; jobz: char; range: char;
                uplo: char; n: LapackInt; ap: ptr cfloat; bp: ptr cfloat; vl: cfloat;
                vu: cfloat; il: LapackInt; iu: LapackInt; abstol: cfloat;
                m: ptr LapackInt; w: ptr cfloat; z: ptr cfloat; ldz: LapackInt;
                work: ptr cfloat; iwork: ptr LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sspgvx_work", dynlib: libName.}
proc ssprfs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr cfloat; afp: ptr cfloat; ipiv: ptr LapackInt; b: ptr cfloat;
            ldb: LapackInt; x: ptr cfloat; ldx: LapackInt; ferr: ptr cfloat;
            berr: ptr cfloat): LapackInt {.importc: "LAPACKE_ssprfs", dynlib: libName.}
proc ssprfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                ap: ptr cfloat; afp: ptr cfloat; ipiv: ptr LapackInt; b: ptr cfloat;
                ldb: LapackInt; x: ptr cfloat; ldx: LapackInt; ferr: ptr cfloat;
                berr: ptr cfloat; work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ssprfs_work", dynlib: libName.}
proc sspsv*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
           ap: ptr cfloat; ipiv: ptr LapackInt; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_sspsv", dynlib: libName.}
proc sspsvWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               ap: ptr cfloat; ipiv: ptr LapackInt; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_sspsv_work", dynlib: libName.}
proc sspsvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr cfloat; afp: ptr cfloat; ipiv: ptr LapackInt; b: ptr cfloat;
            ldb: LapackInt; x: ptr cfloat; ldx: LapackInt; rcond: ptr cfloat;
            ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.importc: "LAPACKE_sspsvx",
    dynlib: libName.}
proc sspsvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; ap: ptr cfloat; afp: ptr cfloat; ipiv: ptr LapackInt;
                b: ptr cfloat; ldb: LapackInt; x: ptr cfloat; ldx: LapackInt;
                rcond: ptr cfloat; ferr: ptr cfloat; berr: ptr cfloat; work: ptr cfloat;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_sspsvx_work",
    dynlib: libName.}
proc ssptrd*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cfloat; d: ptr cfloat;
            e: ptr cfloat; tau: ptr cfloat): LapackInt {.importc: "LAPACKE_ssptrd",
    dynlib: libName.}
proc ssptrdWork*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cfloat;
                d: ptr cfloat; e: ptr cfloat; tau: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssptrd_work", dynlib: libName.}
proc ssptrf*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cfloat;
            ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_ssptrf",
    dynlib: libName.}
proc ssptrfWork*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cfloat;
                ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_ssptrf_work",
    dynlib: libName.}
proc ssptri*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cfloat;
            ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_ssptri",
    dynlib: libName.}
proc ssptriWork*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cfloat;
                ipiv: ptr LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssptri_work", dynlib: libName.}
proc ssptrs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr cfloat; ipiv: ptr LapackInt; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ssptrs", dynlib: libName.}
proc ssptrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                ap: ptr cfloat; ipiv: ptr LapackInt; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ssptrs_work", dynlib: libName.}
proc sstebz*(range: char; order: char; n: LapackInt; vl: cfloat; vu: cfloat; il: LapackInt;
            iu: LapackInt; abstol: cfloat; d: ptr cfloat; e: ptr cfloat; m: ptr LapackInt;
            nsplit: ptr LapackInt; w: ptr cfloat; iblock: ptr LapackInt;
            isplit: ptr LapackInt): LapackInt {.importc: "LAPACKE_sstebz",
    dynlib: libName.}
proc sstebzWork*(range: char; order: char; n: LapackInt; vl: cfloat; vu: cfloat;
                il: LapackInt; iu: LapackInt; abstol: cfloat; d: ptr cfloat;
                e: ptr cfloat; m: ptr LapackInt; nsplit: ptr LapackInt; w: ptr cfloat;
                iblock: ptr LapackInt; isplit: ptr LapackInt; work: ptr cfloat;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_sstebz_work",
    dynlib: libName.}
proc sstedc*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cfloat; e: ptr cfloat;
            z: ptr cfloat; ldz: LapackInt): LapackInt {.importc: "LAPACKE_sstedc",
    dynlib: libName.}
proc sstedcWork*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cfloat;
                e: ptr cfloat; z: ptr cfloat; ldz: LapackInt; work: ptr cfloat;
                lwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sstedc_work", dynlib: libName.}
proc sstegr*(matrixLayout: cint; jobz: char; range: char; n: LapackInt; d: ptr cfloat;
            e: ptr cfloat; vl: cfloat; vu: cfloat; il: LapackInt; iu: LapackInt;
            abstol: cfloat; m: ptr LapackInt; w: ptr cfloat; z: ptr cfloat; ldz: LapackInt;
            isuppz: ptr LapackInt): LapackInt {.importc: "LAPACKE_sstegr",
    dynlib: libName.}
proc sstegrWork*(matrixLayout: cint; jobz: char; range: char; n: LapackInt;
                d: ptr cfloat; e: ptr cfloat; vl: cfloat; vu: cfloat; il: LapackInt;
                iu: LapackInt; abstol: cfloat; m: ptr LapackInt; w: ptr cfloat;
                z: ptr cfloat; ldz: LapackInt; isuppz: ptr LapackInt; work: ptr cfloat;
                lwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sstegr_work", dynlib: libName.}
proc sstein*(matrixLayout: cint; n: LapackInt; d: ptr cfloat; e: ptr cfloat; m: LapackInt;
            w: ptr cfloat; iblock: ptr LapackInt; isplit: ptr LapackInt; z: ptr cfloat;
            ldz: LapackInt; ifailv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sstein", dynlib: libName.}
proc ssteinWork*(matrixLayout: cint; n: LapackInt; d: ptr cfloat; e: ptr cfloat;
                m: LapackInt; w: ptr cfloat; iblock: ptr LapackInt;
                isplit: ptr LapackInt; z: ptr cfloat; ldz: LapackInt; work: ptr cfloat;
                iwork: ptr LapackInt; ifailv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_sstein_work", dynlib: libName.}
proc sstemr*(matrixLayout: cint; jobz: char; range: char; n: LapackInt; d: ptr cfloat;
            e: ptr cfloat; vl: cfloat; vu: cfloat; il: LapackInt; iu: LapackInt;
            m: ptr LapackInt; w: ptr cfloat; z: ptr cfloat; ldz: LapackInt; nzc: LapackInt;
            isuppz: ptr LapackInt; tryrac: ptr LapackLogical): LapackInt {.
    importc: "LAPACKE_sstemr", dynlib: libName.}
proc sstemrWork*(matrixLayout: cint; jobz: char; range: char; n: LapackInt;
                d: ptr cfloat; e: ptr cfloat; vl: cfloat; vu: cfloat; il: LapackInt;
                iu: LapackInt; m: ptr LapackInt; w: ptr cfloat; z: ptr cfloat;
                ldz: LapackInt; nzc: LapackInt; isuppz: ptr LapackInt;
                tryrac: ptr LapackLogical; work: ptr cfloat; lwork: LapackInt;
                iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sstemr_work", dynlib: libName.}
proc ssteqr*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cfloat; e: ptr cfloat;
            z: ptr cfloat; ldz: LapackInt): LapackInt {.importc: "LAPACKE_ssteqr",
    dynlib: libName.}
proc ssteqrWork*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cfloat;
                e: ptr cfloat; z: ptr cfloat; ldz: LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssteqr_work", dynlib: libName.}
proc ssterf*(n: LapackInt; d: ptr cfloat; e: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssterf", dynlib: libName.}
proc ssterfWork*(n: LapackInt; d: ptr cfloat; e: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssterf_work", dynlib: libName.}
proc sstev*(matrixLayout: cint; jobz: char; n: LapackInt; d: ptr cfloat; e: ptr cfloat;
           z: ptr cfloat; ldz: LapackInt): LapackInt {.importc: "LAPACKE_sstev",
    dynlib: libName.}
proc sstevWork*(matrixLayout: cint; jobz: char; n: LapackInt; d: ptr cfloat;
               e: ptr cfloat; z: ptr cfloat; ldz: LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_sstev_work", dynlib: libName.}
proc sstevd*(matrixLayout: cint; jobz: char; n: LapackInt; d: ptr cfloat; e: ptr cfloat;
            z: ptr cfloat; ldz: LapackInt): LapackInt {.importc: "LAPACKE_sstevd",
    dynlib: libName.}
proc sstevdWork*(matrixLayout: cint; jobz: char; n: LapackInt; d: ptr cfloat;
                e: ptr cfloat; z: ptr cfloat; ldz: LapackInt; work: ptr cfloat;
                lwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sstevd_work", dynlib: libName.}
proc sstevr*(matrixLayout: cint; jobz: char; range: char; n: LapackInt; d: ptr cfloat;
            e: ptr cfloat; vl: cfloat; vu: cfloat; il: LapackInt; iu: LapackInt;
            abstol: cfloat; m: ptr LapackInt; w: ptr cfloat; z: ptr cfloat; ldz: LapackInt;
            isuppz: ptr LapackInt): LapackInt {.importc: "LAPACKE_sstevr",
    dynlib: libName.}
proc sstevrWork*(matrixLayout: cint; jobz: char; range: char; n: LapackInt;
                d: ptr cfloat; e: ptr cfloat; vl: cfloat; vu: cfloat; il: LapackInt;
                iu: LapackInt; abstol: cfloat; m: ptr LapackInt; w: ptr cfloat;
                z: ptr cfloat; ldz: LapackInt; isuppz: ptr LapackInt; work: ptr cfloat;
                lwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sstevr_work", dynlib: libName.}
proc sstevx*(matrixLayout: cint; jobz: char; range: char; n: LapackInt; d: ptr cfloat;
            e: ptr cfloat; vl: cfloat; vu: cfloat; il: LapackInt; iu: LapackInt;
            abstol: cfloat; m: ptr LapackInt; w: ptr cfloat; z: ptr cfloat; ldz: LapackInt;
            ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_sstevx",
    dynlib: libName.}
proc sstevxWork*(matrixLayout: cint; jobz: char; range: char; n: LapackInt;
                d: ptr cfloat; e: ptr cfloat; vl: cfloat; vu: cfloat; il: LapackInt;
                iu: LapackInt; abstol: cfloat; m: ptr LapackInt; w: ptr cfloat;
                z: ptr cfloat; ldz: LapackInt; work: ptr cfloat; iwork: ptr LapackInt;
                ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_sstevx_work",
    dynlib: libName.}
proc ssycon*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat; lda: LapackInt;
            ipiv: ptr LapackInt; anorm: cfloat; rcond: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssycon", dynlib: libName.}
proc ssyconWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; ipiv: ptr LapackInt; anorm: cfloat; rcond: ptr cfloat;
                work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ssycon_work", dynlib: libName.}
proc ssyconv*(matrixLayout: cint; uplo: char; way: char; n: LapackInt; a: ptr cfloat;
             lda: LapackInt; ipiv: ptr LapackInt; e: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssyconv", dynlib: libName.}
proc ssyconvWork*(matrixLayout: cint; uplo: char; way: char; n: LapackInt; a: ptr cfloat;
                 lda: LapackInt; ipiv: ptr LapackInt; e: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssyconv_work", dynlib: libName.}
proc ssyequb*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat; lda: LapackInt;
             s: ptr cfloat; scond: ptr cfloat; amax: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssyequb", dynlib: libName.}
proc ssyequbWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                 lda: LapackInt; s: ptr cfloat; scond: ptr cfloat; amax: ptr cfloat;
                 work: ptr cfloat): LapackInt {.importc: "LAPACKE_ssyequb_work",
    dynlib: libName.}
proc ssyev*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; a: ptr cfloat;
           lda: LapackInt; w: ptr cfloat): LapackInt {.importc: "LAPACKE_ssyev",
    dynlib: libName.}
proc ssyevWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; a: ptr cfloat;
               lda: LapackInt; w: ptr cfloat; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ssyev_work", dynlib: libName.}
proc ssyevd*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; a: ptr cfloat;
            lda: LapackInt; w: ptr cfloat): LapackInt {.importc: "LAPACKE_ssyevd",
    dynlib: libName.}
proc ssyevdWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; w: ptr cfloat; work: ptr cfloat; lwork: LapackInt;
                iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ssyevd_work", dynlib: libName.}
proc ssyevr*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
            a: ptr cfloat; lda: LapackInt; vl: cfloat; vu: cfloat; il: LapackInt;
            iu: LapackInt; abstol: cfloat; m: ptr LapackInt; w: ptr cfloat; z: ptr cfloat;
            ldz: LapackInt; isuppz: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ssyevr", dynlib: libName.}
proc ssyevrWork*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                a: ptr cfloat; lda: LapackInt; vl: cfloat; vu: cfloat; il: LapackInt;
                iu: LapackInt; abstol: cfloat; m: ptr LapackInt; w: ptr cfloat;
                z: ptr cfloat; ldz: LapackInt; isuppz: ptr LapackInt; work: ptr cfloat;
                lwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ssyevr_work", dynlib: libName.}
proc ssyevx*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
            a: ptr cfloat; lda: LapackInt; vl: cfloat; vu: cfloat; il: LapackInt;
            iu: LapackInt; abstol: cfloat; m: ptr LapackInt; w: ptr cfloat; z: ptr cfloat;
            ldz: LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ssyevx", dynlib: libName.}
proc ssyevxWork*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                a: ptr cfloat; lda: LapackInt; vl: cfloat; vu: cfloat; il: LapackInt;
                iu: LapackInt; abstol: cfloat; m: ptr LapackInt; w: ptr cfloat;
                z: ptr cfloat; ldz: LapackInt; work: ptr cfloat; lwork: LapackInt;
                iwork: ptr LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ssyevx_work", dynlib: libName.}
proc ssygst*(matrixLayout: cint; itype: LapackInt; uplo: char; n: LapackInt;
            a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ssygst", dynlib: libName.}
proc ssygstWork*(matrixLayout: cint; itype: LapackInt; uplo: char; n: LapackInt;
                a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ssygst_work", dynlib: libName.}
proc ssygv*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char; n: LapackInt;
           a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt; w: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssygv", dynlib: libName.}
proc ssygvWork*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
               n: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
               w: ptr cfloat; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ssygv_work", dynlib: libName.}
proc ssygvd*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char; n: LapackInt;
            a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt; w: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssygvd", dynlib: libName.}
proc ssygvdWork*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
                n: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat;
                ldb: LapackInt; w: ptr cfloat; work: ptr cfloat; lwork: LapackInt;
                iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ssygvd_work", dynlib: libName.}
proc ssygvx*(matrixLayout: cint; itype: LapackInt; jobz: char; range: char; uplo: char;
            n: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
            vl: cfloat; vu: cfloat; il: LapackInt; iu: LapackInt; abstol: cfloat;
            m: ptr LapackInt; w: ptr cfloat; z: ptr cfloat; ldz: LapackInt;
            ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_ssygvx",
    dynlib: libName.}
proc ssygvxWork*(matrixLayout: cint; itype: LapackInt; jobz: char; range: char;
                uplo: char; n: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat;
                ldb: LapackInt; vl: cfloat; vu: cfloat; il: LapackInt; iu: LapackInt;
                abstol: cfloat; m: ptr LapackInt; w: ptr cfloat; z: ptr cfloat;
                ldz: LapackInt; work: ptr cfloat; lwork: LapackInt;
                iwork: ptr LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ssygvx_work", dynlib: libName.}
proc ssyrfs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr cfloat; lda: LapackInt; af: ptr cfloat; ldaf: LapackInt;
            ipiv: ptr LapackInt; b: ptr cfloat; ldb: LapackInt; x: ptr cfloat;
            ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssyrfs", dynlib: libName.}
proc ssyrfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr cfloat; lda: LapackInt; af: ptr cfloat; ldaf: LapackInt;
                ipiv: ptr LapackInt; b: ptr cfloat; ldb: LapackInt; x: ptr cfloat;
                ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat; work: ptr cfloat;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_ssyrfs_work",
    dynlib: libName.}
proc ssyrfsx*(matrixLayout: cint; uplo: char; equed: char; n: LapackInt; nrhs: LapackInt;
             a: ptr cfloat; lda: LapackInt; af: ptr cfloat; ldaf: LapackInt;
             ipiv: ptr LapackInt; s: ptr cfloat; b: ptr cfloat; ldb: LapackInt;
             x: ptr cfloat; ldx: LapackInt; rcond: ptr cfloat; berr: ptr cfloat;
             nErrBnds: LapackInt; errBndsNorm: ptr cfloat; errBndsComp: ptr cfloat;
             nparams: LapackInt; params: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssyrfsx", dynlib: libName.}
proc ssyrfsxWork*(matrixLayout: cint; uplo: char; equed: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr cfloat; lda: LapackInt; af: ptr cfloat;
                 ldaf: LapackInt; ipiv: ptr LapackInt; s: ptr cfloat; b: ptr cfloat;
                 ldb: LapackInt; x: ptr cfloat; ldx: LapackInt; rcond: ptr cfloat;
                 berr: ptr cfloat; nErrBnds: LapackInt; errBndsNorm: ptr cfloat;
                 errBndsComp: ptr cfloat; nparams: LapackInt; params: ptr cfloat;
                 work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ssyrfsx_work", dynlib: libName.}
proc ssysv*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt; a: ptr cfloat;
           lda: LapackInt; ipiv: ptr LapackInt; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ssysv", dynlib: libName.}
proc ssysvRook*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               a: ptr cfloat; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cfloat;
               ldb: LapackInt): LapackInt {.importc: "LAPACKE_ssysv_rook",
    dynlib: libName.}
proc ssysvRookWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                   a: ptr cfloat; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cfloat;
                   ldb: LapackInt; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ssysv_rook_work", dynlib: libName.}
proc ssysvWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               a: ptr cfloat; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cfloat;
               ldb: LapackInt; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ssysv_work", dynlib: libName.}
proc ssysvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr cfloat; lda: LapackInt; af: ptr cfloat; ldaf: LapackInt;
            ipiv: ptr LapackInt; b: ptr cfloat; ldb: LapackInt; x: ptr cfloat;
            ldx: LapackInt; rcond: ptr cfloat; ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssysvx", dynlib: libName.}
proc ssysvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; a: ptr cfloat; lda: LapackInt; af: ptr cfloat;
                ldaf: LapackInt; ipiv: ptr LapackInt; b: ptr cfloat; ldb: LapackInt;
                x: ptr cfloat; ldx: LapackInt; rcond: ptr cfloat; ferr: ptr cfloat;
                berr: ptr cfloat; work: ptr cfloat; lwork: LapackInt;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_ssysvx_work",
    dynlib: libName.}
proc ssysvxx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr cfloat; lda: LapackInt; af: ptr cfloat; ldaf: LapackInt;
             ipiv: ptr LapackInt; equed: cstring; s: ptr cfloat; b: ptr cfloat;
             ldb: LapackInt; x: ptr cfloat; ldx: LapackInt; rcond: ptr cfloat;
             rpvgrw: ptr cfloat; berr: ptr cfloat; nErrBnds: LapackInt;
             errBndsNorm: ptr cfloat; errBndsComp: ptr cfloat; nparams: LapackInt;
             params: ptr cfloat): LapackInt {.importc: "LAPACKE_ssysvxx",
    dynlib: libName.}
proc ssysvxxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr cfloat; lda: LapackInt; af: ptr cfloat;
                 ldaf: LapackInt; ipiv: ptr LapackInt; equed: cstring; s: ptr cfloat;
                 b: ptr cfloat; ldb: LapackInt; x: ptr cfloat; ldx: LapackInt;
                 rcond: ptr cfloat; rpvgrw: ptr cfloat; berr: ptr cfloat;
                 nErrBnds: LapackInt; errBndsNorm: ptr cfloat;
                 errBndsComp: ptr cfloat; nparams: LapackInt; params: ptr cfloat;
                 work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ssysvxx_work", dynlib: libName.}
proc ssyswapr*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
              lda: LapackInt; i1: LapackInt; i2: LapackInt): LapackInt {.
    importc: "LAPACKE_ssyswapr", dynlib: libName.}
proc ssyswaprWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                  lda: LapackInt; i1: LapackInt; i2: LapackInt): LapackInt {.
    importc: "LAPACKE_ssyswapr_work", dynlib: libName.}
proc ssytrd*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat; lda: LapackInt;
            d: ptr cfloat; e: ptr cfloat; tau: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssytrd", dynlib: libName.}
proc ssytrdWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; d: ptr cfloat; e: ptr cfloat; tau: ptr cfloat;
                work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ssytrd_work", dynlib: libName.}
proc ssytrf*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat; lda: LapackInt;
            ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_ssytrf",
    dynlib: libName.}
proc ssytrfRook*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ssytrf_rook", dynlib: libName.}
proc ssytrfRookWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                    lda: LapackInt; ipiv: ptr LapackInt; work: ptr cfloat;
                    lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ssytrf_rook_work", dynlib: libName.}
proc ssytrfWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; ipiv: ptr LapackInt; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ssytrf_work", dynlib: libName.}
proc ssytri*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat; lda: LapackInt;
            ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_ssytri",
    dynlib: libName.}
proc ssytri2*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat; lda: LapackInt;
             ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_ssytri2",
    dynlib: libName.}
proc ssytri2Work*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                 lda: LapackInt; ipiv: ptr LapackInt; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ssytri2_work", dynlib: libName.}
proc ssytri2x*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
              lda: LapackInt; ipiv: ptr LapackInt; nb: LapackInt): LapackInt {.
    importc: "LAPACKE_ssytri2x", dynlib: libName.}
proc ssytri2xWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                  lda: LapackInt; ipiv: ptr LapackInt; work: ptr cfloat; nb: LapackInt): LapackInt {.
    importc: "LAPACKE_ssytri2x_work", dynlib: libName.}
proc ssytriWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; ipiv: ptr LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssytri_work", dynlib: libName.}
proc ssytrs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr cfloat; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cfloat;
            ldb: LapackInt): LapackInt {.importc: "LAPACKE_ssytrs", dynlib: libName.}
proc ssytrs2*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr cfloat; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cfloat;
             ldb: LapackInt): LapackInt {.importc: "LAPACKE_ssytrs2", dynlib: libName.}
proc ssytrs2Work*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr cfloat; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cfloat;
                 ldb: LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssytrs2_work", dynlib: libName.}
proc ssytrsRook*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr cfloat; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cfloat;
                ldb: LapackInt): LapackInt {.importc: "LAPACKE_ssytrs_rook",
    dynlib: libName.}
proc ssytrsRookWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                    a: ptr cfloat; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cfloat;
                    ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ssytrs_rook_work", dynlib: libName.}
proc ssytrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr cfloat; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cfloat;
                ldb: LapackInt): LapackInt {.importc: "LAPACKE_ssytrs_work",
    dynlib: libName.}
proc stbcon*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
            kd: LapackInt; ab: ptr cfloat; ldab: LapackInt; rcond: ptr cfloat): LapackInt {.
    importc: "LAPACKE_stbcon", dynlib: libName.}
proc stbconWork*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
                kd: LapackInt; ab: ptr cfloat; ldab: LapackInt; rcond: ptr cfloat;
                work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_stbcon_work", dynlib: libName.}
proc stbrfs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            kd: LapackInt; nrhs: LapackInt; ab: ptr cfloat; ldab: LapackInt;
            b: ptr cfloat; ldb: LapackInt; x: ptr cfloat; ldx: LapackInt;
            ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.importc: "LAPACKE_stbrfs",
    dynlib: libName.}
proc stbrfsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                kd: LapackInt; nrhs: LapackInt; ab: ptr cfloat; ldab: LapackInt;
                b: ptr cfloat; ldb: LapackInt; x: ptr cfloat; ldx: LapackInt;
                ferr: ptr cfloat; berr: ptr cfloat; work: ptr cfloat;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_stbrfs_work",
    dynlib: libName.}
proc stbtrs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            kd: LapackInt; nrhs: LapackInt; ab: ptr cfloat; ldab: LapackInt;
            b: ptr cfloat; ldb: LapackInt): LapackInt {.importc: "LAPACKE_stbtrs",
    dynlib: libName.}
proc stbtrsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                kd: LapackInt; nrhs: LapackInt; ab: ptr cfloat; ldab: LapackInt;
                b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_stbtrs_work", dynlib: libName.}
proc stfsm*(matrixLayout: cint; transr: char; side: char; uplo: char; trans: char;
           diag: char; m: LapackInt; n: LapackInt; alpha: cfloat; a: ptr cfloat;
           b: ptr cfloat; ldb: LapackInt): LapackInt {.importc: "LAPACKE_stfsm",
    dynlib: libName.}
proc stfsmWork*(matrixLayout: cint; transr: char; side: char; uplo: char; trans: char;
               diag: char; m: LapackInt; n: LapackInt; alpha: cfloat; a: ptr cfloat;
               b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_stfsm_work", dynlib: libName.}
proc stftri*(matrixLayout: cint; transr: char; uplo: char; diag: char; n: LapackInt;
            a: ptr cfloat): LapackInt {.importc: "LAPACKE_stftri", dynlib: libName.}
proc stftriWork*(matrixLayout: cint; transr: char; uplo: char; diag: char; n: LapackInt;
                a: ptr cfloat): LapackInt {.importc: "LAPACKE_stftri_work",
                                        dynlib: libName.}
proc stfttp*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt; arf: ptr cfloat;
            ap: ptr cfloat): LapackInt {.importc: "LAPACKE_stfttp", dynlib: libName.}
proc stfttpWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                arf: ptr cfloat; ap: ptr cfloat): LapackInt {.
    importc: "LAPACKE_stfttp_work", dynlib: libName.}
proc stfttr*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt; arf: ptr cfloat;
            a: ptr cfloat; lda: LapackInt): LapackInt {.importc: "LAPACKE_stfttr",
    dynlib: libName.}
proc stfttrWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                arf: ptr cfloat; a: ptr cfloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_stfttr_work", dynlib: libName.}
proc stgevc*(matrixLayout: cint; side: char; howmny: char; select: ptr LapackLogical;
            n: LapackInt; s: ptr cfloat; lds: LapackInt; p: ptr cfloat; ldp: LapackInt;
            vl: ptr cfloat; ldvl: LapackInt; vr: ptr cfloat; ldvr: LapackInt;
            mm: LapackInt; m: ptr LapackInt): LapackInt {.importc: "LAPACKE_stgevc",
    dynlib: libName.}
proc stgevcWork*(matrixLayout: cint; side: char; howmny: char;
                select: ptr LapackLogical; n: LapackInt; s: ptr cfloat; lds: LapackInt;
                p: ptr cfloat; ldp: LapackInt; vl: ptr cfloat; ldvl: LapackInt;
                vr: ptr cfloat; ldvr: LapackInt; mm: LapackInt; m: ptr LapackInt;
                work: ptr cfloat): LapackInt {.importc: "LAPACKE_stgevc_work",
    dynlib: libName.}
proc stgexc*(matrixLayout: cint; wantq: LapackLogical; wantz: LapackLogical;
            n: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
            q: ptr cfloat; ldq: LapackInt; z: ptr cfloat; ldz: LapackInt;
            ifst: ptr LapackInt; ilst: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_stgexc", dynlib: libName.}
proc stgexcWork*(matrixLayout: cint; wantq: LapackLogical; wantz: LapackLogical;
                n: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat;
                ldb: LapackInt; q: ptr cfloat; ldq: LapackInt; z: ptr cfloat;
                ldz: LapackInt; ifst: ptr LapackInt; ilst: ptr LapackInt;
                work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_stgexc_work", dynlib: libName.}
proc stgsen*(matrixLayout: cint; ijob: LapackInt; wantq: LapackLogical;
            wantz: LapackLogical; select: ptr LapackLogical; n: LapackInt;
            a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
            alphar: ptr cfloat; alphai: ptr cfloat; beta: ptr cfloat; q: ptr cfloat;
            ldq: LapackInt; z: ptr cfloat; ldz: LapackInt; m: ptr LapackInt;
            pl: ptr cfloat; pr: ptr cfloat; dif: ptr cfloat): LapackInt {.
    importc: "LAPACKE_stgsen", dynlib: libName.}
proc stgsenWork*(matrixLayout: cint; ijob: LapackInt; wantq: LapackLogical;
                wantz: LapackLogical; select: ptr LapackLogical; n: LapackInt;
                a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
                alphar: ptr cfloat; alphai: ptr cfloat; beta: ptr cfloat; q: ptr cfloat;
                ldq: LapackInt; z: ptr cfloat; ldz: LapackInt; m: ptr LapackInt;
                pl: ptr cfloat; pr: ptr cfloat; dif: ptr cfloat; work: ptr cfloat;
                lwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_stgsen_work", dynlib: libName.}
proc stgsja*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
            p: LapackInt; n: LapackInt; k: LapackInt; l: LapackInt; a: ptr cfloat;
            lda: LapackInt; b: ptr cfloat; ldb: LapackInt; tola: cfloat; tolb: cfloat;
            alpha: ptr cfloat; beta: ptr cfloat; u: ptr cfloat; ldu: LapackInt;
            v: ptr cfloat; ldv: LapackInt; q: ptr cfloat; ldq: LapackInt;
            ncycle: ptr LapackInt): LapackInt {.importc: "LAPACKE_stgsja",
    dynlib: libName.}
proc stgsjaWork*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
                p: LapackInt; n: LapackInt; k: LapackInt; l: LapackInt; a: ptr cfloat;
                lda: LapackInt; b: ptr cfloat; ldb: LapackInt; tola: cfloat;
                tolb: cfloat; alpha: ptr cfloat; beta: ptr cfloat; u: ptr cfloat;
                ldu: LapackInt; v: ptr cfloat; ldv: LapackInt; q: ptr cfloat;
                ldq: LapackInt; work: ptr cfloat; ncycle: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_stgsja_work", dynlib: libName.}
proc stgsna*(matrixLayout: cint; job: char; howmny: char; select: ptr LapackLogical;
            n: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
            vl: ptr cfloat; ldvl: LapackInt; vr: ptr cfloat; ldvr: LapackInt;
            s: ptr cfloat; dif: ptr cfloat; mm: LapackInt; m: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_stgsna", dynlib: libName.}
proc stgsnaWork*(matrixLayout: cint; job: char; howmny: char;
                select: ptr LapackLogical; n: LapackInt; a: ptr cfloat; lda: LapackInt;
                b: ptr cfloat; ldb: LapackInt; vl: ptr cfloat; ldvl: LapackInt;
                vr: ptr cfloat; ldvr: LapackInt; s: ptr cfloat; dif: ptr cfloat;
                mm: LapackInt; m: ptr LapackInt; work: ptr cfloat; lwork: LapackInt;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_stgsna_work",
    dynlib: libName.}
proc stgsyl*(matrixLayout: cint; trans: char; ijob: LapackInt; m: LapackInt;
            n: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
            c: ptr cfloat; ldc: LapackInt; d: ptr cfloat; ldd: LapackInt; e: ptr cfloat;
            lde: LapackInt; f: ptr cfloat; ldf: LapackInt; scale: ptr cfloat;
            dif: ptr cfloat): LapackInt {.importc: "LAPACKE_stgsyl", dynlib: libName.}
proc stgsylWork*(matrixLayout: cint; trans: char; ijob: LapackInt; m: LapackInt;
                n: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat;
                ldb: LapackInt; c: ptr cfloat; ldc: LapackInt; d: ptr cfloat;
                ldd: LapackInt; e: ptr cfloat; lde: LapackInt; f: ptr cfloat;
                ldf: LapackInt; scale: ptr cfloat; dif: ptr cfloat; work: ptr cfloat;
                lwork: LapackInt; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_stgsyl_work", dynlib: libName.}
proc stpcon*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
            ap: ptr cfloat; rcond: ptr cfloat): LapackInt {.importc: "LAPACKE_stpcon",
    dynlib: libName.}
proc stpconWork*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
                ap: ptr cfloat; rcond: ptr cfloat; work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_stpcon_work", dynlib: libName.}
proc stpmqrt*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
             k: LapackInt; l: LapackInt; nb: LapackInt; v: ptr cfloat; ldv: LapackInt;
             t: ptr cfloat; ldt: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat;
             ldb: LapackInt): LapackInt {.importc: "LAPACKE_stpmqrt", dynlib: libName.}
proc stpmqrtWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt;
                 n: LapackInt; k: LapackInt; l: LapackInt; nb: LapackInt; v: ptr cfloat;
                 ldv: LapackInt; t: ptr cfloat; ldt: LapackInt; a: ptr cfloat;
                 lda: LapackInt; b: ptr cfloat; ldb: LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_stpmqrt_work", dynlib: libName.}
proc stpqrt*(matrixLayout: cint; m: LapackInt; n: LapackInt; l: LapackInt; nb: LapackInt;
            a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt; t: ptr cfloat;
            ldt: LapackInt): LapackInt {.importc: "LAPACKE_stpqrt", dynlib: libName.}
proc stpqrt2*(matrixLayout: cint; m: LapackInt; n: LapackInt; l: LapackInt;
             a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt; t: ptr cfloat;
             ldt: LapackInt): LapackInt {.importc: "LAPACKE_stpqrt2", dynlib: libName.}
proc stpqrt2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt; l: LapackInt;
                 a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
                 t: ptr cfloat; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_stpqrt2_work", dynlib: libName.}
proc stpqrtWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; l: LapackInt;
                nb: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat;
                ldb: LapackInt; t: ptr cfloat; ldt: LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_stpqrt_work", dynlib: libName.}
proc stprfb*(matrixLayout: cint; side: char; trans: char; direct: char; storev: char;
            m: LapackInt; n: LapackInt; k: LapackInt; l: LapackInt; v: ptr cfloat;
            ldv: LapackInt; t: ptr cfloat; ldt: LapackInt; a: ptr cfloat; lda: LapackInt;
            b: ptr cfloat; ldb: LapackInt): LapackInt {.importc: "LAPACKE_stprfb",
    dynlib: libName.}
proc stprfbWork*(matrixLayout: cint; side: char; trans: char; direct: char; storev: char;
                m: LapackInt; n: LapackInt; k: LapackInt; l: LapackInt; v: ptr cfloat;
                ldv: LapackInt; t: ptr cfloat; ldt: LapackInt; a: ptr cfloat;
                lda: LapackInt; b: ptr cfloat; ldb: LapackInt; work: ptr cfloat;
                ldwork: LapackInt): LapackInt {.importc: "LAPACKE_stprfb_work",
    dynlib: libName.}
proc stprfs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            nrhs: LapackInt; ap: ptr cfloat; b: ptr cfloat; ldb: LapackInt; x: ptr cfloat;
            ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.
    importc: "LAPACKE_stprfs", dynlib: libName.}
proc stprfsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                nrhs: LapackInt; ap: ptr cfloat; b: ptr cfloat; ldb: LapackInt;
                x: ptr cfloat; ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat;
                work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_stprfs_work", dynlib: libName.}
proc stptri*(matrixLayout: cint; uplo: char; diag: char; n: LapackInt; ap: ptr cfloat): LapackInt {.
    importc: "LAPACKE_stptri", dynlib: libName.}
proc stptriWork*(matrixLayout: cint; uplo: char; diag: char; n: LapackInt; ap: ptr cfloat): LapackInt {.
    importc: "LAPACKE_stptri_work", dynlib: libName.}
proc stptrs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            nrhs: LapackInt; ap: ptr cfloat; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_stptrs", dynlib: libName.}
proc stptrsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                nrhs: LapackInt; ap: ptr cfloat; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_stptrs_work", dynlib: libName.}
proc stpttf*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt; ap: ptr cfloat;
            arf: ptr cfloat): LapackInt {.importc: "LAPACKE_stpttf", dynlib: libName.}
proc stpttfWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                ap: ptr cfloat; arf: ptr cfloat): LapackInt {.
    importc: "LAPACKE_stpttf_work", dynlib: libName.}
proc stpttr*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cfloat; a: ptr cfloat;
            lda: LapackInt): LapackInt {.importc: "LAPACKE_stpttr", dynlib: libName.}
proc stpttrWork*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr cfloat;
                a: ptr cfloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_stpttr_work", dynlib: libName.}
proc strcon*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
            a: ptr cfloat; lda: LapackInt; rcond: ptr cfloat): LapackInt {.
    importc: "LAPACKE_strcon", dynlib: libName.}
proc strconWork*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
                a: ptr cfloat; lda: LapackInt; rcond: ptr cfloat; work: ptr cfloat;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_strcon_work",
    dynlib: libName.}
proc strevc*(matrixLayout: cint; side: char; howmny: char; select: ptr LapackLogical;
            n: LapackInt; t: ptr cfloat; ldt: LapackInt; vl: ptr cfloat; ldvl: LapackInt;
            vr: ptr cfloat; ldvr: LapackInt; mm: LapackInt; m: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_strevc", dynlib: libName.}
proc strevcWork*(matrixLayout: cint; side: char; howmny: char;
                select: ptr LapackLogical; n: LapackInt; t: ptr cfloat; ldt: LapackInt;
                vl: ptr cfloat; ldvl: LapackInt; vr: ptr cfloat; ldvr: LapackInt;
                mm: LapackInt; m: ptr LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_strevc_work", dynlib: libName.}
proc strexc*(matrixLayout: cint; compq: char; n: LapackInt; t: ptr cfloat; ldt: LapackInt;
            q: ptr cfloat; ldq: LapackInt; ifst: ptr LapackInt; ilst: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_strexc", dynlib: libName.}
proc strexcWork*(matrixLayout: cint; compq: char; n: LapackInt; t: ptr cfloat;
                ldt: LapackInt; q: ptr cfloat; ldq: LapackInt; ifst: ptr LapackInt;
                ilst: ptr LapackInt; work: ptr cfloat): LapackInt {.
    importc: "LAPACKE_strexc_work", dynlib: libName.}
proc strrfs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            nrhs: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
            x: ptr cfloat; ldx: LapackInt; ferr: ptr cfloat; berr: ptr cfloat): LapackInt {.
    importc: "LAPACKE_strrfs", dynlib: libName.}
proc strrfsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                nrhs: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat;
                ldb: LapackInt; x: ptr cfloat; ldx: LapackInt; ferr: ptr cfloat;
                berr: ptr cfloat; work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_strrfs_work", dynlib: libName.}
proc strsen*(matrixLayout: cint; job: char; compq: char; select: ptr LapackLogical;
            n: LapackInt; t: ptr cfloat; ldt: LapackInt; q: ptr cfloat; ldq: LapackInt;
            wr: ptr cfloat; wi: ptr cfloat; m: ptr LapackInt; s: ptr cfloat; sep: ptr cfloat): LapackInt {.
    importc: "LAPACKE_strsen", dynlib: libName.}
proc strsenWork*(matrixLayout: cint; job: char; compq: char; select: ptr LapackLogical;
                n: LapackInt; t: ptr cfloat; ldt: LapackInt; q: ptr cfloat;
                ldq: LapackInt; wr: ptr cfloat; wi: ptr cfloat; m: ptr LapackInt;
                s: ptr cfloat; sep: ptr cfloat; work: ptr cfloat; lwork: LapackInt;
                iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_strsen_work", dynlib: libName.}
proc strsna*(matrixLayout: cint; job: char; howmny: char; select: ptr LapackLogical;
            n: LapackInt; t: ptr cfloat; ldt: LapackInt; vl: ptr cfloat; ldvl: LapackInt;
            vr: ptr cfloat; ldvr: LapackInt; s: ptr cfloat; sep: ptr cfloat; mm: LapackInt;
            m: ptr LapackInt): LapackInt {.importc: "LAPACKE_strsna", dynlib: libName.}
proc strsnaWork*(matrixLayout: cint; job: char; howmny: char;
                select: ptr LapackLogical; n: LapackInt; t: ptr cfloat; ldt: LapackInt;
                vl: ptr cfloat; ldvl: LapackInt; vr: ptr cfloat; ldvr: LapackInt;
                s: ptr cfloat; sep: ptr cfloat; mm: LapackInt; m: ptr LapackInt;
                work: ptr cfloat; ldwork: LapackInt; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_strsna_work", dynlib: libName.}
proc strsyl*(matrixLayout: cint; trana: char; tranb: char; isgn: LapackInt; m: LapackInt;
            n: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt;
            c: ptr cfloat; ldc: LapackInt; scale: ptr cfloat): LapackInt {.
    importc: "LAPACKE_strsyl", dynlib: libName.}
proc strsylWork*(matrixLayout: cint; trana: char; tranb: char; isgn: LapackInt;
                m: LapackInt; n: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat;
                ldb: LapackInt; c: ptr cfloat; ldc: LapackInt; scale: ptr cfloat): LapackInt {.
    importc: "LAPACKE_strsyl_work", dynlib: libName.}
proc strtri*(matrixLayout: cint; uplo: char; diag: char; n: LapackInt; a: ptr cfloat;
            lda: LapackInt): LapackInt {.importc: "LAPACKE_strtri", dynlib: libName.}
proc strtriWork*(matrixLayout: cint; uplo: char; diag: char; n: LapackInt; a: ptr cfloat;
                lda: LapackInt): LapackInt {.importc: "LAPACKE_strtri_work",
    dynlib: libName.}
proc strtrs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            nrhs: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_strtrs", dynlib: libName.}
proc strtrsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                nrhs: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat;
                ldb: LapackInt): LapackInt {.importc: "LAPACKE_strtrs_work",
    dynlib: libName.}
proc strttf*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt; a: ptr cfloat;
            lda: LapackInt; arf: ptr cfloat): LapackInt {.importc: "LAPACKE_strttf",
    dynlib: libName.}
proc strttfWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                a: ptr cfloat; lda: LapackInt; arf: ptr cfloat): LapackInt {.
    importc: "LAPACKE_strttf_work", dynlib: libName.}
proc strttp*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat; lda: LapackInt;
            ap: ptr cfloat): LapackInt {.importc: "LAPACKE_strttp", dynlib: libName.}
proc strttpWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; ap: ptr cfloat): LapackInt {.
    importc: "LAPACKE_strttp_work", dynlib: libName.}
proc stzrzf*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
            lda: LapackInt; tau: ptr cfloat): LapackInt {.importc: "LAPACKE_stzrzf",
    dynlib: libName.}
proc stzrzfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
                lda: LapackInt; tau: ptr cfloat; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_stzrzf_work", dynlib: libName.}
proc zbbcsd*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char; jobv2t: char;
            trans: char; m: LapackInt; p: LapackInt; q: LapackInt; theta: ptr cdouble;
            phi: ptr cdouble; u1: ptr LapackComplexDouble; ldu1: LapackInt;
            u2: ptr LapackComplexDouble; ldu2: LapackInt;
            v1t: ptr LapackComplexDouble; ldv1t: LapackInt;
            v2t: ptr LapackComplexDouble; ldv2t: LapackInt; b11d: ptr cdouble;
            b11e: ptr cdouble; b12d: ptr cdouble; b12e: ptr cdouble; b21d: ptr cdouble;
            b21e: ptr cdouble; b22d: ptr cdouble; b22e: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zbbcsd", dynlib: libName.}
proc zbbcsdWork*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char;
                jobv2t: char; trans: char; m: LapackInt; p: LapackInt; q: LapackInt;
                theta: ptr cdouble; phi: ptr cdouble; u1: ptr LapackComplexDouble;
                ldu1: LapackInt; u2: ptr LapackComplexDouble; ldu2: LapackInt;
                v1t: ptr LapackComplexDouble; ldv1t: LapackInt;
                v2t: ptr LapackComplexDouble; ldv2t: LapackInt; b11d: ptr cdouble;
                b11e: ptr cdouble; b12d: ptr cdouble; b12e: ptr cdouble;
                b21d: ptr cdouble; b21e: ptr cdouble; b22d: ptr cdouble;
                b22e: ptr cdouble; rwork: ptr cdouble; lrwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zbbcsd_work", dynlib: libName.}
proc zbdsqr*(matrixLayout: cint; uplo: char; n: LapackInt; ncvt: LapackInt;
            nru: LapackInt; ncc: LapackInt; d: ptr cdouble; e: ptr cdouble;
            vt: ptr LapackComplexDouble; ldvt: LapackInt; u: ptr LapackComplexDouble;
            ldu: LapackInt; c: ptr LapackComplexDouble; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_zbdsqr", dynlib: libName.}
proc zbdsqrWork*(matrixLayout: cint; uplo: char; n: LapackInt; ncvt: LapackInt;
                nru: LapackInt; ncc: LapackInt; d: ptr cdouble; e: ptr cdouble;
                vt: ptr LapackComplexDouble; ldvt: LapackInt;
                u: ptr LapackComplexDouble; ldu: LapackInt;
                c: ptr LapackComplexDouble; ldc: LapackInt; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zbdsqr_work", dynlib: libName.}
proc zcgesv*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
            b: ptr LapackComplexDouble; ldb: LapackInt; x: ptr LapackComplexDouble;
            ldx: LapackInt; iter: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zcgesv", dynlib: libName.}
proc zcgesvWork*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt;
                work: ptr LapackComplexDouble; swork: ptr LapackComplexFloat;
                rwork: ptr cdouble; iter: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zcgesv_work", dynlib: libName.}
proc zcposv*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; b: ptr LapackComplexDouble;
            ldb: LapackInt; x: ptr LapackComplexDouble; ldx: LapackInt;
            iter: ptr LapackInt): LapackInt {.importc: "LAPACKE_zcposv",
    dynlib: libName.}
proc zcposvWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt;
                work: ptr LapackComplexDouble; swork: ptr LapackComplexFloat;
                rwork: ptr cdouble; iter: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zcposv_work", dynlib: libName.}
proc zgbbrd*(matrixLayout: cint; vect: char; m: LapackInt; n: LapackInt; ncc: LapackInt;
            kl: LapackInt; ku: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
            d: ptr cdouble; e: ptr cdouble; q: ptr LapackComplexDouble; ldq: LapackInt;
            pt: ptr LapackComplexDouble; ldpt: LapackInt; c: ptr LapackComplexDouble;
            ldc: LapackInt): LapackInt {.importc: "LAPACKE_zgbbrd", dynlib: libName.}
proc zgbbrdWork*(matrixLayout: cint; vect: char; m: LapackInt; n: LapackInt;
                ncc: LapackInt; kl: LapackInt; ku: LapackInt;
                ab: ptr LapackComplexDouble; ldab: LapackInt; d: ptr cdouble;
                e: ptr cdouble; q: ptr LapackComplexDouble; ldq: LapackInt;
                pt: ptr LapackComplexDouble; ldpt: LapackInt;
                c: ptr LapackComplexDouble; ldc: LapackInt;
                work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgbbrd_work", dynlib: libName.}
proc zgbcon*(matrixLayout: cint; norm: char; n: LapackInt; kl: LapackInt; ku: LapackInt;
            ab: ptr LapackComplexDouble; ldab: LapackInt; ipiv: ptr LapackInt;
            anorm: cdouble; rcond: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgbcon", dynlib: libName.}
proc zgbconWork*(matrixLayout: cint; norm: char; n: LapackInt; kl: LapackInt;
                ku: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
                ipiv: ptr LapackInt; anorm: cdouble; rcond: ptr cdouble;
                work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgbcon_work", dynlib: libName.}
proc zgbequ*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
            ku: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
            r: ptr cdouble; c: ptr cdouble; rowcnd: ptr cdouble; colcnd: ptr cdouble;
            amax: ptr cdouble): LapackInt {.importc: "LAPACKE_zgbequ", dynlib: libName.}
proc zgbequWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
                ku: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
                r: ptr cdouble; c: ptr cdouble; rowcnd: ptr cdouble; colcnd: ptr cdouble;
                amax: ptr cdouble): LapackInt {.importc: "LAPACKE_zgbequ_work",
    dynlib: libName.}
proc zgbequb*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
             ku: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
             r: ptr cdouble; c: ptr cdouble; rowcnd: ptr cdouble; colcnd: ptr cdouble;
             amax: ptr cdouble): LapackInt {.importc: "LAPACKE_zgbequb",
    dynlib: libName.}
proc zgbequbWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
                 ku: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
                 r: ptr cdouble; c: ptr cdouble; rowcnd: ptr cdouble;
                 colcnd: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgbequb_work", dynlib: libName.}
proc zgbrfs*(matrixLayout: cint; trans: char; n: LapackInt; kl: LapackInt; ku: LapackInt;
            nrhs: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
            afb: ptr LapackComplexDouble; ldafb: LapackInt; ipiv: ptr LapackInt;
            b: ptr LapackComplexDouble; ldb: LapackInt; x: ptr LapackComplexDouble;
            ldx: LapackInt; ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgbrfs", dynlib: libName.}
proc zgbrfsWork*(matrixLayout: cint; trans: char; n: LapackInt; kl: LapackInt;
                ku: LapackInt; nrhs: LapackInt; ab: ptr LapackComplexDouble;
                ldab: LapackInt; afb: ptr LapackComplexDouble; ldafb: LapackInt;
                ipiv: ptr LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; ferr: ptr cdouble;
                berr: ptr cdouble; work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgbrfs_work", dynlib: libName.}
proc zgbrfsx*(matrixLayout: cint; trans: char; equed: char; n: LapackInt; kl: LapackInt;
             ku: LapackInt; nrhs: LapackInt; ab: ptr LapackComplexDouble;
             ldab: LapackInt; afb: ptr LapackComplexDouble; ldafb: LapackInt;
             ipiv: ptr LapackInt; r: ptr cdouble; c: ptr cdouble;
             b: ptr LapackComplexDouble; ldb: LapackInt; x: ptr LapackComplexDouble;
             ldx: LapackInt; rcond: ptr cdouble; berr: ptr cdouble; nErrBnds: LapackInt;
             errBndsNorm: ptr cdouble; errBndsComp: ptr cdouble; nparams: LapackInt;
             params: ptr cdouble): LapackInt {.importc: "LAPACKE_zgbrfsx",
    dynlib: libName.}
proc zgbrfsxWork*(matrixLayout: cint; trans: char; equed: char; n: LapackInt;
                 kl: LapackInt; ku: LapackInt; nrhs: LapackInt;
                 ab: ptr LapackComplexDouble; ldab: LapackInt;
                 afb: ptr LapackComplexDouble; ldafb: LapackInt; ipiv: ptr LapackInt;
                 r: ptr cdouble; c: ptr cdouble; b: ptr LapackComplexDouble;
                 ldb: LapackInt; x: ptr LapackComplexDouble; ldx: LapackInt;
                 rcond: ptr cdouble; berr: ptr cdouble; nErrBnds: LapackInt;
                 errBndsNorm: ptr cdouble; errBndsComp: ptr cdouble;
                 nparams: LapackInt; params: ptr cdouble;
                 work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgbrfsx_work", dynlib: libName.}
proc zgbsv*(matrixLayout: cint; n: LapackInt; kl: LapackInt; ku: LapackInt;
           nrhs: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
           ipiv: ptr LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zgbsv", dynlib: libName.}
proc zgbsvWork*(matrixLayout: cint; n: LapackInt; kl: LapackInt; ku: LapackInt;
               nrhs: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
               ipiv: ptr LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zgbsv_work", dynlib: libName.}
proc zgbsvx*(matrixLayout: cint; fact: char; trans: char; n: LapackInt; kl: LapackInt;
            ku: LapackInt; nrhs: LapackInt; ab: ptr LapackComplexDouble;
            ldab: LapackInt; afb: ptr LapackComplexDouble; ldafb: LapackInt;
            ipiv: ptr LapackInt; equed: cstring; r: ptr cdouble; c: ptr cdouble;
            b: ptr LapackComplexDouble; ldb: LapackInt; x: ptr LapackComplexDouble;
            ldx: LapackInt; rcond: ptr cdouble; ferr: ptr cdouble; berr: ptr cdouble;
            rpivot: ptr cdouble): LapackInt {.importc: "LAPACKE_zgbsvx",
    dynlib: libName.}
proc zgbsvxWork*(matrixLayout: cint; fact: char; trans: char; n: LapackInt;
                kl: LapackInt; ku: LapackInt; nrhs: LapackInt;
                ab: ptr LapackComplexDouble; ldab: LapackInt;
                afb: ptr LapackComplexDouble; ldafb: LapackInt; ipiv: ptr LapackInt;
                equed: cstring; r: ptr cdouble; c: ptr cdouble;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
                ferr: ptr cdouble; berr: ptr cdouble; work: ptr LapackComplexDouble;
                rwork: ptr cdouble): LapackInt {.importc: "LAPACKE_zgbsvx_work",
    dynlib: libName.}
proc zgbsvxx*(matrixLayout: cint; fact: char; trans: char; n: LapackInt; kl: LapackInt;
             ku: LapackInt; nrhs: LapackInt; ab: ptr LapackComplexDouble;
             ldab: LapackInt; afb: ptr LapackComplexDouble; ldafb: LapackInt;
             ipiv: ptr LapackInt; equed: cstring; r: ptr cdouble; c: ptr cdouble;
             b: ptr LapackComplexDouble; ldb: LapackInt; x: ptr LapackComplexDouble;
             ldx: LapackInt; rcond: ptr cdouble; rpvgrw: ptr cdouble; berr: ptr cdouble;
             nErrBnds: LapackInt; errBndsNorm: ptr cdouble; errBndsComp: ptr cdouble;
             nparams: LapackInt; params: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgbsvxx", dynlib: libName.}
proc zgbsvxxWork*(matrixLayout: cint; fact: char; trans: char; n: LapackInt;
                 kl: LapackInt; ku: LapackInt; nrhs: LapackInt;
                 ab: ptr LapackComplexDouble; ldab: LapackInt;
                 afb: ptr LapackComplexDouble; ldafb: LapackInt; ipiv: ptr LapackInt;
                 equed: cstring; r: ptr cdouble; c: ptr cdouble;
                 b: ptr LapackComplexDouble; ldb: LapackInt;
                 x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
                 rpvgrw: ptr cdouble; berr: ptr cdouble; nErrBnds: LapackInt;
                 errBndsNorm: ptr cdouble; errBndsComp: ptr cdouble;
                 nparams: LapackInt; params: ptr cdouble;
                 work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgbsvxx_work", dynlib: libName.}
proc zgbtrf*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
            ku: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
            ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_zgbtrf",
    dynlib: libName.}
proc zgbtrfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
                ku: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
                ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_zgbtrf_work",
    dynlib: libName.}
proc zgbtrs*(matrixLayout: cint; trans: char; n: LapackInt; kl: LapackInt; ku: LapackInt;
            nrhs: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
            ipiv: ptr LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zgbtrs", dynlib: libName.}
proc zgbtrsWork*(matrixLayout: cint; trans: char; n: LapackInt; kl: LapackInt;
                ku: LapackInt; nrhs: LapackInt; ab: ptr LapackComplexDouble;
                ldab: LapackInt; ipiv: ptr LapackInt; b: ptr LapackComplexDouble;
                ldb: LapackInt): LapackInt {.importc: "LAPACKE_zgbtrs_work",
    dynlib: libName.}
proc zgebak*(matrixLayout: cint; job: char; side: char; n: LapackInt; ilo: LapackInt;
            ihi: LapackInt; scale: ptr cdouble; m: LapackInt;
            v: ptr LapackComplexDouble; ldv: LapackInt): LapackInt {.
    importc: "LAPACKE_zgebak", dynlib: libName.}
proc zgebakWork*(matrixLayout: cint; job: char; side: char; n: LapackInt; ilo: LapackInt;
                ihi: LapackInt; scale: ptr cdouble; m: LapackInt;
                v: ptr LapackComplexDouble; ldv: LapackInt): LapackInt {.
    importc: "LAPACKE_zgebak_work", dynlib: libName.}
proc zgebal*(matrixLayout: cint; job: char; n: LapackInt; a: ptr LapackComplexDouble;
            lda: LapackInt; ilo: ptr LapackInt; ihi: ptr LapackInt; scale: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgebal", dynlib: libName.}
proc zgebalWork*(matrixLayout: cint; job: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; ilo: ptr LapackInt;
                ihi: ptr LapackInt; scale: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgebal_work", dynlib: libName.}
proc zgebrd*(matrixLayout: cint; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; d: ptr cdouble; e: ptr cdouble;
            tauq: ptr LapackComplexDouble; taup: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zgebrd", dynlib: libName.}
proc zgebrdWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; d: ptr cdouble;
                e: ptr cdouble; tauq: ptr LapackComplexDouble;
                taup: ptr LapackComplexDouble; work: ptr LapackComplexDouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_zgebrd_work",
    dynlib: libName.}
proc zgecon*(matrixLayout: cint; norm: char; n: LapackInt; a: ptr LapackComplexDouble;
            lda: LapackInt; anorm: cdouble; rcond: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgecon", dynlib: libName.}
proc zgeconWork*(matrixLayout: cint; norm: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; anorm: cdouble;
                rcond: ptr cdouble; work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgecon_work", dynlib: libName.}
proc zgeequ*(matrixLayout: cint; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; r: ptr cdouble; c: ptr cdouble;
            rowcnd: ptr cdouble; colcnd: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgeequ", dynlib: libName.}
proc zgeequWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; r: ptr cdouble;
                c: ptr cdouble; rowcnd: ptr cdouble; colcnd: ptr cdouble;
                amax: ptr cdouble): LapackInt {.importc: "LAPACKE_zgeequ_work",
    dynlib: libName.}
proc zgeequb*(matrixLayout: cint; m: LapackInt; n: LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; r: ptr cdouble; c: ptr cdouble;
             rowcnd: ptr cdouble; colcnd: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgeequb", dynlib: libName.}
proc zgeequbWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt; r: ptr cdouble;
                 c: ptr cdouble; rowcnd: ptr cdouble; colcnd: ptr cdouble;
                 amax: ptr cdouble): LapackInt {.importc: "LAPACKE_zgeequb_work",
    dynlib: libName.}
proc zgees*(matrixLayout: cint; jobvs: char; sort: char; select: Lapack_Z_Select1;
           n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
           sdim: ptr LapackInt; w: ptr LapackComplexDouble;
           vs: ptr LapackComplexDouble; ldvs: LapackInt): LapackInt {.
    importc: "LAPACKE_zgees", dynlib: libName.}
proc zgeesWork*(matrixLayout: cint; jobvs: char; sort: char; select: Lapack_Z_Select1;
               n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
               sdim: ptr LapackInt; w: ptr LapackComplexDouble;
               vs: ptr LapackComplexDouble; ldvs: LapackInt;
               work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble;
               bwork: ptr LapackLogical): LapackInt {.importc: "LAPACKE_zgees_work",
    dynlib: libName.}
proc zgeesx*(matrixLayout: cint; jobvs: char; sort: char; select: Lapack_Z_Select1;
            sense: char; n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            sdim: ptr LapackInt; w: ptr LapackComplexDouble;
            vs: ptr LapackComplexDouble; ldvs: LapackInt; rconde: ptr cdouble;
            rcondv: ptr cdouble): LapackInt {.importc: "LAPACKE_zgeesx",
    dynlib: libName.}
proc zgeesxWork*(matrixLayout: cint; jobvs: char; sort: char; select: Lapack_Z_Select1;
                sense: char; n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                sdim: ptr LapackInt; w: ptr LapackComplexDouble;
                vs: ptr LapackComplexDouble; ldvs: LapackInt; rconde: ptr cdouble;
                rcondv: ptr cdouble; work: ptr LapackComplexDouble; lwork: LapackInt;
                rwork: ptr cdouble; bwork: ptr LapackLogical): LapackInt {.
    importc: "LAPACKE_zgeesx_work", dynlib: libName.}
proc zgeev*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt;
           a: ptr LapackComplexDouble; lda: LapackInt; w: ptr LapackComplexDouble;
           vl: ptr LapackComplexDouble; ldvl: LapackInt; vr: ptr LapackComplexDouble;
           ldvr: LapackInt): LapackInt {.importc: "LAPACKE_zgeev", dynlib: libName.}
proc zgeevWork*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt;
               a: ptr LapackComplexDouble; lda: LapackInt;
               w: ptr LapackComplexDouble; vl: ptr LapackComplexDouble;
               ldvl: LapackInt; vr: ptr LapackComplexDouble; ldvr: LapackInt;
               work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgeev_work", dynlib: libName.}
proc zgeevx*(matrixLayout: cint; balanc: char; jobvl: char; jobvr: char; sense: char;
            n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            w: ptr LapackComplexDouble; vl: ptr LapackComplexDouble; ldvl: LapackInt;
            vr: ptr LapackComplexDouble; ldvr: LapackInt; ilo: ptr LapackInt;
            ihi: ptr LapackInt; scale: ptr cdouble; abnrm: ptr cdouble;
            rconde: ptr cdouble; rcondv: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgeevx", dynlib: libName.}
proc zgeevxWork*(matrixLayout: cint; balanc: char; jobvl: char; jobvr: char; sense: char;
                n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                w: ptr LapackComplexDouble; vl: ptr LapackComplexDouble;
                ldvl: LapackInt; vr: ptr LapackComplexDouble; ldvr: LapackInt;
                ilo: ptr LapackInt; ihi: ptr LapackInt; scale: ptr cdouble;
                abnrm: ptr cdouble; rconde: ptr cdouble; rcondv: ptr cdouble;
                work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgeevx_work", dynlib: libName.}
proc zgehrd*(matrixLayout: cint; n: LapackInt; ilo: LapackInt; ihi: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; tau: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zgehrd", dynlib: libName.}
proc zgehrdWork*(matrixLayout: cint; n: LapackInt; ilo: LapackInt; ihi: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                tau: ptr LapackComplexDouble; work: ptr LapackComplexDouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_zgehrd_work",
    dynlib: libName.}
proc zgejsv*(matrixLayout: cint; joba: char; jobu: char; jobv: char; jobr: char;
            jobt: char; jobp: char; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; sva: ptr cdouble;
            u: ptr LapackComplexDouble; ldu: LapackInt; v: ptr LapackComplexDouble;
            ldv: LapackInt; stat: ptr cdouble; istat: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zgejsv", dynlib: libName.}
proc zgejsvWork*(matrixLayout: cint; joba: char; jobu: char; jobv: char; jobr: char;
                jobt: char; jobp: char; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; sva: ptr cdouble;
                u: ptr LapackComplexDouble; ldu: LapackInt;
                v: ptr LapackComplexDouble; ldv: LapackInt;
                cwork: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble;
                lrwork: LapackInt; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zgejsv_work", dynlib: libName.}
proc zgelq2*(matrixLayout: cint; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; tau: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zgelq2", dynlib: libName.}
proc zgelq2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                tau: ptr LapackComplexDouble; work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zgelq2_work", dynlib: libName.}
proc zgelqf*(matrixLayout: cint; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; tau: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zgelqf", dynlib: libName.}
proc zgelqfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                tau: ptr LapackComplexDouble; work: ptr LapackComplexDouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_zgelqf_work",
    dynlib: libName.}
proc zgels*(matrixLayout: cint; trans: char; m: LapackInt; n: LapackInt; nrhs: LapackInt;
           a: ptr LapackComplexDouble; lda: LapackInt; b: ptr LapackComplexDouble;
           ldb: LapackInt): LapackInt {.importc: "LAPACKE_zgels", dynlib: libName.}
proc zgelsWork*(matrixLayout: cint; trans: char; m: LapackInt; n: LapackInt;
               nrhs: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
               b: ptr LapackComplexDouble; ldb: LapackInt;
               work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zgels_work", dynlib: libName.}
proc zgelsd*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; b: ptr LapackComplexDouble;
            ldb: LapackInt; s: ptr cdouble; rcond: cdouble; rank: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zgelsd", dynlib: libName.}
proc zgelsdWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt; s: ptr cdouble;
                rcond: cdouble; rank: ptr LapackInt; work: ptr LapackComplexDouble;
                lwork: LapackInt; rwork: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zgelsd_work", dynlib: libName.}
proc zgelss*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; b: ptr LapackComplexDouble;
            ldb: LapackInt; s: ptr cdouble; rcond: cdouble; rank: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zgelss", dynlib: libName.}
proc zgelssWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt; s: ptr cdouble;
                rcond: cdouble; rank: ptr LapackInt; work: ptr LapackComplexDouble;
                lwork: LapackInt; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgelss_work", dynlib: libName.}
proc zgelsy*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; b: ptr LapackComplexDouble;
            ldb: LapackInt; jpvt: ptr LapackInt; rcond: cdouble; rank: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zgelsy", dynlib: libName.}
proc zgelsyWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt; jpvt: ptr LapackInt;
                rcond: cdouble; rank: ptr LapackInt; work: ptr LapackComplexDouble;
                lwork: LapackInt; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgelsy_work", dynlib: libName.}
proc zgemqrt*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
             k: LapackInt; nb: LapackInt; v: ptr LapackComplexDouble; ldv: LapackInt;
             t: ptr LapackComplexDouble; ldt: LapackInt; c: ptr LapackComplexDouble;
             ldc: LapackInt): LapackInt {.importc: "LAPACKE_zgemqrt", dynlib: libName.}
proc zgemqrtWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt;
                 n: LapackInt; k: LapackInt; nb: LapackInt;
                 v: ptr LapackComplexDouble; ldv: LapackInt;
                 t: ptr LapackComplexDouble; ldt: LapackInt;
                 c: ptr LapackComplexDouble; ldc: LapackInt;
                 work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zgemqrt_work", dynlib: libName.}
proc zgeqlf*(matrixLayout: cint; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; tau: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zgeqlf", dynlib: libName.}
proc zgeqlfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                tau: ptr LapackComplexDouble; work: ptr LapackComplexDouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_zgeqlf_work",
    dynlib: libName.}
proc zgeqp3*(matrixLayout: cint; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; jpvt: ptr LapackInt;
            tau: ptr LapackComplexDouble): LapackInt {.importc: "LAPACKE_zgeqp3",
    dynlib: libName.}
proc zgeqp3Work*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; jpvt: ptr LapackInt;
                tau: ptr LapackComplexDouble; work: ptr LapackComplexDouble;
                lwork: LapackInt; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgeqp3_work", dynlib: libName.}
proc zgeqpf*(matrixLayout: cint; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; jpvt: ptr LapackInt;
            tau: ptr LapackComplexDouble): LapackInt {.importc: "LAPACKE_zgeqpf",
    dynlib: libName.}
proc zgeqpfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; jpvt: ptr LapackInt;
                tau: ptr LapackComplexDouble; work: ptr LapackComplexDouble;
                rwork: ptr cdouble): LapackInt {.importc: "LAPACKE_zgeqpf_work",
    dynlib: libName.}
proc zgeqr2*(matrixLayout: cint; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; tau: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zgeqr2", dynlib: libName.}
proc zgeqr2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                tau: ptr LapackComplexDouble; work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zgeqr2_work", dynlib: libName.}
proc zgeqrf*(matrixLayout: cint; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; tau: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zgeqrf", dynlib: libName.}
proc zgeqrfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                tau: ptr LapackComplexDouble; work: ptr LapackComplexDouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_zgeqrf_work",
    dynlib: libName.}
proc zgeqrfp*(matrixLayout: cint; m: LapackInt; n: LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; tau: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zgeqrfp", dynlib: libName.}
proc zgeqrfpWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt;
                 tau: ptr LapackComplexDouble; work: ptr LapackComplexDouble;
                 lwork: LapackInt): LapackInt {.importc: "LAPACKE_zgeqrfp_work",
    dynlib: libName.}
proc zgeqrt*(matrixLayout: cint; m: LapackInt; n: LapackInt; nb: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; t: ptr LapackComplexDouble;
            ldt: LapackInt): LapackInt {.importc: "LAPACKE_zgeqrt", dynlib: libName.}
proc zgeqrt2*(matrixLayout: cint; m: LapackInt; n: LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; t: ptr LapackComplexDouble;
             ldt: LapackInt): LapackInt {.importc: "LAPACKE_zgeqrt2", dynlib: libName.}
proc zgeqrt2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt;
                 t: ptr LapackComplexDouble; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_zgeqrt2_work", dynlib: libName.}
proc zgeqrt3*(matrixLayout: cint; m: LapackInt; n: LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; t: ptr LapackComplexDouble;
             ldt: LapackInt): LapackInt {.importc: "LAPACKE_zgeqrt3", dynlib: libName.}
proc zgeqrt3Work*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt;
                 t: ptr LapackComplexDouble; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_zgeqrt3_work", dynlib: libName.}
proc zgeqrtWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; nb: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                t: ptr LapackComplexDouble; ldt: LapackInt;
                work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zgeqrt_work", dynlib: libName.}
proc zgerfs*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; af: ptr LapackComplexDouble;
            ldaf: LapackInt; ipiv: ptr LapackInt; b: ptr LapackComplexDouble;
            ldb: LapackInt; x: ptr LapackComplexDouble; ldx: LapackInt;
            ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgerfs", dynlib: libName.}
proc zgerfsWork*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                af: ptr LapackComplexDouble; ldaf: LapackInt; ipiv: ptr LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; ferr: ptr cdouble;
                berr: ptr cdouble; work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgerfs_work", dynlib: libName.}
proc zgerfsx*(matrixLayout: cint; trans: char; equed: char; n: LapackInt;
             nrhs: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
             af: ptr LapackComplexDouble; ldaf: LapackInt; ipiv: ptr LapackInt;
             r: ptr cdouble; c: ptr cdouble; b: ptr LapackComplexDouble; ldb: LapackInt;
             x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
             berr: ptr cdouble; nErrBnds: LapackInt; errBndsNorm: ptr cdouble;
             errBndsComp: ptr cdouble; nparams: LapackInt; params: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgerfsx", dynlib: libName.}
proc zgerfsxWork*(matrixLayout: cint; trans: char; equed: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                 af: ptr LapackComplexDouble; ldaf: LapackInt; ipiv: ptr LapackInt;
                 r: ptr cdouble; c: ptr cdouble; b: ptr LapackComplexDouble;
                 ldb: LapackInt; x: ptr LapackComplexDouble; ldx: LapackInt;
                 rcond: ptr cdouble; berr: ptr cdouble; nErrBnds: LapackInt;
                 errBndsNorm: ptr cdouble; errBndsComp: ptr cdouble;
                 nparams: LapackInt; params: ptr cdouble;
                 work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgerfsx_work", dynlib: libName.}
proc zgerqf*(matrixLayout: cint; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; tau: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zgerqf", dynlib: libName.}
proc zgerqfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                tau: ptr LapackComplexDouble; work: ptr LapackComplexDouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_zgerqf_work",
    dynlib: libName.}
proc zgesdd*(matrixLayout: cint; jobz: char; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; s: ptr cdouble;
            u: ptr LapackComplexDouble; ldu: LapackInt; vt: ptr LapackComplexDouble;
            ldvt: LapackInt): LapackInt {.importc: "LAPACKE_zgesdd", dynlib: libName.}
proc zgesddWork*(matrixLayout: cint; jobz: char; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; s: ptr cdouble;
                u: ptr LapackComplexDouble; ldu: LapackInt;
                vt: ptr LapackComplexDouble; ldvt: LapackInt;
                work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_zgesdd_work",
    dynlib: libName.}
proc zgesv*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt;
           a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
           b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zgesv", dynlib: libName.}
proc zgesvWork*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt;
               a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
               b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zgesv_work", dynlib: libName.}
proc zgesvd*(matrixLayout: cint; jobu: char; jobvt: char; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; s: ptr cdouble;
            u: ptr LapackComplexDouble; ldu: LapackInt; vt: ptr LapackComplexDouble;
            ldvt: LapackInt; superb: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgesvd", dynlib: libName.}
proc zgesvdWork*(matrixLayout: cint; jobu: char; jobvt: char; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; s: ptr cdouble;
                u: ptr LapackComplexDouble; ldu: LapackInt;
                vt: ptr LapackComplexDouble; ldvt: LapackInt;
                work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgesvd_work", dynlib: libName.}
proc zgesvdx*(matrixLayout: cint; jobu: char; jobvt: char; range: char; m: LapackInt;
             n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt; vl: cdouble;
             vu: cdouble; il: LapackInt; iu: LapackInt; ns: ptr LapackInt;
             s: ptr cdouble; u: ptr LapackComplexDouble; ldu: LapackInt;
             vt: ptr LapackComplexDouble; ldvt: LapackInt; superb: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zgesvdx", dynlib: libName.}
proc zgesvdxWork*(matrixLayout: cint; jobu: char; jobvt: char; range: char; m: LapackInt;
                 n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt; vl: cdouble;
                 vu: cdouble; il: LapackInt; iu: LapackInt; ns: ptr LapackInt;
                 s: ptr cdouble; u: ptr LapackComplexDouble; ldu: LapackInt;
                 vt: ptr LapackComplexDouble; ldvt: LapackInt;
                 work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble;
                 iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_zgesvdx_work",
    dynlib: libName.}
proc zgesvj*(matrixLayout: cint; joba: char; jobu: char; jobv: char; m: LapackInt;
            n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt; sva: ptr cdouble;
            mv: LapackInt; v: ptr LapackComplexDouble; ldv: LapackInt;
            stat: ptr cdouble): LapackInt {.importc: "LAPACKE_zgesvj", dynlib: libName.}
proc zgesvjWork*(matrixLayout: cint; joba: char; jobu: char; jobv: char; m: LapackInt;
                n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                sva: ptr cdouble; mv: LapackInt; v: ptr LapackComplexDouble;
                ldv: LapackInt; cwork: ptr LapackComplexDouble; lwork: LapackInt;
                rwork: ptr cdouble; lrwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zgesvj_work", dynlib: libName.}
proc zgesvx*(matrixLayout: cint; fact: char; trans: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; af: ptr LapackComplexDouble;
            ldaf: LapackInt; ipiv: ptr LapackInt; equed: cstring; r: ptr cdouble;
            c: ptr cdouble; b: ptr LapackComplexDouble; ldb: LapackInt;
            x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
            ferr: ptr cdouble; berr: ptr cdouble; rpivot: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgesvx", dynlib: libName.}
proc zgesvxWork*(matrixLayout: cint; fact: char; trans: char; n: LapackInt;
                nrhs: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                af: ptr LapackComplexDouble; ldaf: LapackInt; ipiv: ptr LapackInt;
                equed: cstring; r: ptr cdouble; c: ptr cdouble;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
                ferr: ptr cdouble; berr: ptr cdouble; work: ptr LapackComplexDouble;
                rwork: ptr cdouble): LapackInt {.importc: "LAPACKE_zgesvx_work",
    dynlib: libName.}
proc zgesvxx*(matrixLayout: cint; fact: char; trans: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; af: ptr LapackComplexDouble;
             ldaf: LapackInt; ipiv: ptr LapackInt; equed: cstring; r: ptr cdouble;
             c: ptr cdouble; b: ptr LapackComplexDouble; ldb: LapackInt;
             x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
             rpvgrw: ptr cdouble; berr: ptr cdouble; nErrBnds: LapackInt;
             errBndsNorm: ptr cdouble; errBndsComp: ptr cdouble; nparams: LapackInt;
             params: ptr cdouble): LapackInt {.importc: "LAPACKE_zgesvxx",
    dynlib: libName.}
proc zgesvxxWork*(matrixLayout: cint; fact: char; trans: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                 af: ptr LapackComplexDouble; ldaf: LapackInt; ipiv: ptr LapackInt;
                 equed: cstring; r: ptr cdouble; c: ptr cdouble;
                 b: ptr LapackComplexDouble; ldb: LapackInt;
                 x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
                 rpvgrw: ptr cdouble; berr: ptr cdouble; nErrBnds: LapackInt;
                 errBndsNorm: ptr cdouble; errBndsComp: ptr cdouble;
                 nparams: LapackInt; params: ptr cdouble;
                 work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgesvxx_work", dynlib: libName.}
proc zgetf2*(matrixLayout: cint; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zgetf2", dynlib: libName.}
proc zgetf2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zgetf2_work", dynlib: libName.}
proc zgetrf*(matrixLayout: cint; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zgetrf", dynlib: libName.}
proc zgetrf2*(matrixLayout: cint; m: LapackInt; n: LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zgetrf2", dynlib: libName.}
proc zgetrf2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zgetrf2_work", dynlib: libName.}
proc zgetrfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zgetrf_work", dynlib: libName.}
proc zgetri*(matrixLayout: cint; n: LapackInt; a: ptr LapackComplexDouble;
            lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zgetri", dynlib: libName.}
proc zgetriWork*(matrixLayout: cint; n: LapackInt; a: ptr LapackComplexDouble;
                lda: LapackInt; ipiv: ptr LapackInt; work: ptr LapackComplexDouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_zgetri_work",
    dynlib: libName.}
proc zgetrs*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
            b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zgetrs", dynlib: libName.}
proc zgetrsWork*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zgetrs_work", dynlib: libName.}
proc zggbak*(matrixLayout: cint; job: char; side: char; n: LapackInt; ilo: LapackInt;
            ihi: LapackInt; lscale: ptr cdouble; rscale: ptr cdouble; m: LapackInt;
            v: ptr LapackComplexDouble; ldv: LapackInt): LapackInt {.
    importc: "LAPACKE_zggbak", dynlib: libName.}
proc zggbakWork*(matrixLayout: cint; job: char; side: char; n: LapackInt; ilo: LapackInt;
                ihi: LapackInt; lscale: ptr cdouble; rscale: ptr cdouble; m: LapackInt;
                v: ptr LapackComplexDouble; ldv: LapackInt): LapackInt {.
    importc: "LAPACKE_zggbak_work", dynlib: libName.}
proc zggbal*(matrixLayout: cint; job: char; n: LapackInt; a: ptr LapackComplexDouble;
            lda: LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt;
            ilo: ptr LapackInt; ihi: ptr LapackInt; lscale: ptr cdouble;
            rscale: ptr cdouble): LapackInt {.importc: "LAPACKE_zggbal",
    dynlib: libName.}
proc zggbalWork*(matrixLayout: cint; job: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt; ilo: ptr LapackInt;
                ihi: ptr LapackInt; lscale: ptr cdouble; rscale: ptr cdouble;
                work: ptr cdouble): LapackInt {.importc: "LAPACKE_zggbal_work",
    dynlib: libName.}
proc zgges*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
           selctg: Lapack_Z_Select2; n: LapackInt; a: ptr LapackComplexDouble;
           lda: LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt;
           sdim: ptr LapackInt; alpha: ptr LapackComplexDouble;
           beta: ptr LapackComplexDouble; vsl: ptr LapackComplexDouble;
           ldvsl: LapackInt; vsr: ptr LapackComplexDouble; ldvsr: LapackInt): LapackInt {.
    importc: "LAPACKE_zgges", dynlib: libName.}
proc zgges3*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
            selctg: Lapack_Z_Select2; n: LapackInt; a: ptr LapackComplexDouble;
            lda: LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt;
            sdim: ptr LapackInt; alpha: ptr LapackComplexDouble;
            beta: ptr LapackComplexDouble; vsl: ptr LapackComplexDouble;
            ldvsl: LapackInt; vsr: ptr LapackComplexDouble; ldvsr: LapackInt): LapackInt {.
    importc: "LAPACKE_zgges3", dynlib: libName.}
proc zgges3Work*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
                selctg: Lapack_Z_Select2; n: LapackInt; a: ptr LapackComplexDouble;
                lda: LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt;
                sdim: ptr LapackInt; alpha: ptr LapackComplexDouble;
                beta: ptr LapackComplexDouble; vsl: ptr LapackComplexDouble;
                ldvsl: LapackInt; vsr: ptr LapackComplexDouble; ldvsr: LapackInt;
                work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble;
                bwork: ptr LapackLogical): LapackInt {.
    importc: "LAPACKE_zgges3_work", dynlib: libName.}
proc zggesWork*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
               selctg: Lapack_Z_Select2; n: LapackInt; a: ptr LapackComplexDouble;
               lda: LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt;
               sdim: ptr LapackInt; alpha: ptr LapackComplexDouble;
               beta: ptr LapackComplexDouble; vsl: ptr LapackComplexDouble;
               ldvsl: LapackInt; vsr: ptr LapackComplexDouble; ldvsr: LapackInt;
               work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble;
               bwork: ptr LapackLogical): LapackInt {.importc: "LAPACKE_zgges_work",
    dynlib: libName.}
proc zggesx*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
            selctg: Lapack_Z_Select2; sense: char; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; b: ptr LapackComplexDouble;
            ldb: LapackInt; sdim: ptr LapackInt; alpha: ptr LapackComplexDouble;
            beta: ptr LapackComplexDouble; vsl: ptr LapackComplexDouble;
            ldvsl: LapackInt; vsr: ptr LapackComplexDouble; ldvsr: LapackInt;
            rconde: ptr cdouble; rcondv: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zggesx", dynlib: libName.}
proc zggesxWork*(matrixLayout: cint; jobvsl: char; jobvsr: char; sort: char;
                selctg: Lapack_Z_Select2; sense: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt; sdim: ptr LapackInt;
                alpha: ptr LapackComplexDouble; beta: ptr LapackComplexDouble;
                vsl: ptr LapackComplexDouble; ldvsl: LapackInt;
                vsr: ptr LapackComplexDouble; ldvsr: LapackInt; rconde: ptr cdouble;
                rcondv: ptr cdouble; work: ptr LapackComplexDouble; lwork: LapackInt;
                rwork: ptr cdouble; iwork: ptr LapackInt; liwork: LapackInt;
                bwork: ptr LapackLogical): LapackInt {.
    importc: "LAPACKE_zggesx_work", dynlib: libName.}
proc zggev*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt;
           a: ptr LapackComplexDouble; lda: LapackInt; b: ptr LapackComplexDouble;
           ldb: LapackInt; alpha: ptr LapackComplexDouble;
           beta: ptr LapackComplexDouble; vl: ptr LapackComplexDouble;
           ldvl: LapackInt; vr: ptr LapackComplexDouble; ldvr: LapackInt): LapackInt {.
    importc: "LAPACKE_zggev", dynlib: libName.}
proc zggev3*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; b: ptr LapackComplexDouble;
            ldb: LapackInt; alpha: ptr LapackComplexDouble;
            beta: ptr LapackComplexDouble; vl: ptr LapackComplexDouble;
            ldvl: LapackInt; vr: ptr LapackComplexDouble; ldvr: LapackInt): LapackInt {.
    importc: "LAPACKE_zggev3", dynlib: libName.}
proc zggev3Work*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                alpha: ptr LapackComplexDouble; beta: ptr LapackComplexDouble;
                vl: ptr LapackComplexDouble; ldvl: LapackInt;
                vr: ptr LapackComplexDouble; ldvr: LapackInt;
                work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zggev3_work", dynlib: libName.}
proc zggevWork*(matrixLayout: cint; jobvl: char; jobvr: char; n: LapackInt;
               a: ptr LapackComplexDouble; lda: LapackInt;
               b: ptr LapackComplexDouble; ldb: LapackInt;
               alpha: ptr LapackComplexDouble; beta: ptr LapackComplexDouble;
               vl: ptr LapackComplexDouble; ldvl: LapackInt;
               vr: ptr LapackComplexDouble; ldvr: LapackInt;
               work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zggev_work", dynlib: libName.}
proc zggevx*(matrixLayout: cint; balanc: char; jobvl: char; jobvr: char; sense: char;
            n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            b: ptr LapackComplexDouble; ldb: LapackInt;
            alpha: ptr LapackComplexDouble; beta: ptr LapackComplexDouble;
            vl: ptr LapackComplexDouble; ldvl: LapackInt;
            vr: ptr LapackComplexDouble; ldvr: LapackInt; ilo: ptr LapackInt;
            ihi: ptr LapackInt; lscale: ptr cdouble; rscale: ptr cdouble;
            abnrm: ptr cdouble; bbnrm: ptr cdouble; rconde: ptr cdouble;
            rcondv: ptr cdouble): LapackInt {.importc: "LAPACKE_zggevx",
    dynlib: libName.}
proc zggevxWork*(matrixLayout: cint; balanc: char; jobvl: char; jobvr: char; sense: char;
                n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                alpha: ptr LapackComplexDouble; beta: ptr LapackComplexDouble;
                vl: ptr LapackComplexDouble; ldvl: LapackInt;
                vr: ptr LapackComplexDouble; ldvr: LapackInt; ilo: ptr LapackInt;
                ihi: ptr LapackInt; lscale: ptr cdouble; rscale: ptr cdouble;
                abnrm: ptr cdouble; bbnrm: ptr cdouble; rconde: ptr cdouble;
                rcondv: ptr cdouble; work: ptr LapackComplexDouble; lwork: LapackInt;
                rwork: ptr cdouble; iwork: ptr LapackInt; bwork: ptr LapackLogical): LapackInt {.
    importc: "LAPACKE_zggevx_work", dynlib: libName.}
proc zggglm*(matrixLayout: cint; n: LapackInt; m: LapackInt; p: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; b: ptr LapackComplexDouble;
            ldb: LapackInt; d: ptr LapackComplexDouble; x: ptr LapackComplexDouble;
            y: ptr LapackComplexDouble): LapackInt {.importc: "LAPACKE_zggglm",
    dynlib: libName.}
proc zggglmWork*(matrixLayout: cint; n: LapackInt; m: LapackInt; p: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                d: ptr LapackComplexDouble; x: ptr LapackComplexDouble;
                y: ptr LapackComplexDouble; work: ptr LapackComplexDouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_zggglm_work",
    dynlib: libName.}
proc zgghd3*(matrixLayout: cint; compq: char; compz: char; n: LapackInt; ilo: LapackInt;
            ihi: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            b: ptr LapackComplexDouble; ldb: LapackInt; q: ptr LapackComplexDouble;
            ldq: LapackInt; z: ptr LapackComplexDouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_zgghd3", dynlib: libName.}
proc zgghd3Work*(matrixLayout: cint; compq: char; compz: char; n: LapackInt;
                ilo: LapackInt; ihi: LapackInt; a: ptr LapackComplexDouble;
                lda: LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt;
                q: ptr LapackComplexDouble; ldq: LapackInt;
                z: ptr LapackComplexDouble; ldz: LapackInt;
                work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zgghd3_work", dynlib: libName.}
proc zgghrd*(matrixLayout: cint; compq: char; compz: char; n: LapackInt; ilo: LapackInt;
            ihi: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            b: ptr LapackComplexDouble; ldb: LapackInt; q: ptr LapackComplexDouble;
            ldq: LapackInt; z: ptr LapackComplexDouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_zgghrd", dynlib: libName.}
proc zgghrdWork*(matrixLayout: cint; compq: char; compz: char; n: LapackInt;
                ilo: LapackInt; ihi: LapackInt; a: ptr LapackComplexDouble;
                lda: LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt;
                q: ptr LapackComplexDouble; ldq: LapackInt;
                z: ptr LapackComplexDouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_zgghrd_work", dynlib: libName.}
proc zgglse*(matrixLayout: cint; m: LapackInt; n: LapackInt; p: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; b: ptr LapackComplexDouble;
            ldb: LapackInt; c: ptr LapackComplexDouble; d: ptr LapackComplexDouble;
            x: ptr LapackComplexDouble): LapackInt {.importc: "LAPACKE_zgglse",
    dynlib: libName.}
proc zgglseWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; p: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                c: ptr LapackComplexDouble; d: ptr LapackComplexDouble;
                x: ptr LapackComplexDouble; work: ptr LapackComplexDouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_zgglse_work",
    dynlib: libName.}
proc zggqrf*(matrixLayout: cint; n: LapackInt; m: LapackInt; p: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt;
            taua: ptr LapackComplexDouble; b: ptr LapackComplexDouble; ldb: LapackInt;
            taub: ptr LapackComplexDouble): LapackInt {.importc: "LAPACKE_zggqrf",
    dynlib: libName.}
proc zggqrfWork*(matrixLayout: cint; n: LapackInt; m: LapackInt; p: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                taua: ptr LapackComplexDouble; b: ptr LapackComplexDouble;
                ldb: LapackInt; taub: ptr LapackComplexDouble;
                work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zggqrf_work", dynlib: libName.}
proc zggrqf*(matrixLayout: cint; m: LapackInt; p: LapackInt; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt;
            taua: ptr LapackComplexDouble; b: ptr LapackComplexDouble; ldb: LapackInt;
            taub: ptr LapackComplexDouble): LapackInt {.importc: "LAPACKE_zggrqf",
    dynlib: libName.}
proc zggrqfWork*(matrixLayout: cint; m: LapackInt; p: LapackInt; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                taua: ptr LapackComplexDouble; b: ptr LapackComplexDouble;
                ldb: LapackInt; taub: ptr LapackComplexDouble;
                work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zggrqf_work", dynlib: libName.}
proc zggsvd*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
            n: LapackInt; p: LapackInt; k: ptr LapackInt; l: ptr LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; b: ptr LapackComplexDouble;
            ldb: LapackInt; alpha: ptr cdouble; beta: ptr cdouble;
            u: ptr LapackComplexDouble; ldu: LapackInt; v: ptr LapackComplexDouble;
            ldv: LapackInt; q: ptr LapackComplexDouble; ldq: LapackInt;
            iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_zggsvd",
    dynlib: libName.}
proc zggsvd3*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
             n: LapackInt; p: LapackInt; k: ptr LapackInt; l: ptr LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; b: ptr LapackComplexDouble;
             ldb: LapackInt; alpha: ptr cdouble; beta: ptr cdouble;
             u: ptr LapackComplexDouble; ldu: LapackInt; v: ptr LapackComplexDouble;
             ldv: LapackInt; q: ptr LapackComplexDouble; ldq: LapackInt;
             iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_zggsvd3",
    dynlib: libName.}
proc zggsvd3Work*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
                 n: LapackInt; p: LapackInt; k: ptr LapackInt; l: ptr LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt;
                 b: ptr LapackComplexDouble; ldb: LapackInt; alpha: ptr cdouble;
                 beta: ptr cdouble; u: ptr LapackComplexDouble; ldu: LapackInt;
                 v: ptr LapackComplexDouble; ldv: LapackInt;
                 q: ptr LapackComplexDouble; ldq: LapackInt;
                 work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble;
                 iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_zggsvd3_work",
    dynlib: libName.}
proc zggsvdWork*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
                n: LapackInt; p: LapackInt; k: ptr LapackInt; l: ptr LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt; alpha: ptr cdouble;
                beta: ptr cdouble; u: ptr LapackComplexDouble; ldu: LapackInt;
                v: ptr LapackComplexDouble; ldv: LapackInt;
                q: ptr LapackComplexDouble; ldq: LapackInt;
                work: ptr LapackComplexDouble; rwork: ptr cdouble;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_zggsvd_work",
    dynlib: libName.}
proc zggsvp*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
            p: LapackInt; n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            b: ptr LapackComplexDouble; ldb: LapackInt; tola: cdouble; tolb: cdouble;
            k: ptr LapackInt; l: ptr LapackInt; u: ptr LapackComplexDouble;
            ldu: LapackInt; v: ptr LapackComplexDouble; ldv: LapackInt;
            q: ptr LapackComplexDouble; ldq: LapackInt): LapackInt {.
    importc: "LAPACKE_zggsvp", dynlib: libName.}
proc zggsvp3*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
             p: LapackInt; n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
             b: ptr LapackComplexDouble; ldb: LapackInt; tola: cdouble; tolb: cdouble;
             k: ptr LapackInt; l: ptr LapackInt; u: ptr LapackComplexDouble;
             ldu: LapackInt; v: ptr LapackComplexDouble; ldv: LapackInt;
             q: ptr LapackComplexDouble; ldq: LapackInt): LapackInt {.
    importc: "LAPACKE_zggsvp3", dynlib: libName.}
proc zggsvp3Work*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
                 p: LapackInt; n: LapackInt; a: ptr LapackComplexDouble;
                 lda: LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt;
                 tola: cdouble; tolb: cdouble; k: ptr LapackInt; l: ptr LapackInt;
                 u: ptr LapackComplexDouble; ldu: LapackInt;
                 v: ptr LapackComplexDouble; ldv: LapackInt;
                 q: ptr LapackComplexDouble; ldq: LapackInt; iwork: ptr LapackInt;
                 rwork: ptr cdouble; tau: ptr LapackComplexDouble;
                 work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zggsvp3_work", dynlib: libName.}
proc zggsvpWork*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
                p: LapackInt; n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt; tola: cdouble;
                tolb: cdouble; k: ptr LapackInt; l: ptr LapackInt;
                u: ptr LapackComplexDouble; ldu: LapackInt;
                v: ptr LapackComplexDouble; ldv: LapackInt;
                q: ptr LapackComplexDouble; ldq: LapackInt; iwork: ptr LapackInt;
                rwork: ptr cdouble; tau: ptr LapackComplexDouble;
                work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zggsvp_work", dynlib: libName.}
proc zgtcon*(norm: char; n: LapackInt; dl: ptr LapackComplexDouble;
            d: ptr LapackComplexDouble; du: ptr LapackComplexDouble;
            du2: ptr LapackComplexDouble; ipiv: ptr LapackInt; anorm: cdouble;
            rcond: ptr cdouble): LapackInt {.importc: "LAPACKE_zgtcon",
    dynlib: libName.}
proc zgtconWork*(norm: char; n: LapackInt; dl: ptr LapackComplexDouble;
                d: ptr LapackComplexDouble; du: ptr LapackComplexDouble;
                du2: ptr LapackComplexDouble; ipiv: ptr LapackInt; anorm: cdouble;
                rcond: ptr cdouble; work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zgtcon_work", dynlib: libName.}
proc zgtrfs*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
            dl: ptr LapackComplexDouble; d: ptr LapackComplexDouble;
            du: ptr LapackComplexDouble; dlf: ptr LapackComplexDouble;
            df: ptr LapackComplexDouble; duf: ptr LapackComplexDouble;
            du2: ptr LapackComplexDouble; ipiv: ptr LapackInt;
            b: ptr LapackComplexDouble; ldb: LapackInt; x: ptr LapackComplexDouble;
            ldx: LapackInt; ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgtrfs", dynlib: libName.}
proc zgtrfsWork*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
                dl: ptr LapackComplexDouble; d: ptr LapackComplexDouble;
                du: ptr LapackComplexDouble; dlf: ptr LapackComplexDouble;
                df: ptr LapackComplexDouble; duf: ptr LapackComplexDouble;
                du2: ptr LapackComplexDouble; ipiv: ptr LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; ferr: ptr cdouble;
                berr: ptr cdouble; work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgtrfs_work", dynlib: libName.}
proc zgtsv*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt;
           dl: ptr LapackComplexDouble; d: ptr LapackComplexDouble;
           du: ptr LapackComplexDouble; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zgtsv", dynlib: libName.}
proc zgtsvWork*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt;
               dl: ptr LapackComplexDouble; d: ptr LapackComplexDouble;
               du: ptr LapackComplexDouble; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zgtsv_work", dynlib: libName.}
proc zgtsvx*(matrixLayout: cint; fact: char; trans: char; n: LapackInt; nrhs: LapackInt;
            dl: ptr LapackComplexDouble; d: ptr LapackComplexDouble;
            du: ptr LapackComplexDouble; dlf: ptr LapackComplexDouble;
            df: ptr LapackComplexDouble; duf: ptr LapackComplexDouble;
            du2: ptr LapackComplexDouble; ipiv: ptr LapackInt;
            b: ptr LapackComplexDouble; ldb: LapackInt; x: ptr LapackComplexDouble;
            ldx: LapackInt; rcond: ptr cdouble; ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zgtsvx", dynlib: libName.}
proc zgtsvxWork*(matrixLayout: cint; fact: char; trans: char; n: LapackInt;
                nrhs: LapackInt; dl: ptr LapackComplexDouble;
                d: ptr LapackComplexDouble; du: ptr LapackComplexDouble;
                dlf: ptr LapackComplexDouble; df: ptr LapackComplexDouble;
                duf: ptr LapackComplexDouble; du2: ptr LapackComplexDouble;
                ipiv: ptr LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
                ferr: ptr cdouble; berr: ptr cdouble; work: ptr LapackComplexDouble;
                rwork: ptr cdouble): LapackInt {.importc: "LAPACKE_zgtsvx_work",
    dynlib: libName.}
proc zgttrf*(n: LapackInt; dl: ptr LapackComplexDouble; d: ptr LapackComplexDouble;
            du: ptr LapackComplexDouble; du2: ptr LapackComplexDouble;
            ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_zgttrf",
    dynlib: libName.}
proc zgttrfWork*(n: LapackInt; dl: ptr LapackComplexDouble;
                d: ptr LapackComplexDouble; du: ptr LapackComplexDouble;
                du2: ptr LapackComplexDouble; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zgttrf_work", dynlib: libName.}
proc zgttrs*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
            dl: ptr LapackComplexDouble; d: ptr LapackComplexDouble;
            du: ptr LapackComplexDouble; du2: ptr LapackComplexDouble;
            ipiv: ptr LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zgttrs", dynlib: libName.}
proc zgttrsWork*(matrixLayout: cint; trans: char; n: LapackInt; nrhs: LapackInt;
                dl: ptr LapackComplexDouble; d: ptr LapackComplexDouble;
                du: ptr LapackComplexDouble; du2: ptr LapackComplexDouble;
                ipiv: ptr LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zgttrs_work", dynlib: libName.}
proc zhbev*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; kd: LapackInt;
           ab: ptr LapackComplexDouble; ldab: LapackInt; w: ptr cdouble;
           z: ptr LapackComplexDouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_zhbev", dynlib: libName.}
proc zhbevWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; kd: LapackInt;
               ab: ptr LapackComplexDouble; ldab: LapackInt; w: ptr cdouble;
               z: ptr LapackComplexDouble; ldz: LapackInt;
               work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zhbev_work", dynlib: libName.}
proc zhbevd*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; kd: LapackInt;
            ab: ptr LapackComplexDouble; ldab: LapackInt; w: ptr cdouble;
            z: ptr LapackComplexDouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_zhbevd", dynlib: libName.}
proc zhbevdWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; kd: LapackInt;
                ab: ptr LapackComplexDouble; ldab: LapackInt; w: ptr cdouble;
                z: ptr LapackComplexDouble; ldz: LapackInt;
                work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble;
                lrwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zhbevd_work", dynlib: libName.}
proc zhbevx*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
            kd: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
            q: ptr LapackComplexDouble; ldq: LapackInt; vl: cdouble; vu: cdouble;
            il: LapackInt; iu: LapackInt; abstol: cdouble; m: ptr LapackInt;
            w: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt;
            ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_zhbevx",
    dynlib: libName.}
proc zhbevxWork*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                kd: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
                q: ptr LapackComplexDouble; ldq: LapackInt; vl: cdouble; vu: cdouble;
                il: LapackInt; iu: LapackInt; abstol: cdouble; m: ptr LapackInt;
                w: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt;
                work: ptr LapackComplexDouble; rwork: ptr cdouble;
                iwork: ptr LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zhbevx_work", dynlib: libName.}
proc zhbgst*(matrixLayout: cint; vect: char; uplo: char; n: LapackInt; ka: LapackInt;
            kb: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
            bb: ptr LapackComplexDouble; ldbb: LapackInt; x: ptr LapackComplexDouble;
            ldx: LapackInt): LapackInt {.importc: "LAPACKE_zhbgst", dynlib: libName.}
proc zhbgstWork*(matrixLayout: cint; vect: char; uplo: char; n: LapackInt; ka: LapackInt;
                kb: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
                bb: ptr LapackComplexDouble; ldbb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt;
                work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zhbgst_work", dynlib: libName.}
proc zhbgv*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; ka: LapackInt;
           kb: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
           bb: ptr LapackComplexDouble; ldbb: LapackInt; w: ptr cdouble;
           z: ptr LapackComplexDouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_zhbgv", dynlib: libName.}
proc zhbgvWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; ka: LapackInt;
               kb: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
               bb: ptr LapackComplexDouble; ldbb: LapackInt; w: ptr cdouble;
               z: ptr LapackComplexDouble; ldz: LapackInt;
               work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zhbgv_work", dynlib: libName.}
proc zhbgvd*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; ka: LapackInt;
            kb: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
            bb: ptr LapackComplexDouble; ldbb: LapackInt; w: ptr cdouble;
            z: ptr LapackComplexDouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_zhbgvd", dynlib: libName.}
proc zhbgvdWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt; ka: LapackInt;
                kb: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
                bb: ptr LapackComplexDouble; ldbb: LapackInt; w: ptr cdouble;
                z: ptr LapackComplexDouble; ldz: LapackInt;
                work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble;
                lrwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zhbgvd_work", dynlib: libName.}
proc zhbgvx*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
            ka: LapackInt; kb: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
            bb: ptr LapackComplexDouble; ldbb: LapackInt; q: ptr LapackComplexDouble;
            ldq: LapackInt; vl: cdouble; vu: cdouble; il: LapackInt; iu: LapackInt;
            abstol: cdouble; m: ptr LapackInt; w: ptr cdouble;
            z: ptr LapackComplexDouble; ldz: LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zhbgvx", dynlib: libName.}
proc zhbgvxWork*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                ka: LapackInt; kb: LapackInt; ab: ptr LapackComplexDouble;
                ldab: LapackInt; bb: ptr LapackComplexDouble; ldbb: LapackInt;
                q: ptr LapackComplexDouble; ldq: LapackInt; vl: cdouble; vu: cdouble;
                il: LapackInt; iu: LapackInt; abstol: cdouble; m: ptr LapackInt;
                w: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt;
                work: ptr LapackComplexDouble; rwork: ptr cdouble;
                iwork: ptr LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zhbgvx_work", dynlib: libName.}
proc zhbtrd*(matrixLayout: cint; vect: char; uplo: char; n: LapackInt; kd: LapackInt;
            ab: ptr LapackComplexDouble; ldab: LapackInt; d: ptr cdouble;
            e: ptr cdouble; q: ptr LapackComplexDouble; ldq: LapackInt): LapackInt {.
    importc: "LAPACKE_zhbtrd", dynlib: libName.}
proc zhbtrdWork*(matrixLayout: cint; vect: char; uplo: char; n: LapackInt; kd: LapackInt;
                ab: ptr LapackComplexDouble; ldab: LapackInt; d: ptr cdouble;
                e: ptr cdouble; q: ptr LapackComplexDouble; ldq: LapackInt;
                work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zhbtrd_work", dynlib: libName.}
proc zhecon*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexDouble;
            lda: LapackInt; ipiv: ptr LapackInt; anorm: cdouble; rcond: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zhecon", dynlib: libName.}
proc zheconWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                anorm: cdouble; rcond: ptr cdouble; work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zhecon_work", dynlib: libName.}
proc zheequb*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexDouble;
             lda: LapackInt; s: ptr cdouble; scond: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zheequb", dynlib: libName.}
proc zheequbWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt; s: ptr cdouble;
                 scond: ptr cdouble; amax: ptr cdouble; work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zheequb_work", dynlib: libName.}
proc zheev*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
           a: ptr LapackComplexDouble; lda: LapackInt; w: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zheev", dynlib: libName.}
proc zheevWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
               a: ptr LapackComplexDouble; lda: LapackInt; w: ptr cdouble;
               work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zheev_work", dynlib: libName.}
proc zheevd*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; w: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zheevd", dynlib: libName.}
proc zheevdWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; w: ptr cdouble;
                work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble;
                lrwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zheevd_work", dynlib: libName.}
proc zheevr*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; vl: cdouble; vu: cdouble;
            il: LapackInt; iu: LapackInt; abstol: cdouble; m: ptr LapackInt;
            w: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt;
            isuppz: ptr LapackInt): LapackInt {.importc: "LAPACKE_zheevr",
    dynlib: libName.}
proc zheevrWork*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; vl: cdouble; vu: cdouble;
                il: LapackInt; iu: LapackInt; abstol: cdouble; m: ptr LapackInt;
                w: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt;
                isuppz: ptr LapackInt; work: ptr LapackComplexDouble;
                lwork: LapackInt; rwork: ptr cdouble; lrwork: LapackInt;
                iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zheevr_work", dynlib: libName.}
proc zheevx*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; vl: cdouble; vu: cdouble;
            il: LapackInt; iu: LapackInt; abstol: cdouble; m: ptr LapackInt;
            w: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt;
            ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_zheevx",
    dynlib: libName.}
proc zheevxWork*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; vl: cdouble; vu: cdouble;
                il: LapackInt; iu: LapackInt; abstol: cdouble; m: ptr LapackInt;
                w: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt;
                work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble;
                iwork: ptr LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zheevx_work", dynlib: libName.}
proc zhegst*(matrixLayout: cint; itype: LapackInt; uplo: char; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; b: ptr LapackComplexDouble;
            ldb: LapackInt): LapackInt {.importc: "LAPACKE_zhegst", dynlib: libName.}
proc zhegstWork*(matrixLayout: cint; itype: LapackInt; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zhegst_work", dynlib: libName.}
proc zhegv*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char; n: LapackInt;
           a: ptr LapackComplexDouble; lda: LapackInt; b: ptr LapackComplexDouble;
           ldb: LapackInt; w: ptr cdouble): LapackInt {.importc: "LAPACKE_zhegv",
    dynlib: libName.}
proc zhegvWork*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
               n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
               b: ptr LapackComplexDouble; ldb: LapackInt; w: ptr cdouble;
               work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zhegv_work", dynlib: libName.}
proc zhegvd*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; b: ptr LapackComplexDouble;
            ldb: LapackInt; w: ptr cdouble): LapackInt {.importc: "LAPACKE_zhegvd",
    dynlib: libName.}
proc zhegvdWork*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
                n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt; w: ptr cdouble;
                work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble;
                lrwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zhegvd_work", dynlib: libName.}
proc zhegvx*(matrixLayout: cint; itype: LapackInt; jobz: char; range: char; uplo: char;
            n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            b: ptr LapackComplexDouble; ldb: LapackInt; vl: cdouble; vu: cdouble;
            il: LapackInt; iu: LapackInt; abstol: cdouble; m: ptr LapackInt;
            w: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt;
            ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_zhegvx",
    dynlib: libName.}
proc zhegvxWork*(matrixLayout: cint; itype: LapackInt; jobz: char; range: char;
                uplo: char; n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt; vl: cdouble; vu: cdouble;
                il: LapackInt; iu: LapackInt; abstol: cdouble; m: ptr LapackInt;
                w: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt;
                work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble;
                iwork: ptr LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zhegvx_work", dynlib: libName.}
proc zherfs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; af: ptr LapackComplexDouble;
            ldaf: LapackInt; ipiv: ptr LapackInt; b: ptr LapackComplexDouble;
            ldb: LapackInt; x: ptr LapackComplexDouble; ldx: LapackInt;
            ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zherfs", dynlib: libName.}
proc zherfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                af: ptr LapackComplexDouble; ldaf: LapackInt; ipiv: ptr LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; ferr: ptr cdouble;
                berr: ptr cdouble; work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zherfs_work", dynlib: libName.}
proc zherfsx*(matrixLayout: cint; uplo: char; equed: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; af: ptr LapackComplexDouble;
             ldaf: LapackInt; ipiv: ptr LapackInt; s: ptr cdouble;
             b: ptr LapackComplexDouble; ldb: LapackInt; x: ptr LapackComplexDouble;
             ldx: LapackInt; rcond: ptr cdouble; berr: ptr cdouble; nErrBnds: LapackInt;
             errBndsNorm: ptr cdouble; errBndsComp: ptr cdouble; nparams: LapackInt;
             params: ptr cdouble): LapackInt {.importc: "LAPACKE_zherfsx",
    dynlib: libName.}
proc zherfsxWork*(matrixLayout: cint; uplo: char; equed: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                 af: ptr LapackComplexDouble; ldaf: LapackInt; ipiv: ptr LapackInt;
                 s: ptr cdouble; b: ptr LapackComplexDouble; ldb: LapackInt;
                 x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
                 berr: ptr cdouble; nErrBnds: LapackInt; errBndsNorm: ptr cdouble;
                 errBndsComp: ptr cdouble; nparams: LapackInt; params: ptr cdouble;
                 work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zherfsx_work", dynlib: libName.}
proc zhesv*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
           a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
           b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zhesv", dynlib: libName.}
proc zhesvWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
               b: ptr LapackComplexDouble; ldb: LapackInt;
               work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zhesv_work", dynlib: libName.}
proc zhesvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; af: ptr LapackComplexDouble;
            ldaf: LapackInt; ipiv: ptr LapackInt; b: ptr LapackComplexDouble;
            ldb: LapackInt; x: ptr LapackComplexDouble; ldx: LapackInt;
            rcond: ptr cdouble; ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zhesvx", dynlib: libName.}
proc zhesvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                af: ptr LapackComplexDouble; ldaf: LapackInt; ipiv: ptr LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
                ferr: ptr cdouble; berr: ptr cdouble; work: ptr LapackComplexDouble;
                lwork: LapackInt; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zhesvx_work", dynlib: libName.}
proc zhesvxx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; af: ptr LapackComplexDouble;
             ldaf: LapackInt; ipiv: ptr LapackInt; equed: cstring; s: ptr cdouble;
             b: ptr LapackComplexDouble; ldb: LapackInt; x: ptr LapackComplexDouble;
             ldx: LapackInt; rcond: ptr cdouble; rpvgrw: ptr cdouble; berr: ptr cdouble;
             nErrBnds: LapackInt; errBndsNorm: ptr cdouble; errBndsComp: ptr cdouble;
             nparams: LapackInt; params: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zhesvxx", dynlib: libName.}
proc zhesvxxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                 af: ptr LapackComplexDouble; ldaf: LapackInt; ipiv: ptr LapackInt;
                 equed: cstring; s: ptr cdouble; b: ptr LapackComplexDouble;
                 ldb: LapackInt; x: ptr LapackComplexDouble; ldx: LapackInt;
                 rcond: ptr cdouble; rpvgrw: ptr cdouble; berr: ptr cdouble;
                 nErrBnds: LapackInt; errBndsNorm: ptr cdouble;
                 errBndsComp: ptr cdouble; nparams: LapackInt; params: ptr cdouble;
                 work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zhesvxx_work", dynlib: libName.}
proc zheswapr*(matrixLayout: cint; uplo: char; n: LapackInt;
              a: ptr LapackComplexDouble; lda: LapackInt; i1: LapackInt; i2: LapackInt): LapackInt {.
    importc: "LAPACKE_zheswapr", dynlib: libName.}
proc zheswaprWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                  a: ptr LapackComplexDouble; lda: LapackInt; i1: LapackInt;
                  i2: LapackInt): LapackInt {.importc: "LAPACKE_zheswapr_work",
    dynlib: libName.}
proc zhetrd*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexDouble;
            lda: LapackInt; d: ptr cdouble; e: ptr cdouble; tau: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zhetrd", dynlib: libName.}
proc zhetrdWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; d: ptr cdouble;
                e: ptr cdouble; tau: ptr LapackComplexDouble;
                work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zhetrd_work", dynlib: libName.}
proc zhetrf*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexDouble;
            lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zhetrf", dynlib: libName.}
proc zhetrfRook*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zhetrf_rook", dynlib: libName.}
proc zhetrfRookWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                    a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                    work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zhetrf_rook_work", dynlib: libName.}
proc zhetrfWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zhetrf_work", dynlib: libName.}
proc zhetri*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexDouble;
            lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zhetri", dynlib: libName.}
proc zhetri2*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexDouble;
             lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zhetri2", dynlib: libName.}
proc zhetri2Work*(matrixLayout: cint; uplo: char; n: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                 work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zhetri2_work", dynlib: libName.}
proc zhetri2x*(matrixLayout: cint; uplo: char; n: LapackInt;
              a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
              nb: LapackInt): LapackInt {.importc: "LAPACKE_zhetri2x",
                                       dynlib: libName.}
proc zhetri2xWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                  a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                  work: ptr LapackComplexDouble; nb: LapackInt): LapackInt {.
    importc: "LAPACKE_zhetri2x_work", dynlib: libName.}
proc zhetriWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zhetri_work", dynlib: libName.}
proc zhetrs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
            b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zhetrs", dynlib: libName.}
proc zhetrs2*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
             b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zhetrs2", dynlib: libName.}
proc zhetrs2Work*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                 b: ptr LapackComplexDouble; ldb: LapackInt;
                 work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zhetrs2_work", dynlib: libName.}
proc zhetrsRook*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zhetrs_rook", dynlib: libName.}
proc zhetrsRookWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                    a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                    b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zhetrs_rook_work", dynlib: libName.}
proc zhetrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zhetrs_work", dynlib: libName.}
proc zhfrk*(matrixLayout: cint; transr: char; uplo: char; trans: char; n: LapackInt;
           k: LapackInt; alpha: cdouble; a: ptr LapackComplexDouble; lda: LapackInt;
           beta: cdouble; c: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zhfrk", dynlib: libName.}
proc zhfrkWork*(matrixLayout: cint; transr: char; uplo: char; trans: char; n: LapackInt;
               k: LapackInt; alpha: cdouble; a: ptr LapackComplexDouble;
               lda: LapackInt; beta: cdouble; c: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zhfrk_work", dynlib: libName.}
proc zhgeqz*(matrixLayout: cint; job: char; compq: char; compz: char; n: LapackInt;
            ilo: LapackInt; ihi: LapackInt; h: ptr LapackComplexDouble; ldh: LapackInt;
            t: ptr LapackComplexDouble; ldt: LapackInt;
            alpha: ptr LapackComplexDouble; beta: ptr LapackComplexDouble;
            q: ptr LapackComplexDouble; ldq: LapackInt; z: ptr LapackComplexDouble;
            ldz: LapackInt): LapackInt {.importc: "LAPACKE_zhgeqz", dynlib: libName.}
proc zhgeqzWork*(matrixLayout: cint; job: char; compq: char; compz: char; n: LapackInt;
                ilo: LapackInt; ihi: LapackInt; h: ptr LapackComplexDouble;
                ldh: LapackInt; t: ptr LapackComplexDouble; ldt: LapackInt;
                alpha: ptr LapackComplexDouble; beta: ptr LapackComplexDouble;
                q: ptr LapackComplexDouble; ldq: LapackInt;
                z: ptr LapackComplexDouble; ldz: LapackInt;
                work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zhgeqz_work", dynlib: libName.}
proc zhpcon*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexDouble;
            ipiv: ptr LapackInt; anorm: cdouble; rcond: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zhpcon", dynlib: libName.}
proc zhpconWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexDouble; ipiv: ptr LapackInt; anorm: cdouble;
                rcond: ptr cdouble; work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zhpcon_work", dynlib: libName.}
proc zhpev*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
           ap: ptr LapackComplexDouble; w: ptr cdouble; z: ptr LapackComplexDouble;
           ldz: LapackInt): LapackInt {.importc: "LAPACKE_zhpev", dynlib: libName.}
proc zhpevWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
               ap: ptr LapackComplexDouble; w: ptr cdouble;
               z: ptr LapackComplexDouble; ldz: LapackInt;
               work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zhpev_work", dynlib: libName.}
proc zhpevd*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
            ap: ptr LapackComplexDouble; w: ptr cdouble; z: ptr LapackComplexDouble;
            ldz: LapackInt): LapackInt {.importc: "LAPACKE_zhpevd", dynlib: libName.}
proc zhpevdWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                ap: ptr LapackComplexDouble; w: ptr cdouble;
                z: ptr LapackComplexDouble; ldz: LapackInt;
                work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble;
                lrwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zhpevd_work", dynlib: libName.}
proc zhpevx*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
            ap: ptr LapackComplexDouble; vl: cdouble; vu: cdouble; il: LapackInt;
            iu: LapackInt; abstol: cdouble; m: ptr LapackInt; w: ptr cdouble;
            z: ptr LapackComplexDouble; ldz: LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zhpevx", dynlib: libName.}
proc zhpevxWork*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                ap: ptr LapackComplexDouble; vl: cdouble; vu: cdouble; il: LapackInt;
                iu: LapackInt; abstol: cdouble; m: ptr LapackInt; w: ptr cdouble;
                z: ptr LapackComplexDouble; ldz: LapackInt;
                work: ptr LapackComplexDouble; rwork: ptr cdouble;
                iwork: ptr LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zhpevx_work", dynlib: libName.}
proc zhpgst*(matrixLayout: cint; itype: LapackInt; uplo: char; n: LapackInt;
            ap: ptr LapackComplexDouble; bp: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zhpgst", dynlib: libName.}
proc zhpgstWork*(matrixLayout: cint; itype: LapackInt; uplo: char; n: LapackInt;
                ap: ptr LapackComplexDouble; bp: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zhpgst_work", dynlib: libName.}
proc zhpgv*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char; n: LapackInt;
           ap: ptr LapackComplexDouble; bp: ptr LapackComplexDouble; w: ptr cdouble;
           z: ptr LapackComplexDouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_zhpgv", dynlib: libName.}
proc zhpgvWork*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
               n: LapackInt; ap: ptr LapackComplexDouble;
               bp: ptr LapackComplexDouble; w: ptr cdouble;
               z: ptr LapackComplexDouble; ldz: LapackInt;
               work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zhpgv_work", dynlib: libName.}
proc zhpgvd*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char; n: LapackInt;
            ap: ptr LapackComplexDouble; bp: ptr LapackComplexDouble; w: ptr cdouble;
            z: ptr LapackComplexDouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_zhpgvd", dynlib: libName.}
proc zhpgvdWork*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
                n: LapackInt; ap: ptr LapackComplexDouble;
                bp: ptr LapackComplexDouble; w: ptr cdouble;
                z: ptr LapackComplexDouble; ldz: LapackInt;
                work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble;
                lrwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zhpgvd_work", dynlib: libName.}
proc zhpgvx*(matrixLayout: cint; itype: LapackInt; jobz: char; range: char; uplo: char;
            n: LapackInt; ap: ptr LapackComplexDouble; bp: ptr LapackComplexDouble;
            vl: cdouble; vu: cdouble; il: LapackInt; iu: LapackInt; abstol: cdouble;
            m: ptr LapackInt; w: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt;
            ifail: ptr LapackInt): LapackInt {.importc: "LAPACKE_zhpgvx",
    dynlib: libName.}
proc zhpgvxWork*(matrixLayout: cint; itype: LapackInt; jobz: char; range: char;
                uplo: char; n: LapackInt; ap: ptr LapackComplexDouble;
                bp: ptr LapackComplexDouble; vl: cdouble; vu: cdouble; il: LapackInt;
                iu: LapackInt; abstol: cdouble; m: ptr LapackInt; w: ptr cdouble;
                z: ptr LapackComplexDouble; ldz: LapackInt;
                work: ptr LapackComplexDouble; rwork: ptr cdouble;
                iwork: ptr LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zhpgvx_work", dynlib: libName.}
proc zhprfs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr LapackComplexDouble; afp: ptr LapackComplexDouble;
            ipiv: ptr LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt;
            x: ptr LapackComplexDouble; ldx: LapackInt; ferr: ptr cdouble;
            berr: ptr cdouble): LapackInt {.importc: "LAPACKE_zhprfs", dynlib: libName.}
proc zhprfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                ap: ptr LapackComplexDouble; afp: ptr LapackComplexDouble;
                ipiv: ptr LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; ferr: ptr cdouble;
                berr: ptr cdouble; work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zhprfs_work", dynlib: libName.}
proc zhpsv*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
           ap: ptr LapackComplexDouble; ipiv: ptr LapackInt;
           b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zhpsv", dynlib: libName.}
proc zhpsvWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               ap: ptr LapackComplexDouble; ipiv: ptr LapackInt;
               b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zhpsv_work", dynlib: libName.}
proc zhpsvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr LapackComplexDouble; afp: ptr LapackComplexDouble;
            ipiv: ptr LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt;
            x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
            ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zhpsvx", dynlib: libName.}
proc zhpsvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; ap: ptr LapackComplexDouble;
                afp: ptr LapackComplexDouble; ipiv: ptr LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
                ferr: ptr cdouble; berr: ptr cdouble; work: ptr LapackComplexDouble;
                rwork: ptr cdouble): LapackInt {.importc: "LAPACKE_zhpsvx_work",
    dynlib: libName.}
proc zhptrd*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexDouble;
            d: ptr cdouble; e: ptr cdouble; tau: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zhptrd", dynlib: libName.}
proc zhptrdWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexDouble; d: ptr cdouble; e: ptr cdouble;
                tau: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zhptrd_work", dynlib: libName.}
proc zhptrf*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexDouble;
            ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_zhptrf",
    dynlib: libName.}
proc zhptrfWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexDouble; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zhptrf_work", dynlib: libName.}
proc zhptri*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexDouble;
            ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_zhptri",
    dynlib: libName.}
proc zhptriWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexDouble; ipiv: ptr LapackInt;
                work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zhptri_work", dynlib: libName.}
proc zhptrs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr LapackComplexDouble; ipiv: ptr LapackInt;
            b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zhptrs", dynlib: libName.}
proc zhptrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                ap: ptr LapackComplexDouble; ipiv: ptr LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zhptrs_work", dynlib: libName.}
proc zhsein*(matrixLayout: cint; job: char; eigsrc: char; initv: char;
            select: ptr LapackLogical; n: LapackInt; h: ptr LapackComplexDouble;
            ldh: LapackInt; w: ptr LapackComplexDouble; vl: ptr LapackComplexDouble;
            ldvl: LapackInt; vr: ptr LapackComplexDouble; ldvr: LapackInt;
            mm: LapackInt; m: ptr LapackInt; ifaill: ptr LapackInt;
            ifailr: ptr LapackInt): LapackInt {.importc: "LAPACKE_zhsein",
    dynlib: libName.}
proc zhseinWork*(matrixLayout: cint; job: char; eigsrc: char; initv: char;
                select: ptr LapackLogical; n: LapackInt; h: ptr LapackComplexDouble;
                ldh: LapackInt; w: ptr LapackComplexDouble;
                vl: ptr LapackComplexDouble; ldvl: LapackInt;
                vr: ptr LapackComplexDouble; ldvr: LapackInt; mm: LapackInt;
                m: ptr LapackInt; work: ptr LapackComplexDouble; rwork: ptr cdouble;
                ifaill: ptr LapackInt; ifailr: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zhsein_work", dynlib: libName.}
proc zhseqr*(matrixLayout: cint; job: char; compz: char; n: LapackInt; ilo: LapackInt;
            ihi: LapackInt; h: ptr LapackComplexDouble; ldh: LapackInt;
            w: ptr LapackComplexDouble; z: ptr LapackComplexDouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_zhseqr", dynlib: libName.}
proc zhseqrWork*(matrixLayout: cint; job: char; compz: char; n: LapackInt;
                ilo: LapackInt; ihi: LapackInt; h: ptr LapackComplexDouble;
                ldh: LapackInt; w: ptr LapackComplexDouble;
                z: ptr LapackComplexDouble; ldz: LapackInt;
                work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zhseqr_work", dynlib: libName.}
proc zlacgv*(n: LapackInt; x: ptr LapackComplexDouble; incx: LapackInt): LapackInt {.
    importc: "LAPACKE_zlacgv", dynlib: libName.}
proc zlacgvWork*(n: LapackInt; x: ptr LapackComplexDouble; incx: LapackInt): LapackInt {.
    importc: "LAPACKE_zlacgv_work", dynlib: libName.}
proc zlacn2*(n: LapackInt; v: ptr LapackComplexDouble; x: ptr LapackComplexDouble;
            est: ptr cdouble; kase: ptr LapackInt; isave: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zlacn2", dynlib: libName.}
proc zlacn2Work*(n: LapackInt; v: ptr LapackComplexDouble; x: ptr LapackComplexDouble;
                est: ptr cdouble; kase: ptr LapackInt; isave: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zlacn2_work", dynlib: libName.}
proc zlacp2*(matrixLayout: cint; uplo: char; m: LapackInt; n: LapackInt; a: ptr cdouble;
            lda: LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zlacp2", dynlib: libName.}
proc zlacp2Work*(matrixLayout: cint; uplo: char; m: LapackInt; n: LapackInt;
                a: ptr cdouble; lda: LapackInt; b: ptr LapackComplexDouble;
                ldb: LapackInt): LapackInt {.importc: "LAPACKE_zlacp2_work",
    dynlib: libName.}
proc zlacpy*(matrixLayout: cint; uplo: char; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; b: ptr LapackComplexDouble;
            ldb: LapackInt): LapackInt {.importc: "LAPACKE_zlacpy", dynlib: libName.}
proc zlacpyWork*(matrixLayout: cint; uplo: char; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zlacpy_work", dynlib: libName.}
proc zlag2c*(matrixLayout: cint; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; sa: ptr LapackComplexFloat;
            ldsa: LapackInt): LapackInt {.importc: "LAPACKE_zlag2c", dynlib: libName.}
proc zlag2cWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                sa: ptr LapackComplexFloat; ldsa: LapackInt): LapackInt {.
    importc: "LAPACKE_zlag2c_work", dynlib: libName.}
proc zlagge*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
            ku: LapackInt; d: ptr cdouble; a: ptr LapackComplexDouble; lda: LapackInt;
            iseed: ptr LapackInt): LapackInt {.importc: "LAPACKE_zlagge",
    dynlib: libName.}
proc zlaggeWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; kl: LapackInt;
                ku: LapackInt; d: ptr cdouble; a: ptr LapackComplexDouble;
                lda: LapackInt; iseed: ptr LapackInt; work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zlagge_work", dynlib: libName.}
proc zlaghe*(matrixLayout: cint; n: LapackInt; k: LapackInt; d: ptr cdouble;
            a: ptr LapackComplexDouble; lda: LapackInt; iseed: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zlaghe", dynlib: libName.}
proc zlagheWork*(matrixLayout: cint; n: LapackInt; k: LapackInt; d: ptr cdouble;
                a: ptr LapackComplexDouble; lda: LapackInt; iseed: ptr LapackInt;
                work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zlaghe_work", dynlib: libName.}
proc zlagsy*(matrixLayout: cint; n: LapackInt; k: LapackInt; d: ptr cdouble;
            a: ptr LapackComplexDouble; lda: LapackInt; iseed: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zlagsy", dynlib: libName.}
proc zlagsyWork*(matrixLayout: cint; n: LapackInt; k: LapackInt; d: ptr cdouble;
                a: ptr LapackComplexDouble; lda: LapackInt; iseed: ptr LapackInt;
                work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zlagsy_work", dynlib: libName.}
proc zlapmr*(matrixLayout: cint; forwrd: LapackLogical; m: LapackInt; n: LapackInt;
            x: ptr LapackComplexDouble; ldx: LapackInt; k: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zlapmr", dynlib: libName.}
proc zlapmrWork*(matrixLayout: cint; forwrd: LapackLogical; m: LapackInt; n: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; k: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zlapmr_work", dynlib: libName.}
proc zlapmt*(matrixLayout: cint; forwrd: LapackLogical; m: LapackInt; n: LapackInt;
            x: ptr LapackComplexDouble; ldx: LapackInt; k: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zlapmt", dynlib: libName.}
proc zlapmtWork*(matrixLayout: cint; forwrd: LapackLogical; m: LapackInt; n: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; k: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zlapmt_work", dynlib: libName.}
proc zlarfb*(matrixLayout: cint; side: char; trans: char; direct: char; storev: char;
            m: LapackInt; n: LapackInt; k: LapackInt; v: ptr LapackComplexDouble;
            ldv: LapackInt; t: ptr LapackComplexDouble; ldt: LapackInt;
            c: ptr LapackComplexDouble; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_zlarfb", dynlib: libName.}
proc zlarfbWork*(matrixLayout: cint; side: char; trans: char; direct: char; storev: char;
                m: LapackInt; n: LapackInt; k: LapackInt; v: ptr LapackComplexDouble;
                ldv: LapackInt; t: ptr LapackComplexDouble; ldt: LapackInt;
                c: ptr LapackComplexDouble; ldc: LapackInt;
                work: ptr LapackComplexDouble; ldwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zlarfb_work", dynlib: libName.}
proc zlarfg*(n: LapackInt; alpha: ptr LapackComplexDouble; x: ptr LapackComplexDouble;
            incx: LapackInt; tau: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zlarfg", dynlib: libName.}
proc zlarfgWork*(n: LapackInt; alpha: ptr LapackComplexDouble;
                x: ptr LapackComplexDouble; incx: LapackInt;
                tau: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zlarfg_work", dynlib: libName.}
proc zlarft*(matrixLayout: cint; direct: char; storev: char; n: LapackInt; k: LapackInt;
            v: ptr LapackComplexDouble; ldv: LapackInt; tau: ptr LapackComplexDouble;
            t: ptr LapackComplexDouble; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_zlarft", dynlib: libName.}
proc zlarftWork*(matrixLayout: cint; direct: char; storev: char; n: LapackInt;
                k: LapackInt; v: ptr LapackComplexDouble; ldv: LapackInt;
                tau: ptr LapackComplexDouble; t: ptr LapackComplexDouble;
                ldt: LapackInt): LapackInt {.importc: "LAPACKE_zlarft_work",
    dynlib: libName.}
proc zlarfx*(matrixLayout: cint; side: char; m: LapackInt; n: LapackInt;
            v: ptr LapackComplexDouble; tau: LapackComplexDouble;
            c: ptr LapackComplexDouble; ldc: LapackInt; work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zlarfx", dynlib: libName.}
proc zlarfxWork*(matrixLayout: cint; side: char; m: LapackInt; n: LapackInt;
                v: ptr LapackComplexDouble; tau: LapackComplexDouble;
                c: ptr LapackComplexDouble; ldc: LapackInt;
                work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zlarfx_work", dynlib: libName.}
proc zlarnv*(idist: LapackInt; iseed: ptr LapackInt; n: LapackInt;
            x: ptr LapackComplexDouble): LapackInt {.importc: "LAPACKE_zlarnv",
    dynlib: libName.}
proc zlarnvWork*(idist: LapackInt; iseed: ptr LapackInt; n: LapackInt;
                x: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zlarnv_work", dynlib: libName.}
proc zlascl*(matrixLayout: cint; `type`: char; kl: LapackInt; ku: LapackInt;
            cfrom: cdouble; cto: cdouble; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_zlascl", dynlib: libName.}
proc zlasclWork*(matrixLayout: cint; `type`: char; kl: LapackInt; ku: LapackInt;
                cfrom: cdouble; cto: cdouble; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_zlascl_work", dynlib: libName.}
proc zlaset*(matrixLayout: cint; uplo: char; m: LapackInt; n: LapackInt;
            alpha: LapackComplexDouble; beta: LapackComplexDouble;
            a: ptr LapackComplexDouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_zlaset", dynlib: libName.}
proc zlasetWork*(matrixLayout: cint; uplo: char; m: LapackInt; n: LapackInt;
                alpha: LapackComplexDouble; beta: LapackComplexDouble;
                a: ptr LapackComplexDouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_zlaset_work", dynlib: libName.}
proc zlaswp*(matrixLayout: cint; n: LapackInt; a: ptr LapackComplexDouble;
            lda: LapackInt; k1: LapackInt; k2: LapackInt; ipiv: ptr LapackInt;
            incx: LapackInt): LapackInt {.importc: "LAPACKE_zlaswp", dynlib: libName.}
proc zlaswpWork*(matrixLayout: cint; n: LapackInt; a: ptr LapackComplexDouble;
                lda: LapackInt; k1: LapackInt; k2: LapackInt; ipiv: ptr LapackInt;
                incx: LapackInt): LapackInt {.importc: "LAPACKE_zlaswp_work",
    dynlib: libName.}
proc zlatms*(matrixLayout: cint; m: LapackInt; n: LapackInt; dist: char;
            iseed: ptr LapackInt; sym: char; d: ptr cdouble; mode: LapackInt;
            cond: cdouble; dmax: cdouble; kl: LapackInt; ku: LapackInt; pack: char;
            a: ptr LapackComplexDouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_zlatms", dynlib: libName.}
proc zlatmsWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; dist: char;
                iseed: ptr LapackInt; sym: char; d: ptr cdouble; mode: LapackInt;
                cond: cdouble; dmax: cdouble; kl: LapackInt; ku: LapackInt; pack: char;
                a: ptr LapackComplexDouble; lda: LapackInt;
                work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zlatms_work", dynlib: libName.}
proc zlauum*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexDouble;
            lda: LapackInt): LapackInt {.importc: "LAPACKE_zlauum", dynlib: libName.}
proc zlauumWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_zlauum_work", dynlib: libName.}
proc zpbcon*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
            ab: ptr LapackComplexDouble; ldab: LapackInt; anorm: cdouble;
            rcond: ptr cdouble): LapackInt {.importc: "LAPACKE_zpbcon",
    dynlib: libName.}
proc zpbconWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
                ab: ptr LapackComplexDouble; ldab: LapackInt; anorm: cdouble;
                rcond: ptr cdouble; work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zpbcon_work", dynlib: libName.}
proc zpbequ*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
            ab: ptr LapackComplexDouble; ldab: LapackInt; s: ptr cdouble;
            scond: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zpbequ", dynlib: libName.}
proc zpbequWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
                ab: ptr LapackComplexDouble; ldab: LapackInt; s: ptr cdouble;
                scond: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zpbequ_work", dynlib: libName.}
proc zpbrfs*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
            nrhs: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
            afb: ptr LapackComplexDouble; ldafb: LapackInt;
            b: ptr LapackComplexDouble; ldb: LapackInt; x: ptr LapackComplexDouble;
            ldx: LapackInt; ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zpbrfs", dynlib: libName.}
proc zpbrfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
                nrhs: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
                afb: ptr LapackComplexDouble; ldafb: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; ferr: ptr cdouble;
                berr: ptr cdouble; work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zpbrfs_work", dynlib: libName.}
proc zpbstf*(matrixLayout: cint; uplo: char; n: LapackInt; kb: LapackInt;
            bb: ptr LapackComplexDouble; ldbb: LapackInt): LapackInt {.
    importc: "LAPACKE_zpbstf", dynlib: libName.}
proc zpbstfWork*(matrixLayout: cint; uplo: char; n: LapackInt; kb: LapackInt;
                bb: ptr LapackComplexDouble; ldbb: LapackInt): LapackInt {.
    importc: "LAPACKE_zpbstf_work", dynlib: libName.}
proc zpbsv*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt; nrhs: LapackInt;
           ab: ptr LapackComplexDouble; ldab: LapackInt; b: ptr LapackComplexDouble;
           ldb: LapackInt): LapackInt {.importc: "LAPACKE_zpbsv", dynlib: libName.}
proc zpbsvWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
               nrhs: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
               b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zpbsv_work", dynlib: libName.}
proc zpbsvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; kd: LapackInt;
            nrhs: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
            afb: ptr LapackComplexDouble; ldafb: LapackInt; equed: cstring;
            s: ptr cdouble; b: ptr LapackComplexDouble; ldb: LapackInt;
            x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
            ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zpbsvx", dynlib: libName.}
proc zpbsvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; kd: LapackInt;
                nrhs: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
                afb: ptr LapackComplexDouble; ldafb: LapackInt; equed: cstring;
                s: ptr cdouble; b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
                ferr: ptr cdouble; berr: ptr cdouble; work: ptr LapackComplexDouble;
                rwork: ptr cdouble): LapackInt {.importc: "LAPACKE_zpbsvx_work",
    dynlib: libName.}
proc zpbtrf*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
            ab: ptr LapackComplexDouble; ldab: LapackInt): LapackInt {.
    importc: "LAPACKE_zpbtrf", dynlib: libName.}
proc zpbtrfWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
                ab: ptr LapackComplexDouble; ldab: LapackInt): LapackInt {.
    importc: "LAPACKE_zpbtrf_work", dynlib: libName.}
proc zpbtrs*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
            nrhs: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
            b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zpbtrs", dynlib: libName.}
proc zpbtrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; kd: LapackInt;
                nrhs: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zpbtrs_work", dynlib: libName.}
proc zpftrf*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
            a: ptr LapackComplexDouble): LapackInt {.importc: "LAPACKE_zpftrf",
    dynlib: libName.}
proc zpftrfWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zpftrf_work", dynlib: libName.}
proc zpftri*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
            a: ptr LapackComplexDouble): LapackInt {.importc: "LAPACKE_zpftri",
    dynlib: libName.}
proc zpftriWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zpftri_work", dynlib: libName.}
proc zpftrs*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexDouble; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zpftrs", dynlib: libName.}
proc zpftrsWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; a: ptr LapackComplexDouble;
                b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zpftrs_work", dynlib: libName.}
proc zpocon*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexDouble;
            lda: LapackInt; anorm: cdouble; rcond: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zpocon", dynlib: libName.}
proc zpoconWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; anorm: cdouble;
                rcond: ptr cdouble; work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zpocon_work", dynlib: libName.}
proc zpoequ*(matrixLayout: cint; n: LapackInt; a: ptr LapackComplexDouble;
            lda: LapackInt; s: ptr cdouble; scond: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zpoequ", dynlib: libName.}
proc zpoequWork*(matrixLayout: cint; n: LapackInt; a: ptr LapackComplexDouble;
                lda: LapackInt; s: ptr cdouble; scond: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zpoequ_work", dynlib: libName.}
proc zpoequb*(matrixLayout: cint; n: LapackInt; a: ptr LapackComplexDouble;
             lda: LapackInt; s: ptr cdouble; scond: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zpoequb", dynlib: libName.}
proc zpoequbWork*(matrixLayout: cint; n: LapackInt; a: ptr LapackComplexDouble;
                 lda: LapackInt; s: ptr cdouble; scond: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zpoequb_work", dynlib: libName.}
proc zporfs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; af: ptr LapackComplexDouble;
            ldaf: LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt;
            x: ptr LapackComplexDouble; ldx: LapackInt; ferr: ptr cdouble;
            berr: ptr cdouble): LapackInt {.importc: "LAPACKE_zporfs", dynlib: libName.}
proc zporfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                af: ptr LapackComplexDouble; ldaf: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; ferr: ptr cdouble;
                berr: ptr cdouble; work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zporfs_work", dynlib: libName.}
proc zporfsx*(matrixLayout: cint; uplo: char; equed: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; af: ptr LapackComplexDouble;
             ldaf: LapackInt; s: ptr cdouble; b: ptr LapackComplexDouble;
             ldb: LapackInt; x: ptr LapackComplexDouble; ldx: LapackInt;
             rcond: ptr cdouble; berr: ptr cdouble; nErrBnds: LapackInt;
             errBndsNorm: ptr cdouble; errBndsComp: ptr cdouble; nparams: LapackInt;
             params: ptr cdouble): LapackInt {.importc: "LAPACKE_zporfsx",
    dynlib: libName.}
proc zporfsxWork*(matrixLayout: cint; uplo: char; equed: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                 af: ptr LapackComplexDouble; ldaf: LapackInt; s: ptr cdouble;
                 b: ptr LapackComplexDouble; ldb: LapackInt;
                 x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
                 berr: ptr cdouble; nErrBnds: LapackInt; errBndsNorm: ptr cdouble;
                 errBndsComp: ptr cdouble; nparams: LapackInt; params: ptr cdouble;
                 work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zporfsx_work", dynlib: libName.}
proc zposv*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
           a: ptr LapackComplexDouble; lda: LapackInt; b: ptr LapackComplexDouble;
           ldb: LapackInt): LapackInt {.importc: "LAPACKE_zposv", dynlib: libName.}
proc zposvWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               a: ptr LapackComplexDouble; lda: LapackInt;
               b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zposv_work", dynlib: libName.}
proc zposvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; af: ptr LapackComplexDouble;
            ldaf: LapackInt; equed: cstring; s: ptr cdouble;
            b: ptr LapackComplexDouble; ldb: LapackInt; x: ptr LapackComplexDouble;
            ldx: LapackInt; rcond: ptr cdouble; ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zposvx", dynlib: libName.}
proc zposvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                af: ptr LapackComplexDouble; ldaf: LapackInt; equed: cstring;
                s: ptr cdouble; b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
                ferr: ptr cdouble; berr: ptr cdouble; work: ptr LapackComplexDouble;
                rwork: ptr cdouble): LapackInt {.importc: "LAPACKE_zposvx_work",
    dynlib: libName.}
proc zposvxx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; af: ptr LapackComplexDouble;
             ldaf: LapackInt; equed: cstring; s: ptr cdouble;
             b: ptr LapackComplexDouble; ldb: LapackInt; x: ptr LapackComplexDouble;
             ldx: LapackInt; rcond: ptr cdouble; rpvgrw: ptr cdouble; berr: ptr cdouble;
             nErrBnds: LapackInt; errBndsNorm: ptr cdouble; errBndsComp: ptr cdouble;
             nparams: LapackInt; params: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zposvxx", dynlib: libName.}
proc zposvxxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                 af: ptr LapackComplexDouble; ldaf: LapackInt; equed: cstring;
                 s: ptr cdouble; b: ptr LapackComplexDouble; ldb: LapackInt;
                 x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
                 rpvgrw: ptr cdouble; berr: ptr cdouble; nErrBnds: LapackInt;
                 errBndsNorm: ptr cdouble; errBndsComp: ptr cdouble;
                 nparams: LapackInt; params: ptr cdouble;
                 work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zposvxx_work", dynlib: libName.}
proc zpotrf*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexDouble;
            lda: LapackInt): LapackInt {.importc: "LAPACKE_zpotrf", dynlib: libName.}
proc zpotrf2*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexDouble;
             lda: LapackInt): LapackInt {.importc: "LAPACKE_zpotrf2", dynlib: libName.}
proc zpotrf2Work*(matrixLayout: cint; uplo: char; n: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_zpotrf2_work", dynlib: libName.}
proc zpotrfWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_zpotrf_work", dynlib: libName.}
proc zpotri*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexDouble;
            lda: LapackInt): LapackInt {.importc: "LAPACKE_zpotri", dynlib: libName.}
proc zpotriWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_zpotri_work", dynlib: libName.}
proc zpotrs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; b: ptr LapackComplexDouble;
            ldb: LapackInt): LapackInt {.importc: "LAPACKE_zpotrs", dynlib: libName.}
proc zpotrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zpotrs_work", dynlib: libName.}
proc zppcon*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexDouble;
            anorm: cdouble; rcond: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zppcon", dynlib: libName.}
proc zppconWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexDouble; anorm: cdouble; rcond: ptr cdouble;
                work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zppcon_work", dynlib: libName.}
proc zppequ*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexDouble;
            s: ptr cdouble; scond: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zppequ", dynlib: libName.}
proc zppequWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexDouble; s: ptr cdouble; scond: ptr cdouble;
                amax: ptr cdouble): LapackInt {.importc: "LAPACKE_zppequ_work",
    dynlib: libName.}
proc zpprfs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr LapackComplexDouble; afp: ptr LapackComplexDouble;
            b: ptr LapackComplexDouble; ldb: LapackInt; x: ptr LapackComplexDouble;
            ldx: LapackInt; ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zpprfs", dynlib: libName.}
proc zpprfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                ap: ptr LapackComplexDouble; afp: ptr LapackComplexDouble;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; ferr: ptr cdouble;
                berr: ptr cdouble; work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zpprfs_work", dynlib: libName.}
proc zppsv*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
           ap: ptr LapackComplexDouble; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zppsv", dynlib: libName.}
proc zppsvWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               ap: ptr LapackComplexDouble; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zppsv_work", dynlib: libName.}
proc zppsvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr LapackComplexDouble; afp: ptr LapackComplexDouble; equed: cstring;
            s: ptr cdouble; b: ptr LapackComplexDouble; ldb: LapackInt;
            x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
            ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zppsvx", dynlib: libName.}
proc zppsvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; ap: ptr LapackComplexDouble;
                afp: ptr LapackComplexDouble; equed: cstring; s: ptr cdouble;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
                ferr: ptr cdouble; berr: ptr cdouble; work: ptr LapackComplexDouble;
                rwork: ptr cdouble): LapackInt {.importc: "LAPACKE_zppsvx_work",
    dynlib: libName.}
proc zpptrf*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zpptrf", dynlib: libName.}
proc zpptrfWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zpptrf_work", dynlib: libName.}
proc zpptri*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zpptri", dynlib: libName.}
proc zpptriWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zpptri_work", dynlib: libName.}
proc zpptrs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr LapackComplexDouble; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zpptrs", dynlib: libName.}
proc zpptrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                ap: ptr LapackComplexDouble; b: ptr LapackComplexDouble;
                ldb: LapackInt): LapackInt {.importc: "LAPACKE_zpptrs_work",
    dynlib: libName.}
proc zpstrf*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexDouble;
            lda: LapackInt; piv: ptr LapackInt; rank: ptr LapackInt; tol: cdouble): LapackInt {.
    importc: "LAPACKE_zpstrf", dynlib: libName.}
proc zpstrfWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; piv: ptr LapackInt;
                rank: ptr LapackInt; tol: cdouble; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zpstrf_work", dynlib: libName.}
proc zptcon*(n: LapackInt; d: ptr cdouble; e: ptr LapackComplexDouble; anorm: cdouble;
            rcond: ptr cdouble): LapackInt {.importc: "LAPACKE_zptcon",
    dynlib: libName.}
proc zptconWork*(n: LapackInt; d: ptr cdouble; e: ptr LapackComplexDouble;
                anorm: cdouble; rcond: ptr cdouble; work: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zptcon_work", dynlib: libName.}
proc zpteqr*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cdouble;
            e: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_zpteqr", dynlib: libName.}
proc zpteqrWork*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cdouble;
                e: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt;
                work: ptr cdouble): LapackInt {.importc: "LAPACKE_zpteqr_work",
    dynlib: libName.}
proc zptrfs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            d: ptr cdouble; e: ptr LapackComplexDouble; df: ptr cdouble;
            ef: ptr LapackComplexDouble; b: ptr LapackComplexDouble; ldb: LapackInt;
            x: ptr LapackComplexDouble; ldx: LapackInt; ferr: ptr cdouble;
            berr: ptr cdouble): LapackInt {.importc: "LAPACKE_zptrfs", dynlib: libName.}
proc zptrfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                d: ptr cdouble; e: ptr LapackComplexDouble; df: ptr cdouble;
                ef: ptr LapackComplexDouble; b: ptr LapackComplexDouble;
                ldb: LapackInt; x: ptr LapackComplexDouble; ldx: LapackInt;
                ferr: ptr cdouble; berr: ptr cdouble; work: ptr LapackComplexDouble;
                rwork: ptr cdouble): LapackInt {.importc: "LAPACKE_zptrfs_work",
    dynlib: libName.}
proc zptsv*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; d: ptr cdouble;
           e: ptr LapackComplexDouble; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zptsv", dynlib: libName.}
proc zptsvWork*(matrixLayout: cint; n: LapackInt; nrhs: LapackInt; d: ptr cdouble;
               e: ptr LapackComplexDouble; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zptsv_work", dynlib: libName.}
proc zptsvx*(matrixLayout: cint; fact: char; n: LapackInt; nrhs: LapackInt;
            d: ptr cdouble; e: ptr LapackComplexDouble; df: ptr cdouble;
            ef: ptr LapackComplexDouble; b: ptr LapackComplexDouble; ldb: LapackInt;
            x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
            ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zptsvx", dynlib: libName.}
proc zptsvxWork*(matrixLayout: cint; fact: char; n: LapackInt; nrhs: LapackInt;
                d: ptr cdouble; e: ptr LapackComplexDouble; df: ptr cdouble;
                ef: ptr LapackComplexDouble; b: ptr LapackComplexDouble;
                ldb: LapackInt; x: ptr LapackComplexDouble; ldx: LapackInt;
                rcond: ptr cdouble; ferr: ptr cdouble; berr: ptr cdouble;
                work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zptsvx_work", dynlib: libName.}
proc zpttrf*(n: LapackInt; d: ptr cdouble; e: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zpttrf", dynlib: libName.}
proc zpttrfWork*(n: LapackInt; d: ptr cdouble; e: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zpttrf_work", dynlib: libName.}
proc zpttrs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            d: ptr cdouble; e: ptr LapackComplexDouble; b: ptr LapackComplexDouble;
            ldb: LapackInt): LapackInt {.importc: "LAPACKE_zpttrs", dynlib: libName.}
proc zpttrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                d: ptr cdouble; e: ptr LapackComplexDouble;
                b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zpttrs_work", dynlib: libName.}
proc zspcon*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexDouble;
            ipiv: ptr LapackInt; anorm: cdouble; rcond: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zspcon", dynlib: libName.}
proc zspconWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexDouble; ipiv: ptr LapackInt; anorm: cdouble;
                rcond: ptr cdouble; work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zspcon_work", dynlib: libName.}
proc zsprfs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr LapackComplexDouble; afp: ptr LapackComplexDouble;
            ipiv: ptr LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt;
            x: ptr LapackComplexDouble; ldx: LapackInt; ferr: ptr cdouble;
            berr: ptr cdouble): LapackInt {.importc: "LAPACKE_zsprfs", dynlib: libName.}
proc zsprfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                ap: ptr LapackComplexDouble; afp: ptr LapackComplexDouble;
                ipiv: ptr LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; ferr: ptr cdouble;
                berr: ptr cdouble; work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zsprfs_work", dynlib: libName.}
proc zspsv*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
           ap: ptr LapackComplexDouble; ipiv: ptr LapackInt;
           b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zspsv", dynlib: libName.}
proc zspsvWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               ap: ptr LapackComplexDouble; ipiv: ptr LapackInt;
               b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zspsv_work", dynlib: libName.}
proc zspsvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr LapackComplexDouble; afp: ptr LapackComplexDouble;
            ipiv: ptr LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt;
            x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
            ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zspsvx", dynlib: libName.}
proc zspsvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; ap: ptr LapackComplexDouble;
                afp: ptr LapackComplexDouble; ipiv: ptr LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
                ferr: ptr cdouble; berr: ptr cdouble; work: ptr LapackComplexDouble;
                rwork: ptr cdouble): LapackInt {.importc: "LAPACKE_zspsvx_work",
    dynlib: libName.}
proc zsptrf*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexDouble;
            ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_zsptrf",
    dynlib: libName.}
proc zsptrfWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexDouble; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zsptrf_work", dynlib: libName.}
proc zsptri*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexDouble;
            ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_zsptri",
    dynlib: libName.}
proc zsptriWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexDouble; ipiv: ptr LapackInt;
                work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zsptri_work", dynlib: libName.}
proc zsptrs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            ap: ptr LapackComplexDouble; ipiv: ptr LapackInt;
            b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zsptrs", dynlib: libName.}
proc zsptrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                ap: ptr LapackComplexDouble; ipiv: ptr LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zsptrs_work", dynlib: libName.}
proc zstedc*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cdouble;
            e: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_zstedc", dynlib: libName.}
proc zstedcWork*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cdouble;
                e: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt;
                work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble;
                lrwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zstedc_work", dynlib: libName.}
proc zstegr*(matrixLayout: cint; jobz: char; range: char; n: LapackInt; d: ptr cdouble;
            e: ptr cdouble; vl: cdouble; vu: cdouble; il: LapackInt; iu: LapackInt;
            abstol: cdouble; m: ptr LapackInt; w: ptr cdouble;
            z: ptr LapackComplexDouble; ldz: LapackInt; isuppz: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zstegr", dynlib: libName.}
proc zstegrWork*(matrixLayout: cint; jobz: char; range: char; n: LapackInt;
                d: ptr cdouble; e: ptr cdouble; vl: cdouble; vu: cdouble; il: LapackInt;
                iu: LapackInt; abstol: cdouble; m: ptr LapackInt; w: ptr cdouble;
                z: ptr LapackComplexDouble; ldz: LapackInt; isuppz: ptr LapackInt;
                work: ptr cdouble; lwork: LapackInt; iwork: ptr LapackInt;
                liwork: LapackInt): LapackInt {.importc: "LAPACKE_zstegr_work",
    dynlib: libName.}
proc zstein*(matrixLayout: cint; n: LapackInt; d: ptr cdouble; e: ptr cdouble;
            m: LapackInt; w: ptr cdouble; iblock: ptr LapackInt; isplit: ptr LapackInt;
            z: ptr LapackComplexDouble; ldz: LapackInt; ifailv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zstein", dynlib: libName.}
proc zsteinWork*(matrixLayout: cint; n: LapackInt; d: ptr cdouble; e: ptr cdouble;
                m: LapackInt; w: ptr cdouble; iblock: ptr LapackInt;
                isplit: ptr LapackInt; z: ptr LapackComplexDouble; ldz: LapackInt;
                work: ptr cdouble; iwork: ptr LapackInt; ifailv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zstein_work", dynlib: libName.}
proc zstemr*(matrixLayout: cint; jobz: char; range: char; n: LapackInt; d: ptr cdouble;
            e: ptr cdouble; vl: cdouble; vu: cdouble; il: LapackInt; iu: LapackInt;
            m: ptr LapackInt; w: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt;
            nzc: LapackInt; isuppz: ptr LapackInt; tryrac: ptr LapackLogical): LapackInt {.
    importc: "LAPACKE_zstemr", dynlib: libName.}
proc zstemrWork*(matrixLayout: cint; jobz: char; range: char; n: LapackInt;
                d: ptr cdouble; e: ptr cdouble; vl: cdouble; vu: cdouble; il: LapackInt;
                iu: LapackInt; m: ptr LapackInt; w: ptr cdouble;
                z: ptr LapackComplexDouble; ldz: LapackInt; nzc: LapackInt;
                isuppz: ptr LapackInt; tryrac: ptr LapackLogical; work: ptr cdouble;
                lwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zstemr_work", dynlib: libName.}
proc zsteqr*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cdouble;
            e: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_zsteqr", dynlib: libName.}
proc zsteqrWork*(matrixLayout: cint; compz: char; n: LapackInt; d: ptr cdouble;
                e: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt;
                work: ptr cdouble): LapackInt {.importc: "LAPACKE_zsteqr_work",
    dynlib: libName.}
proc zsycon*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexDouble;
            lda: LapackInt; ipiv: ptr LapackInt; anorm: cdouble; rcond: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zsycon", dynlib: libName.}
proc zsyconWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                anorm: cdouble; rcond: ptr cdouble; work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zsycon_work", dynlib: libName.}
proc zsyconv*(matrixLayout: cint; uplo: char; way: char; n: LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
             e: ptr LapackComplexDouble): LapackInt {.importc: "LAPACKE_zsyconv",
    dynlib: libName.}
proc zsyconvWork*(matrixLayout: cint; uplo: char; way: char; n: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                 e: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zsyconv_work", dynlib: libName.}
proc zsyequb*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexDouble;
             lda: LapackInt; s: ptr cdouble; scond: ptr cdouble; amax: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zsyequb", dynlib: libName.}
proc zsyequbWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt; s: ptr cdouble;
                 scond: ptr cdouble; amax: ptr cdouble; work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zsyequb_work", dynlib: libName.}
proc zsyr*(matrixLayout: cint; uplo: char; n: LapackInt; alpha: LapackComplexDouble;
          x: ptr LapackComplexDouble; incx: LapackInt; a: ptr LapackComplexDouble;
          lda: LapackInt): LapackInt {.importc: "LAPACKE_zsyr", dynlib: libName.}
proc zsyrWork*(matrixLayout: cint; uplo: char; n: LapackInt;
              alpha: LapackComplexDouble; x: ptr LapackComplexDouble;
              incx: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_zsyr_work", dynlib: libName.}
proc zsyrfs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; af: ptr LapackComplexDouble;
            ldaf: LapackInt; ipiv: ptr LapackInt; b: ptr LapackComplexDouble;
            ldb: LapackInt; x: ptr LapackComplexDouble; ldx: LapackInt;
            ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zsyrfs", dynlib: libName.}
proc zsyrfsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                af: ptr LapackComplexDouble; ldaf: LapackInt; ipiv: ptr LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; ferr: ptr cdouble;
                berr: ptr cdouble; work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zsyrfs_work", dynlib: libName.}
proc zsyrfsx*(matrixLayout: cint; uplo: char; equed: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; af: ptr LapackComplexDouble;
             ldaf: LapackInt; ipiv: ptr LapackInt; s: ptr cdouble;
             b: ptr LapackComplexDouble; ldb: LapackInt; x: ptr LapackComplexDouble;
             ldx: LapackInt; rcond: ptr cdouble; berr: ptr cdouble; nErrBnds: LapackInt;
             errBndsNorm: ptr cdouble; errBndsComp: ptr cdouble; nparams: LapackInt;
             params: ptr cdouble): LapackInt {.importc: "LAPACKE_zsyrfsx",
    dynlib: libName.}
proc zsyrfsxWork*(matrixLayout: cint; uplo: char; equed: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                 af: ptr LapackComplexDouble; ldaf: LapackInt; ipiv: ptr LapackInt;
                 s: ptr cdouble; b: ptr LapackComplexDouble; ldb: LapackInt;
                 x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
                 berr: ptr cdouble; nErrBnds: LapackInt; errBndsNorm: ptr cdouble;
                 errBndsComp: ptr cdouble; nparams: LapackInt; params: ptr cdouble;
                 work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zsyrfsx_work", dynlib: libName.}
proc zsysv*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
           a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
           b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zsysv", dynlib: libName.}
proc zsysvRook*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
               b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zsysv_rook", dynlib: libName.}
proc zsysvRookWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                   a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                   b: ptr LapackComplexDouble; ldb: LapackInt;
                   work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zsysv_rook_work", dynlib: libName.}
proc zsysvWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
               a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
               b: ptr LapackComplexDouble; ldb: LapackInt;
               work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zsysv_work", dynlib: libName.}
proc zsysvx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; af: ptr LapackComplexDouble;
            ldaf: LapackInt; ipiv: ptr LapackInt; b: ptr LapackComplexDouble;
            ldb: LapackInt; x: ptr LapackComplexDouble; ldx: LapackInt;
            rcond: ptr cdouble; ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zsysvx", dynlib: libName.}
proc zsysvxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                nrhs: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                af: ptr LapackComplexDouble; ldaf: LapackInt; ipiv: ptr LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; rcond: ptr cdouble;
                ferr: ptr cdouble; berr: ptr cdouble; work: ptr LapackComplexDouble;
                lwork: LapackInt; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zsysvx_work", dynlib: libName.}
proc zsysvxx*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; af: ptr LapackComplexDouble;
             ldaf: LapackInt; ipiv: ptr LapackInt; equed: cstring; s: ptr cdouble;
             b: ptr LapackComplexDouble; ldb: LapackInt; x: ptr LapackComplexDouble;
             ldx: LapackInt; rcond: ptr cdouble; rpvgrw: ptr cdouble; berr: ptr cdouble;
             nErrBnds: LapackInt; errBndsNorm: ptr cdouble; errBndsComp: ptr cdouble;
             nparams: LapackInt; params: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zsysvxx", dynlib: libName.}
proc zsysvxxWork*(matrixLayout: cint; fact: char; uplo: char; n: LapackInt;
                 nrhs: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                 af: ptr LapackComplexDouble; ldaf: LapackInt; ipiv: ptr LapackInt;
                 equed: cstring; s: ptr cdouble; b: ptr LapackComplexDouble;
                 ldb: LapackInt; x: ptr LapackComplexDouble; ldx: LapackInt;
                 rcond: ptr cdouble; rpvgrw: ptr cdouble; berr: ptr cdouble;
                 nErrBnds: LapackInt; errBndsNorm: ptr cdouble;
                 errBndsComp: ptr cdouble; nparams: LapackInt; params: ptr cdouble;
                 work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zsysvxx_work", dynlib: libName.}
proc zsyswapr*(matrixLayout: cint; uplo: char; n: LapackInt;
              a: ptr LapackComplexDouble; lda: LapackInt; i1: LapackInt; i2: LapackInt): LapackInt {.
    importc: "LAPACKE_zsyswapr", dynlib: libName.}
proc zsyswaprWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                  a: ptr LapackComplexDouble; lda: LapackInt; i1: LapackInt;
                  i2: LapackInt): LapackInt {.importc: "LAPACKE_zsyswapr_work",
    dynlib: libName.}
proc zsytrf*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexDouble;
            lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zsytrf", dynlib: libName.}
proc zsytrfRook*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zsytrf_rook", dynlib: libName.}
proc zsytrfRookWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                    a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                    work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zsytrf_rook_work", dynlib: libName.}
proc zsytrfWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zsytrf_work", dynlib: libName.}
proc zsytri*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexDouble;
            lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zsytri", dynlib: libName.}
proc zsytri2*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexDouble;
             lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zsytri2", dynlib: libName.}
proc zsytri2Work*(matrixLayout: cint; uplo: char; n: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                 work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zsytri2_work", dynlib: libName.}
proc zsytri2x*(matrixLayout: cint; uplo: char; n: LapackInt;
              a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
              nb: LapackInt): LapackInt {.importc: "LAPACKE_zsytri2x",
                                       dynlib: libName.}
proc zsytri2xWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                  a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                  work: ptr LapackComplexDouble; nb: LapackInt): LapackInt {.
    importc: "LAPACKE_zsytri2x_work", dynlib: libName.}
proc zsytriWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zsytri_work", dynlib: libName.}
proc zsytrs*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
            b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zsytrs", dynlib: libName.}
proc zsytrs2*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
             b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zsytrs2", dynlib: libName.}
proc zsytrs2Work*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                 b: ptr LapackComplexDouble; ldb: LapackInt;
                 work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zsytrs2_work", dynlib: libName.}
proc zsytrsRook*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zsytrs_rook", dynlib: libName.}
proc zsytrsRookWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                    a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                    b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zsytrs_rook_work", dynlib: libName.}
proc zsytrsWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zsytrs_work", dynlib: libName.}
proc ztbcon*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
            kd: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
            rcond: ptr cdouble): LapackInt {.importc: "LAPACKE_ztbcon",
    dynlib: libName.}
proc ztbconWork*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
                kd: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
                rcond: ptr cdouble; work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_ztbcon_work", dynlib: libName.}
proc ztbrfs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            kd: LapackInt; nrhs: LapackInt; ab: ptr LapackComplexDouble;
            ldab: LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt;
            x: ptr LapackComplexDouble; ldx: LapackInt; ferr: ptr cdouble;
            berr: ptr cdouble): LapackInt {.importc: "LAPACKE_ztbrfs", dynlib: libName.}
proc ztbrfsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                kd: LapackInt; nrhs: LapackInt; ab: ptr LapackComplexDouble;
                ldab: LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; ferr: ptr cdouble;
                berr: ptr cdouble; work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_ztbrfs_work", dynlib: libName.}
proc ztbtrs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            kd: LapackInt; nrhs: LapackInt; ab: ptr LapackComplexDouble;
            ldab: LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ztbtrs", dynlib: libName.}
proc ztbtrsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                kd: LapackInt; nrhs: LapackInt; ab: ptr LapackComplexDouble;
                ldab: LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ztbtrs_work", dynlib: libName.}
proc ztfsm*(matrixLayout: cint; transr: char; side: char; uplo: char; trans: char;
           diag: char; m: LapackInt; n: LapackInt; alpha: LapackComplexDouble;
           a: ptr LapackComplexDouble; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ztfsm", dynlib: libName.}
proc ztfsmWork*(matrixLayout: cint; transr: char; side: char; uplo: char; trans: char;
               diag: char; m: LapackInt; n: LapackInt; alpha: LapackComplexDouble;
               a: ptr LapackComplexDouble; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ztfsm_work", dynlib: libName.}
proc ztftri*(matrixLayout: cint; transr: char; uplo: char; diag: char; n: LapackInt;
            a: ptr LapackComplexDouble): LapackInt {.importc: "LAPACKE_ztftri",
    dynlib: libName.}
proc ztftriWork*(matrixLayout: cint; transr: char; uplo: char; diag: char; n: LapackInt;
                a: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_ztftri_work", dynlib: libName.}
proc ztfttp*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
            arf: ptr LapackComplexDouble; ap: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_ztfttp", dynlib: libName.}
proc ztfttpWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                arf: ptr LapackComplexDouble; ap: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_ztfttp_work", dynlib: libName.}
proc ztfttr*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
            arf: ptr LapackComplexDouble; a: ptr LapackComplexDouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_ztfttr", dynlib: libName.}
proc ztfttrWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                arf: ptr LapackComplexDouble; a: ptr LapackComplexDouble;
                lda: LapackInt): LapackInt {.importc: "LAPACKE_ztfttr_work",
    dynlib: libName.}
proc ztgevc*(matrixLayout: cint; side: char; howmny: char; select: ptr LapackLogical;
            n: LapackInt; s: ptr LapackComplexDouble; lds: LapackInt;
            p: ptr LapackComplexDouble; ldp: LapackInt; vl: ptr LapackComplexDouble;
            ldvl: LapackInt; vr: ptr LapackComplexDouble; ldvr: LapackInt;
            mm: LapackInt; m: ptr LapackInt): LapackInt {.importc: "LAPACKE_ztgevc",
    dynlib: libName.}
proc ztgevcWork*(matrixLayout: cint; side: char; howmny: char;
                select: ptr LapackLogical; n: LapackInt; s: ptr LapackComplexDouble;
                lds: LapackInt; p: ptr LapackComplexDouble; ldp: LapackInt;
                vl: ptr LapackComplexDouble; ldvl: LapackInt;
                vr: ptr LapackComplexDouble; ldvr: LapackInt; mm: LapackInt;
                m: ptr LapackInt; work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_ztgevc_work", dynlib: libName.}
proc ztgexc*(matrixLayout: cint; wantq: LapackLogical; wantz: LapackLogical;
            n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            b: ptr LapackComplexDouble; ldb: LapackInt; q: ptr LapackComplexDouble;
            ldq: LapackInt; z: ptr LapackComplexDouble; ldz: LapackInt;
            ifst: LapackInt; ilst: LapackInt): LapackInt {.importc: "LAPACKE_ztgexc",
    dynlib: libName.}
proc ztgexcWork*(matrixLayout: cint; wantq: LapackLogical; wantz: LapackLogical;
                n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                q: ptr LapackComplexDouble; ldq: LapackInt;
                z: ptr LapackComplexDouble; ldz: LapackInt; ifst: LapackInt;
                ilst: LapackInt): LapackInt {.importc: "LAPACKE_ztgexc_work",
    dynlib: libName.}
proc ztgsen*(matrixLayout: cint; ijob: LapackInt; wantq: LapackLogical;
            wantz: LapackLogical; select: ptr LapackLogical; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; b: ptr LapackComplexDouble;
            ldb: LapackInt; alpha: ptr LapackComplexDouble;
            beta: ptr LapackComplexDouble; q: ptr LapackComplexDouble; ldq: LapackInt;
            z: ptr LapackComplexDouble; ldz: LapackInt; m: ptr LapackInt;
            pl: ptr cdouble; pr: ptr cdouble; dif: ptr cdouble): LapackInt {.
    importc: "LAPACKE_ztgsen", dynlib: libName.}
proc ztgsenWork*(matrixLayout: cint; ijob: LapackInt; wantq: LapackLogical;
                wantz: LapackLogical; select: ptr LapackLogical; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                alpha: ptr LapackComplexDouble; beta: ptr LapackComplexDouble;
                q: ptr LapackComplexDouble; ldq: LapackInt;
                z: ptr LapackComplexDouble; ldz: LapackInt; m: ptr LapackInt;
                pl: ptr cdouble; pr: ptr cdouble; dif: ptr cdouble;
                work: ptr LapackComplexDouble; lwork: LapackInt;
                iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ztgsen_work", dynlib: libName.}
proc ztgsja*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
            p: LapackInt; n: LapackInt; k: LapackInt; l: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; b: ptr LapackComplexDouble;
            ldb: LapackInt; tola: cdouble; tolb: cdouble; alpha: ptr cdouble;
            beta: ptr cdouble; u: ptr LapackComplexDouble; ldu: LapackInt;
            v: ptr LapackComplexDouble; ldv: LapackInt; q: ptr LapackComplexDouble;
            ldq: LapackInt; ncycle: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ztgsja", dynlib: libName.}
proc ztgsjaWork*(matrixLayout: cint; jobu: char; jobv: char; jobq: char; m: LapackInt;
                p: LapackInt; n: LapackInt; k: LapackInt; l: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt; tola: cdouble;
                tolb: cdouble; alpha: ptr cdouble; beta: ptr cdouble;
                u: ptr LapackComplexDouble; ldu: LapackInt;
                v: ptr LapackComplexDouble; ldv: LapackInt;
                q: ptr LapackComplexDouble; ldq: LapackInt;
                work: ptr LapackComplexDouble; ncycle: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ztgsja_work", dynlib: libName.}
proc ztgsna*(matrixLayout: cint; job: char; howmny: char; select: ptr LapackLogical;
            n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            b: ptr LapackComplexDouble; ldb: LapackInt; vl: ptr LapackComplexDouble;
            ldvl: LapackInt; vr: ptr LapackComplexDouble; ldvr: LapackInt;
            s: ptr cdouble; dif: ptr cdouble; mm: LapackInt; m: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ztgsna", dynlib: libName.}
proc ztgsnaWork*(matrixLayout: cint; job: char; howmny: char;
                select: ptr LapackLogical; n: LapackInt; a: ptr LapackComplexDouble;
                lda: LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt;
                vl: ptr LapackComplexDouble; ldvl: LapackInt;
                vr: ptr LapackComplexDouble; ldvr: LapackInt; s: ptr cdouble;
                dif: ptr cdouble; mm: LapackInt; m: ptr LapackInt;
                work: ptr LapackComplexDouble; lwork: LapackInt; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ztgsna_work", dynlib: libName.}
proc ztgsyl*(matrixLayout: cint; trans: char; ijob: LapackInt; m: LapackInt;
            n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            b: ptr LapackComplexDouble; ldb: LapackInt; c: ptr LapackComplexDouble;
            ldc: LapackInt; d: ptr LapackComplexDouble; ldd: LapackInt;
            e: ptr LapackComplexDouble; lde: LapackInt; f: ptr LapackComplexDouble;
            ldf: LapackInt; scale: ptr cdouble; dif: ptr cdouble): LapackInt {.
    importc: "LAPACKE_ztgsyl", dynlib: libName.}
proc ztgsylWork*(matrixLayout: cint; trans: char; ijob: LapackInt; m: LapackInt;
                n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                c: ptr LapackComplexDouble; ldc: LapackInt;
                d: ptr LapackComplexDouble; ldd: LapackInt;
                e: ptr LapackComplexDouble; lde: LapackInt;
                f: ptr LapackComplexDouble; ldf: LapackInt; scale: ptr cdouble;
                dif: ptr cdouble; work: ptr LapackComplexDouble; lwork: LapackInt;
                iwork: ptr LapackInt): LapackInt {.importc: "LAPACKE_ztgsyl_work",
    dynlib: libName.}
proc ztpcon*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
            ap: ptr LapackComplexDouble; rcond: ptr cdouble): LapackInt {.
    importc: "LAPACKE_ztpcon", dynlib: libName.}
proc ztpconWork*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
                ap: ptr LapackComplexDouble; rcond: ptr cdouble;
                work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_ztpcon_work", dynlib: libName.}
proc ztpmqrt*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
             k: LapackInt; l: LapackInt; nb: LapackInt; v: ptr LapackComplexDouble;
             ldv: LapackInt; t: ptr LapackComplexDouble; ldt: LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; b: ptr LapackComplexDouble;
             ldb: LapackInt): LapackInt {.importc: "LAPACKE_ztpmqrt", dynlib: libName.}
proc ztpmqrtWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt;
                 n: LapackInt; k: LapackInt; l: LapackInt; nb: LapackInt;
                 v: ptr LapackComplexDouble; ldv: LapackInt;
                 t: ptr LapackComplexDouble; ldt: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt;
                 b: ptr LapackComplexDouble; ldb: LapackInt;
                 work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_ztpmqrt_work", dynlib: libName.}
proc ztpqrt*(matrixLayout: cint; m: LapackInt; n: LapackInt; l: LapackInt; nb: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; b: ptr LapackComplexDouble;
            ldb: LapackInt; t: ptr LapackComplexDouble; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_ztpqrt", dynlib: libName.}
proc ztpqrt2*(matrixLayout: cint; m: LapackInt; n: LapackInt; l: LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; b: ptr LapackComplexDouble;
             ldb: LapackInt; t: ptr LapackComplexDouble; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_ztpqrt2", dynlib: libName.}
proc ztpqrt2Work*(matrixLayout: cint; m: LapackInt; n: LapackInt; l: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt;
                 b: ptr LapackComplexDouble; ldb: LapackInt;
                 t: ptr LapackComplexDouble; ldt: LapackInt): LapackInt {.
    importc: "LAPACKE_ztpqrt2_work", dynlib: libName.}
proc ztpqrtWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; l: LapackInt;
                nb: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                t: ptr LapackComplexDouble; ldt: LapackInt;
                work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_ztpqrt_work", dynlib: libName.}
proc ztprfb*(matrixLayout: cint; side: char; trans: char; direct: char; storev: char;
            m: LapackInt; n: LapackInt; k: LapackInt; l: LapackInt;
            v: ptr LapackComplexDouble; ldv: LapackInt; t: ptr LapackComplexDouble;
            ldt: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ztprfb", dynlib: libName.}
proc ztprfbWork*(matrixLayout: cint; side: char; trans: char; direct: char; storev: char;
                m: LapackInt; n: LapackInt; k: LapackInt; l: LapackInt;
                v: ptr LapackComplexDouble; ldv: LapackInt;
                t: ptr LapackComplexDouble; ldt: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                work: ptr LapackComplexDouble; ldwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ztprfb_work", dynlib: libName.}
proc ztprfs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            nrhs: LapackInt; ap: ptr LapackComplexDouble; b: ptr LapackComplexDouble;
            ldb: LapackInt; x: ptr LapackComplexDouble; ldx: LapackInt;
            ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_ztprfs", dynlib: libName.}
proc ztprfsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                nrhs: LapackInt; ap: ptr LapackComplexDouble;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; ferr: ptr cdouble;
                berr: ptr cdouble; work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_ztprfs_work", dynlib: libName.}
proc ztptri*(matrixLayout: cint; uplo: char; diag: char; n: LapackInt;
            ap: ptr LapackComplexDouble): LapackInt {.importc: "LAPACKE_ztptri",
    dynlib: libName.}
proc ztptriWork*(matrixLayout: cint; uplo: char; diag: char; n: LapackInt;
                ap: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_ztptri_work", dynlib: libName.}
proc ztptrs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            nrhs: LapackInt; ap: ptr LapackComplexDouble; b: ptr LapackComplexDouble;
            ldb: LapackInt): LapackInt {.importc: "LAPACKE_ztptrs", dynlib: libName.}
proc ztptrsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                nrhs: LapackInt; ap: ptr LapackComplexDouble;
                b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ztptrs_work", dynlib: libName.}
proc ztpttf*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
            ap: ptr LapackComplexDouble; arf: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_ztpttf", dynlib: libName.}
proc ztpttfWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                ap: ptr LapackComplexDouble; arf: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_ztpttf_work", dynlib: libName.}
proc ztpttr*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexDouble;
            a: ptr LapackComplexDouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_ztpttr", dynlib: libName.}
proc ztpttrWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexDouble; a: ptr LapackComplexDouble;
                lda: LapackInt): LapackInt {.importc: "LAPACKE_ztpttr_work",
    dynlib: libName.}
proc ztrcon*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; rcond: ptr cdouble): LapackInt {.
    importc: "LAPACKE_ztrcon", dynlib: libName.}
proc ztrconWork*(matrixLayout: cint; norm: char; uplo: char; diag: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt; rcond: ptr cdouble;
                work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_ztrcon_work", dynlib: libName.}
proc ztrevc*(matrixLayout: cint; side: char; howmny: char; select: ptr LapackLogical;
            n: LapackInt; t: ptr LapackComplexDouble; ldt: LapackInt;
            vl: ptr LapackComplexDouble; ldvl: LapackInt;
            vr: ptr LapackComplexDouble; ldvr: LapackInt; mm: LapackInt;
            m: ptr LapackInt): LapackInt {.importc: "LAPACKE_ztrevc", dynlib: libName.}
proc ztrevcWork*(matrixLayout: cint; side: char; howmny: char;
                select: ptr LapackLogical; n: LapackInt; t: ptr LapackComplexDouble;
                ldt: LapackInt; vl: ptr LapackComplexDouble; ldvl: LapackInt;
                vr: ptr LapackComplexDouble; ldvr: LapackInt; mm: LapackInt;
                m: ptr LapackInt; work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_ztrevc_work", dynlib: libName.}
proc ztrexc*(matrixLayout: cint; compq: char; n: LapackInt; t: ptr LapackComplexDouble;
            ldt: LapackInt; q: ptr LapackComplexDouble; ldq: LapackInt;
            ifst: LapackInt; ilst: LapackInt): LapackInt {.importc: "LAPACKE_ztrexc",
    dynlib: libName.}
proc ztrexcWork*(matrixLayout: cint; compq: char; n: LapackInt;
                t: ptr LapackComplexDouble; ldt: LapackInt;
                q: ptr LapackComplexDouble; ldq: LapackInt; ifst: LapackInt;
                ilst: LapackInt): LapackInt {.importc: "LAPACKE_ztrexc_work",
    dynlib: libName.}
proc ztrrfs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            nrhs: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            b: ptr LapackComplexDouble; ldb: LapackInt; x: ptr LapackComplexDouble;
            ldx: LapackInt; ferr: ptr cdouble; berr: ptr cdouble): LapackInt {.
    importc: "LAPACKE_ztrrfs", dynlib: libName.}
proc ztrrfsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                nrhs: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                x: ptr LapackComplexDouble; ldx: LapackInt; ferr: ptr cdouble;
                berr: ptr cdouble; work: ptr LapackComplexDouble; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_ztrrfs_work", dynlib: libName.}
proc ztrsen*(matrixLayout: cint; job: char; compq: char; select: ptr LapackLogical;
            n: LapackInt; t: ptr LapackComplexDouble; ldt: LapackInt;
            q: ptr LapackComplexDouble; ldq: LapackInt; w: ptr LapackComplexDouble;
            m: ptr LapackInt; s: ptr cdouble; sep: ptr cdouble): LapackInt {.
    importc: "LAPACKE_ztrsen", dynlib: libName.}
proc ztrsenWork*(matrixLayout: cint; job: char; compq: char; select: ptr LapackLogical;
                n: LapackInt; t: ptr LapackComplexDouble; ldt: LapackInt;
                q: ptr LapackComplexDouble; ldq: LapackInt;
                w: ptr LapackComplexDouble; m: ptr LapackInt; s: ptr cdouble;
                sep: ptr cdouble; work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ztrsen_work", dynlib: libName.}
proc ztrsna*(matrixLayout: cint; job: char; howmny: char; select: ptr LapackLogical;
            n: LapackInt; t: ptr LapackComplexDouble; ldt: LapackInt;
            vl: ptr LapackComplexDouble; ldvl: LapackInt;
            vr: ptr LapackComplexDouble; ldvr: LapackInt; s: ptr cdouble;
            sep: ptr cdouble; mm: LapackInt; m: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ztrsna", dynlib: libName.}
proc ztrsnaWork*(matrixLayout: cint; job: char; howmny: char;
                select: ptr LapackLogical; n: LapackInt; t: ptr LapackComplexDouble;
                ldt: LapackInt; vl: ptr LapackComplexDouble; ldvl: LapackInt;
                vr: ptr LapackComplexDouble; ldvr: LapackInt; s: ptr cdouble;
                sep: ptr cdouble; mm: LapackInt; m: ptr LapackInt;
                work: ptr LapackComplexDouble; ldwork: LapackInt; rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_ztrsna_work", dynlib: libName.}
proc ztrsyl*(matrixLayout: cint; trana: char; tranb: char; isgn: LapackInt; m: LapackInt;
            n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            b: ptr LapackComplexDouble; ldb: LapackInt; c: ptr LapackComplexDouble;
            ldc: LapackInt; scale: ptr cdouble): LapackInt {.
    importc: "LAPACKE_ztrsyl", dynlib: libName.}
proc ztrsylWork*(matrixLayout: cint; trana: char; tranb: char; isgn: LapackInt;
                m: LapackInt; n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt;
                c: ptr LapackComplexDouble; ldc: LapackInt; scale: ptr cdouble): LapackInt {.
    importc: "LAPACKE_ztrsyl_work", dynlib: libName.}
proc ztrtri*(matrixLayout: cint; uplo: char; diag: char; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_ztrtri", dynlib: libName.}
proc ztrtriWork*(matrixLayout: cint; uplo: char; diag: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_ztrtri_work", dynlib: libName.}
proc ztrtrs*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
            nrhs: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ztrtrs", dynlib: libName.}
proc ztrtrsWork*(matrixLayout: cint; uplo: char; trans: char; diag: char; n: LapackInt;
                nrhs: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ztrtrs_work", dynlib: libName.}
proc ztrttf*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; arf: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_ztrttf", dynlib: libName.}
proc ztrttfWork*(matrixLayout: cint; transr: char; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                arf: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_ztrttf_work", dynlib: libName.}
proc ztrttp*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexDouble;
            lda: LapackInt; ap: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_ztrttp", dynlib: libName.}
proc ztrttpWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                ap: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_ztrttp_work", dynlib: libName.}
proc ztzrzf*(matrixLayout: cint; m: LapackInt; n: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; tau: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_ztzrzf", dynlib: libName.}
proc ztzrzfWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                tau: ptr LapackComplexDouble; work: ptr LapackComplexDouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_ztzrzf_work",
    dynlib: libName.}
proc zunbdb*(matrixLayout: cint; trans: char; signs: char; m: LapackInt; p: LapackInt;
            q: LapackInt; x11: ptr LapackComplexDouble; ldx11: LapackInt;
            x12: ptr LapackComplexDouble; ldx12: LapackInt;
            x21: ptr LapackComplexDouble; ldx21: LapackInt;
            x22: ptr LapackComplexDouble; ldx22: LapackInt; theta: ptr cdouble;
            phi: ptr cdouble; taup1: ptr LapackComplexDouble;
            taup2: ptr LapackComplexDouble; tauq1: ptr LapackComplexDouble;
            tauq2: ptr LapackComplexDouble): LapackInt {.importc: "LAPACKE_zunbdb",
    dynlib: libName.}
proc zunbdbWork*(matrixLayout: cint; trans: char; signs: char; m: LapackInt;
                p: LapackInt; q: LapackInt; x11: ptr LapackComplexDouble;
                ldx11: LapackInt; x12: ptr LapackComplexDouble; ldx12: LapackInt;
                x21: ptr LapackComplexDouble; ldx21: LapackInt;
                x22: ptr LapackComplexDouble; ldx22: LapackInt; theta: ptr cdouble;
                phi: ptr cdouble; taup1: ptr LapackComplexDouble;
                taup2: ptr LapackComplexDouble; tauq1: ptr LapackComplexDouble;
                tauq2: ptr LapackComplexDouble; work: ptr LapackComplexDouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_zunbdb_work",
    dynlib: libName.}
proc zuncsd*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char; jobv2t: char;
            trans: char; signs: char; m: LapackInt; p: LapackInt; q: LapackInt;
            x11: ptr LapackComplexDouble; ldx11: LapackInt;
            x12: ptr LapackComplexDouble; ldx12: LapackInt;
            x21: ptr LapackComplexDouble; ldx21: LapackInt;
            x22: ptr LapackComplexDouble; ldx22: LapackInt; theta: ptr cdouble;
            u1: ptr LapackComplexDouble; ldu1: LapackInt;
            u2: ptr LapackComplexDouble; ldu2: LapackInt;
            v1t: ptr LapackComplexDouble; ldv1t: LapackInt;
            v2t: ptr LapackComplexDouble; ldv2t: LapackInt): LapackInt {.
    importc: "LAPACKE_zuncsd", dynlib: libName.}
proc zuncsd2by1*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char;
                m: LapackInt; p: LapackInt; q: LapackInt;
                x11: ptr LapackComplexDouble; ldx11: LapackInt;
                x21: ptr LapackComplexDouble; ldx21: LapackInt; theta: ptr cdouble;
                u1: ptr LapackComplexDouble; ldu1: LapackInt;
                u2: ptr LapackComplexDouble; ldu2: LapackInt;
                v1t: ptr LapackComplexDouble; ldv1t: LapackInt): LapackInt {.
    importc: "LAPACKE_zuncsd2by1", dynlib: libName.}
proc zuncsd2by1Work*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char;
                    m: LapackInt; p: LapackInt; q: LapackInt;
                    x11: ptr LapackComplexDouble; ldx11: LapackInt;
                    x21: ptr LapackComplexDouble; ldx21: LapackInt;
                    theta: ptr cdouble; u1: ptr LapackComplexDouble; ldu1: LapackInt;
                    u2: ptr LapackComplexDouble; ldu2: LapackInt;
                    v1t: ptr LapackComplexDouble; ldv1t: LapackInt;
                    work: ptr LapackComplexDouble; lwork: LapackInt;
                    rwork: ptr cdouble; lrwork: LapackInt; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zuncsd2by1_work", dynlib: libName.}
proc zuncsdWork*(matrixLayout: cint; jobu1: char; jobu2: char; jobv1t: char;
                jobv2t: char; trans: char; signs: char; m: LapackInt; p: LapackInt;
                q: LapackInt; x11: ptr LapackComplexDouble; ldx11: LapackInt;
                x12: ptr LapackComplexDouble; ldx12: LapackInt;
                x21: ptr LapackComplexDouble; ldx21: LapackInt;
                x22: ptr LapackComplexDouble; ldx22: LapackInt; theta: ptr cdouble;
                u1: ptr LapackComplexDouble; ldu1: LapackInt;
                u2: ptr LapackComplexDouble; ldu2: LapackInt;
                v1t: ptr LapackComplexDouble; ldv1t: LapackInt;
                v2t: ptr LapackComplexDouble; ldv2t: LapackInt;
                work: ptr LapackComplexDouble; lwork: LapackInt; rwork: ptr cdouble;
                lrwork: LapackInt; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zuncsd_work", dynlib: libName.}
proc zungbr*(matrixLayout: cint; vect: char; m: LapackInt; n: LapackInt; k: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; tau: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zungbr", dynlib: libName.}
proc zungbrWork*(matrixLayout: cint; vect: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                tau: ptr LapackComplexDouble; work: ptr LapackComplexDouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_zungbr_work",
    dynlib: libName.}
proc zunghr*(matrixLayout: cint; n: LapackInt; ilo: LapackInt; ihi: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; tau: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zunghr", dynlib: libName.}
proc zunghrWork*(matrixLayout: cint; n: LapackInt; ilo: LapackInt; ihi: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                tau: ptr LapackComplexDouble; work: ptr LapackComplexDouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_zunghr_work",
    dynlib: libName.}
proc zunglq*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; tau: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zunglq", dynlib: libName.}
proc zunglqWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                tau: ptr LapackComplexDouble; work: ptr LapackComplexDouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_zunglq_work",
    dynlib: libName.}
proc zungql*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; tau: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zungql", dynlib: libName.}
proc zungqlWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                tau: ptr LapackComplexDouble; work: ptr LapackComplexDouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_zungql_work",
    dynlib: libName.}
proc zungqr*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; tau: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zungqr", dynlib: libName.}
proc zungqrWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                tau: ptr LapackComplexDouble; work: ptr LapackComplexDouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_zungqr_work",
    dynlib: libName.}
proc zungrq*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
            a: ptr LapackComplexDouble; lda: LapackInt; tau: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zungrq", dynlib: libName.}
proc zungrqWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; k: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                tau: ptr LapackComplexDouble; work: ptr LapackComplexDouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_zungrq_work",
    dynlib: libName.}
proc zungtr*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexDouble;
            lda: LapackInt; tau: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zungtr", dynlib: libName.}
proc zungtrWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                a: ptr LapackComplexDouble; lda: LapackInt;
                tau: ptr LapackComplexDouble; work: ptr LapackComplexDouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_zungtr_work",
    dynlib: libName.}
proc zunmbr*(matrixLayout: cint; vect: char; side: char; trans: char; m: LapackInt;
            n: LapackInt; k: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            tau: ptr LapackComplexDouble; c: ptr LapackComplexDouble; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_zunmbr", dynlib: libName.}
proc zunmbrWork*(matrixLayout: cint; vect: char; side: char; trans: char; m: LapackInt;
                n: LapackInt; k: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                tau: ptr LapackComplexDouble; c: ptr LapackComplexDouble;
                ldc: LapackInt; work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zunmbr_work", dynlib: libName.}
proc zunmhr*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            ilo: LapackInt; ihi: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            tau: ptr LapackComplexDouble; c: ptr LapackComplexDouble; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_zunmhr", dynlib: libName.}
proc zunmhrWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                ilo: LapackInt; ihi: LapackInt; a: ptr LapackComplexDouble;
                lda: LapackInt; tau: ptr LapackComplexDouble;
                c: ptr LapackComplexDouble; ldc: LapackInt;
                work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zunmhr_work", dynlib: libName.}
proc zunmlq*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            tau: ptr LapackComplexDouble; c: ptr LapackComplexDouble; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_zunmlq", dynlib: libName.}
proc zunmlqWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                tau: ptr LapackComplexDouble; c: ptr LapackComplexDouble;
                ldc: LapackInt; work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zunmlq_work", dynlib: libName.}
proc zunmql*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            tau: ptr LapackComplexDouble; c: ptr LapackComplexDouble; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_zunmql", dynlib: libName.}
proc zunmqlWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                tau: ptr LapackComplexDouble; c: ptr LapackComplexDouble;
                ldc: LapackInt; work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zunmql_work", dynlib: libName.}
proc zunmqr*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            tau: ptr LapackComplexDouble; c: ptr LapackComplexDouble; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_zunmqr", dynlib: libName.}
proc zunmqrWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                tau: ptr LapackComplexDouble; c: ptr LapackComplexDouble;
                ldc: LapackInt; work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zunmqr_work", dynlib: libName.}
proc zunmrq*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            tau: ptr LapackComplexDouble; c: ptr LapackComplexDouble; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_zunmrq", dynlib: libName.}
proc zunmrqWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                tau: ptr LapackComplexDouble; c: ptr LapackComplexDouble;
                ldc: LapackInt; work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zunmrq_work", dynlib: libName.}
proc zunmrz*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; l: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            tau: ptr LapackComplexDouble; c: ptr LapackComplexDouble; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_zunmrz", dynlib: libName.}
proc zunmrzWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; l: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                tau: ptr LapackComplexDouble; c: ptr LapackComplexDouble;
                ldc: LapackInt; work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zunmrz_work", dynlib: libName.}
proc zunmtr*(matrixLayout: cint; side: char; uplo: char; trans: char; m: LapackInt;
            n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            tau: ptr LapackComplexDouble; c: ptr LapackComplexDouble; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_zunmtr", dynlib: libName.}
proc zunmtrWork*(matrixLayout: cint; side: char; uplo: char; trans: char; m: LapackInt;
                n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                tau: ptr LapackComplexDouble; c: ptr LapackComplexDouble;
                ldc: LapackInt; work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zunmtr_work", dynlib: libName.}
proc zupgtr*(matrixLayout: cint; uplo: char; n: LapackInt; ap: ptr LapackComplexDouble;
            tau: ptr LapackComplexDouble; q: ptr LapackComplexDouble; ldq: LapackInt): LapackInt {.
    importc: "LAPACKE_zupgtr", dynlib: libName.}
proc zupgtrWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                ap: ptr LapackComplexDouble; tau: ptr LapackComplexDouble;
                q: ptr LapackComplexDouble; ldq: LapackInt;
                work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zupgtr_work", dynlib: libName.}
proc zupmtr*(matrixLayout: cint; side: char; uplo: char; trans: char; m: LapackInt;
            n: LapackInt; ap: ptr LapackComplexDouble; tau: ptr LapackComplexDouble;
            c: ptr LapackComplexDouble; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_zupmtr", dynlib: libName.}
proc zupmtrWork*(matrixLayout: cint; side: char; uplo: char; trans: char; m: LapackInt;
                n: LapackInt; ap: ptr LapackComplexDouble;
                tau: ptr LapackComplexDouble; c: ptr LapackComplexDouble;
                ldc: LapackInt; work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zupmtr_work", dynlib: libName.}
proc ilaver*(versMajor: ptr LapackInt; versMinor: ptr LapackInt;
            versPatch: ptr LapackInt) {.importc: "LAPACKE_ilaver", dynlib: libName.}
proc xerbla*(name: cstring; info: LapackInt) {.importc: "LAPACKE_xerbla",
    dynlib: libName.}
proc chesvAa*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
             b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_chesv_aa", dynlib: libName.}
proc chesvAaWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                 b: ptr LapackComplexFloat; ldb: LapackInt;
                 work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_chesv_aa_work", dynlib: libName.}
proc chetrfAa*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
              lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_chetrf_aa", dynlib: libName.}
proc chetrfAaWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                  a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                  work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_chetrf_aa_work", dynlib: libName.}
proc chetrsAa*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
              a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
              b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_chetrs_aa", dynlib: libName.}
proc chetrsAaWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                  a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                  b: ptr LapackComplexFloat; ldb: LapackInt;
                  work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_chetrs_aa_work", dynlib: libName.}
proc dsysvAa*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr cdouble; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cdouble;
             ldb: LapackInt): LapackInt {.importc: "LAPACKE_dsysv_aa",
                                       dynlib: libName.}
proc dsysvAaWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr cdouble; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cdouble;
                 ldb: LapackInt; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dsysv_aa_work", dynlib: libName.}
proc dsytrfAa*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
              lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsytrf_aa", dynlib: libName.}
proc dsytrfAaWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                  lda: LapackInt; ipiv: ptr LapackInt; work: ptr cdouble;
                  lwork: LapackInt): LapackInt {.importc: "LAPACKE_dsytrf_aa_work",
    dynlib: libName.}
proc dsytrsAa*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
              a: ptr cdouble; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cdouble;
              ldb: LapackInt): LapackInt {.importc: "LAPACKE_dsytrs_aa",
                                        dynlib: libName.}
proc dsytrsAaWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                  a: ptr cdouble; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cdouble;
                  ldb: LapackInt; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dsytrs_aa_work", dynlib: libName.}
proc ssysvAa*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr cfloat; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cfloat;
             ldb: LapackInt): LapackInt {.importc: "LAPACKE_ssysv_aa",
                                       dynlib: libName.}
proc ssysvAaWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr cfloat; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cfloat;
                 ldb: LapackInt; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ssysv_aa_work", dynlib: libName.}
proc ssytrfAa*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
              lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ssytrf_aa", dynlib: libName.}
proc ssytrfAaWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                  lda: LapackInt; ipiv: ptr LapackInt; work: ptr cfloat;
                  lwork: LapackInt): LapackInt {.importc: "LAPACKE_ssytrf_aa_work",
    dynlib: libName.}
proc ssytrsAa*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
              a: ptr cfloat; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cfloat;
              ldb: LapackInt): LapackInt {.importc: "LAPACKE_ssytrs_aa",
                                        dynlib: libName.}
proc ssytrsAaWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                  a: ptr cfloat; lda: LapackInt; ipiv: ptr LapackInt; b: ptr cfloat;
                  ldb: LapackInt; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ssytrs_aa_work", dynlib: libName.}
proc zhesvAa*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
             b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zhesv_aa", dynlib: libName.}
proc zhesvAaWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                 b: ptr LapackComplexDouble; ldb: LapackInt;
                 work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zhesv_aa_work", dynlib: libName.}
proc zhetrfAa*(matrixLayout: cint; uplo: char; n: LapackInt;
              a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zhetrf_aa", dynlib: libName.}
proc zhetrfAaWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                  a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                  work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zhetrf_aa_work", dynlib: libName.}
proc zhetrsAa*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
              a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
              b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zhetrs_aa", dynlib: libName.}
proc zhetrsAaWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                  a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                  b: ptr LapackComplexDouble; ldb: LapackInt;
                  work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zhetrs_aa_work", dynlib: libName.}
proc cgemqr*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            t: ptr LapackComplexFloat; tsize: LapackInt; c: ptr LapackComplexFloat;
            ldc: LapackInt): LapackInt {.importc: "LAPACKE_cgemqr", dynlib: libName.}
proc cgemqrWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                t: ptr LapackComplexFloat; tsize: LapackInt;
                c: ptr LapackComplexFloat; ldc: LapackInt;
                work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cgemqr_work", dynlib: libName.}
proc cgeqr*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr LapackComplexFloat;
           lda: LapackInt; t: ptr LapackComplexFloat; tsize: LapackInt): LapackInt {.
    importc: "LAPACKE_cgeqr", dynlib: libName.}
proc cgeqrWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
               a: ptr LapackComplexFloat; lda: LapackInt; t: ptr LapackComplexFloat;
               tsize: LapackInt; work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cgeqr_work", dynlib: libName.}
proc dgemqr*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr cdouble; lda: LapackInt; t: ptr cdouble;
            tsize: LapackInt; c: ptr cdouble; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_dgemqr", dynlib: libName.}
proc dgemqrWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr cdouble; lda: LapackInt; t: ptr cdouble;
                tsize: LapackInt; c: ptr cdouble; ldc: LapackInt; work: ptr cdouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_dgemqr_work",
    dynlib: libName.}
proc dgeqr*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
           lda: LapackInt; t: ptr cdouble; tsize: LapackInt): LapackInt {.
    importc: "LAPACKE_dgeqr", dynlib: libName.}
proc dgeqrWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
               lda: LapackInt; t: ptr cdouble; tsize: LapackInt; work: ptr cdouble;
               lwork: LapackInt): LapackInt {.importc: "LAPACKE_dgeqr_work",
    dynlib: libName.}
proc sgemqr*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr cfloat; lda: LapackInt; t: ptr cfloat; tsize: LapackInt;
            c: ptr cfloat; ldc: LapackInt): LapackInt {.importc: "LAPACKE_sgemqr",
    dynlib: libName.}
proc sgemqrWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr cfloat; lda: LapackInt; t: ptr cfloat;
                tsize: LapackInt; c: ptr cfloat; ldc: LapackInt; work: ptr cfloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_sgemqr_work",
    dynlib: libName.}
proc sgeqr*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat; lda: LapackInt;
           t: ptr cfloat; tsize: LapackInt): LapackInt {.importc: "LAPACKE_sgeqr",
    dynlib: libName.}
proc sgeqrWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
               lda: LapackInt; t: ptr cfloat; tsize: LapackInt; work: ptr cfloat;
               lwork: LapackInt): LapackInt {.importc: "LAPACKE_sgeqr_work",
    dynlib: libName.}
proc zgemqr*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            t: ptr LapackComplexDouble; tsize: LapackInt; c: ptr LapackComplexDouble;
            ldc: LapackInt): LapackInt {.importc: "LAPACKE_zgemqr", dynlib: libName.}
proc zgemqrWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                t: ptr LapackComplexDouble; tsize: LapackInt;
                c: ptr LapackComplexDouble; ldc: LapackInt;
                work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zgemqr_work", dynlib: libName.}
proc zgeqr*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr LapackComplexDouble;
           lda: LapackInt; t: ptr LapackComplexDouble; tsize: LapackInt): LapackInt {.
    importc: "LAPACKE_zgeqr", dynlib: libName.}
proc zgeqrWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
               a: ptr LapackComplexDouble; lda: LapackInt;
               t: ptr LapackComplexDouble; tsize: LapackInt;
               work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zgeqr_work", dynlib: libName.}
proc cgelq*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr LapackComplexFloat;
           lda: LapackInt; t: ptr LapackComplexFloat; tsize: LapackInt): LapackInt {.
    importc: "LAPACKE_cgelq", dynlib: libName.}
proc cgelqWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
               a: ptr LapackComplexFloat; lda: LapackInt; t: ptr LapackComplexFloat;
               tsize: LapackInt; work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cgelq_work", dynlib: libName.}
proc cgemlq*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
            t: ptr LapackComplexFloat; tsize: LapackInt; c: ptr LapackComplexFloat;
            ldc: LapackInt): LapackInt {.importc: "LAPACKE_cgemlq", dynlib: libName.}
proc cgemlqWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                t: ptr LapackComplexFloat; tsize: LapackInt;
                c: ptr LapackComplexFloat; ldc: LapackInt;
                work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cgemlq_work", dynlib: libName.}
proc dgelq*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
           lda: LapackInt; t: ptr cdouble; tsize: LapackInt): LapackInt {.
    importc: "LAPACKE_dgelq", dynlib: libName.}
proc dgelqWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
               lda: LapackInt; t: ptr cdouble; tsize: LapackInt; work: ptr cdouble;
               lwork: LapackInt): LapackInt {.importc: "LAPACKE_dgelq_work",
    dynlib: libName.}
proc dgemlq*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr cdouble; lda: LapackInt; t: ptr cdouble;
            tsize: LapackInt; c: ptr cdouble; ldc: LapackInt): LapackInt {.
    importc: "LAPACKE_dgemlq", dynlib: libName.}
proc dgemlqWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr cdouble; lda: LapackInt; t: ptr cdouble;
                tsize: LapackInt; c: ptr cdouble; ldc: LapackInt; work: ptr cdouble;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_dgemlq_work",
    dynlib: libName.}
proc sgelq*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat; lda: LapackInt;
           t: ptr cfloat; tsize: LapackInt): LapackInt {.importc: "LAPACKE_sgelq",
    dynlib: libName.}
proc sgelqWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
               lda: LapackInt; t: ptr cfloat; tsize: LapackInt; work: ptr cfloat;
               lwork: LapackInt): LapackInt {.importc: "LAPACKE_sgelq_work",
    dynlib: libName.}
proc sgemlq*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr cfloat; lda: LapackInt; t: ptr cfloat; tsize: LapackInt;
            c: ptr cfloat; ldc: LapackInt): LapackInt {.importc: "LAPACKE_sgemlq",
    dynlib: libName.}
proc sgemlqWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr cfloat; lda: LapackInt; t: ptr cfloat;
                tsize: LapackInt; c: ptr cfloat; ldc: LapackInt; work: ptr cfloat;
                lwork: LapackInt): LapackInt {.importc: "LAPACKE_sgemlq_work",
    dynlib: libName.}
proc zgelq*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr LapackComplexDouble;
           lda: LapackInt; t: ptr LapackComplexDouble; tsize: LapackInt): LapackInt {.
    importc: "LAPACKE_zgelq", dynlib: libName.}
proc zgelqWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
               a: ptr LapackComplexDouble; lda: LapackInt;
               t: ptr LapackComplexDouble; tsize: LapackInt;
               work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zgelq_work", dynlib: libName.}
proc zgemlq*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
            k: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
            t: ptr LapackComplexDouble; tsize: LapackInt; c: ptr LapackComplexDouble;
            ldc: LapackInt): LapackInt {.importc: "LAPACKE_zgemlq", dynlib: libName.}
proc zgemlqWork*(matrixLayout: cint; side: char; trans: char; m: LapackInt; n: LapackInt;
                k: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                t: ptr LapackComplexDouble; tsize: LapackInt;
                c: ptr LapackComplexDouble; ldc: LapackInt;
                work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zgemlq_work", dynlib: libName.}
proc cgetsls*(matrixLayout: cint; trans: char; m: LapackInt; n: LapackInt;
             nrhs: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
             b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_cgetsls", dynlib: libName.}
proc cgetslsWork*(matrixLayout: cint; trans: char; m: LapackInt; n: LapackInt;
                 nrhs: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                 b: ptr LapackComplexFloat; ldb: LapackInt;
                 work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cgetsls_work", dynlib: libName.}
proc dgetsls*(matrixLayout: cint; trans: char; m: LapackInt; n: LapackInt;
             nrhs: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble;
             ldb: LapackInt): LapackInt {.importc: "LAPACKE_dgetsls", dynlib: libName.}
proc dgetslsWork*(matrixLayout: cint; trans: char; m: LapackInt; n: LapackInt;
                 nrhs: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble;
                 ldb: LapackInt; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dgetsls_work", dynlib: libName.}
proc sgetsls*(matrixLayout: cint; trans: char; m: LapackInt; n: LapackInt;
             nrhs: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_sgetsls", dynlib: libName.}
proc sgetslsWork*(matrixLayout: cint; trans: char; m: LapackInt; n: LapackInt;
                 nrhs: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat;
                 ldb: LapackInt; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_sgetsls_work", dynlib: libName.}
proc zgetsls*(matrixLayout: cint; trans: char; m: LapackInt; n: LapackInt;
             nrhs: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
             b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zgetsls", dynlib: libName.}
proc zgetslsWork*(matrixLayout: cint; trans: char; m: LapackInt; n: LapackInt;
                 nrhs: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                 b: ptr LapackComplexDouble; ldb: LapackInt;
                 work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zgetsls_work", dynlib: libName.}
proc checon3*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
             lda: LapackInt; e: ptr LapackComplexFloat; ipiv: ptr LapackInt;
             anorm: cfloat; rcond: ptr cfloat): LapackInt {.
    importc: "LAPACKE_checon_3", dynlib: libName.}
proc checon3Work*(matrixLayout: cint; uplo: char; n: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt;
                 e: ptr LapackComplexFloat; ipiv: ptr LapackInt; anorm: cfloat;
                 rcond: ptr cfloat; work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_checon_3_work", dynlib: libName.}
proc chesvRk*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; e: ptr LapackComplexFloat;
             ipiv: ptr LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_chesv_rk", dynlib: libName.}
proc chesvRkWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt;
                 e: ptr LapackComplexFloat; ipiv: ptr LapackInt;
                 b: ptr LapackComplexFloat; ldb: LapackInt;
                 work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_chesv_rk_work", dynlib: libName.}
proc chetrfRk*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
              lda: LapackInt; e: ptr LapackComplexFloat; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_chetrf_rk", dynlib: libName.}
proc chetrfRkWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                  a: ptr LapackComplexFloat; lda: LapackInt;
                  e: ptr LapackComplexFloat; ipiv: ptr LapackInt;
                  work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_chetrf_rk_work", dynlib: libName.}
proc chetri3*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
             lda: LapackInt; e: ptr LapackComplexFloat; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_chetri_3", dynlib: libName.}
proc chetri3Work*(matrixLayout: cint; uplo: char; n: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt;
                 e: ptr LapackComplexFloat; ipiv: ptr LapackInt;
                 work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_chetri_3_work", dynlib: libName.}
proc chetrs3*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; e: ptr LapackComplexFloat;
             ipiv: ptr LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_chetrs_3", dynlib: libName.}
proc chetrs3Work*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt;
                 e: ptr LapackComplexFloat; ipiv: ptr LapackInt;
                 b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_chetrs_3_work", dynlib: libName.}
proc csycon3*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
             lda: LapackInt; e: ptr LapackComplexFloat; ipiv: ptr LapackInt;
             anorm: cfloat; rcond: ptr cfloat): LapackInt {.
    importc: "LAPACKE_csycon_3", dynlib: libName.}
proc csycon3Work*(matrixLayout: cint; uplo: char; n: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt;
                 e: ptr LapackComplexFloat; ipiv: ptr LapackInt; anorm: cfloat;
                 rcond: ptr cfloat; work: ptr LapackComplexFloat): LapackInt {.
    importc: "LAPACKE_csycon_3_work", dynlib: libName.}
proc csysvAa*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
             b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_csysv_aa", dynlib: libName.}
proc csysvAaWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                 b: ptr LapackComplexFloat; ldb: LapackInt;
                 work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_csysv_aa_work", dynlib: libName.}
proc csysvRk*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; e: ptr LapackComplexFloat;
             ipiv: ptr LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_csysv_rk", dynlib: libName.}
proc csysvRkWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt;
                 e: ptr LapackComplexFloat; ipiv: ptr LapackInt;
                 b: ptr LapackComplexFloat; ldb: LapackInt;
                 work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_csysv_rk_work", dynlib: libName.}
proc csytrfAa*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
              lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_csytrf_aa", dynlib: libName.}
proc csytrfAaWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                  a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                  work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_csytrf_aa_work", dynlib: libName.}
proc csytrfRk*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
              lda: LapackInt; e: ptr LapackComplexFloat; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_csytrf_rk", dynlib: libName.}
proc csytrfRkWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                  a: ptr LapackComplexFloat; lda: LapackInt;
                  e: ptr LapackComplexFloat; ipiv: ptr LapackInt;
                  work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_csytrf_rk_work", dynlib: libName.}
proc csytri3*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexFloat;
             lda: LapackInt; e: ptr LapackComplexFloat; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_csytri_3", dynlib: libName.}
proc csytri3Work*(matrixLayout: cint; uplo: char; n: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt;
                 e: ptr LapackComplexFloat; ipiv: ptr LapackInt;
                 work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_csytri_3_work", dynlib: libName.}
proc csytrs3*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexFloat; lda: LapackInt; e: ptr LapackComplexFloat;
             ipiv: ptr LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_csytrs_3", dynlib: libName.}
proc csytrs3Work*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt;
                 e: ptr LapackComplexFloat; ipiv: ptr LapackInt;
                 b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_csytrs_3_work", dynlib: libName.}
proc csytrsAa*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
              a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
              b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_csytrs_aa", dynlib: libName.}
proc csytrsAaWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                  a: ptr LapackComplexFloat; lda: LapackInt; ipiv: ptr LapackInt;
                  b: ptr LapackComplexFloat; ldb: LapackInt;
                  work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_csytrs_aa_work", dynlib: libName.}
proc dsycon3*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
             lda: LapackInt; e: ptr cdouble; ipiv: ptr LapackInt; anorm: cdouble;
             rcond: ptr cdouble): LapackInt {.importc: "LAPACKE_dsycon_3",
    dynlib: libName.}
proc dsycon3Work*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                 lda: LapackInt; e: ptr cdouble; ipiv: ptr LapackInt; anorm: cdouble;
                 rcond: ptr cdouble; work: ptr cdouble; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsycon_3_work", dynlib: libName.}
proc dsysvRk*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr cdouble; lda: LapackInt; e: ptr cdouble; ipiv: ptr LapackInt;
             b: ptr cdouble; ldb: LapackInt): LapackInt {.importc: "LAPACKE_dsysv_rk",
    dynlib: libName.}
proc dsysvRkWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr cdouble; lda: LapackInt; e: ptr cdouble; ipiv: ptr LapackInt;
                 b: ptr cdouble; ldb: LapackInt; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dsysv_rk_work", dynlib: libName.}
proc dsytrfRk*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
              lda: LapackInt; e: ptr cdouble; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsytrf_rk", dynlib: libName.}
proc dsytrfRkWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                  lda: LapackInt; e: ptr cdouble; ipiv: ptr LapackInt;
                  work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dsytrf_rk_work", dynlib: libName.}
proc dsytri3*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
             lda: LapackInt; e: ptr cdouble; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsytri_3", dynlib: libName.}
proc dsytri3Work*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                 lda: LapackInt; e: ptr cdouble; ipiv: ptr LapackInt; work: ptr cdouble;
                 lwork: LapackInt): LapackInt {.importc: "LAPACKE_dsytri_3_work",
    dynlib: libName.}
proc dsytrs3*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr cdouble; lda: LapackInt; e: ptr cdouble; ipiv: ptr LapackInt;
             b: ptr cdouble; ldb: LapackInt): LapackInt {.importc: "LAPACKE_dsytrs_3",
    dynlib: libName.}
proc dsytrs3Work*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr cdouble; lda: LapackInt; e: ptr cdouble; ipiv: ptr LapackInt;
                 b: ptr cdouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dsytrs_3_work", dynlib: libName.}
proc ssycon3*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat; lda: LapackInt;
             e: ptr cfloat; ipiv: ptr LapackInt; anorm: cfloat; rcond: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssycon_3", dynlib: libName.}
proc ssycon3Work*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                 lda: LapackInt; e: ptr cfloat; ipiv: ptr LapackInt; anorm: cfloat;
                 rcond: ptr cfloat; work: ptr cfloat; iwork: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ssycon_3_work", dynlib: libName.}
proc ssysvRk*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr cfloat; lda: LapackInt; e: ptr cfloat; ipiv: ptr LapackInt;
             b: ptr cfloat; ldb: LapackInt): LapackInt {.importc: "LAPACKE_ssysv_rk",
    dynlib: libName.}
proc ssysvRkWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr cfloat; lda: LapackInt; e: ptr cfloat; ipiv: ptr LapackInt;
                 b: ptr cfloat; ldb: LapackInt; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ssysv_rk_work", dynlib: libName.}
proc ssytrfRk*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
              lda: LapackInt; e: ptr cfloat; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ssytrf_rk", dynlib: libName.}
proc ssytrfRkWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                  lda: LapackInt; e: ptr cfloat; ipiv: ptr LapackInt; work: ptr cfloat;
                  lwork: LapackInt): LapackInt {.importc: "LAPACKE_ssytrf_rk_work",
    dynlib: libName.}
proc ssytri3*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat; lda: LapackInt;
             e: ptr cfloat; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ssytri_3", dynlib: libName.}
proc ssytri3Work*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                 lda: LapackInt; e: ptr cfloat; ipiv: ptr LapackInt; work: ptr cfloat;
                 lwork: LapackInt): LapackInt {.importc: "LAPACKE_ssytri_3_work",
    dynlib: libName.}
proc ssytrs3*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr cfloat; lda: LapackInt; e: ptr cfloat; ipiv: ptr LapackInt;
             b: ptr cfloat; ldb: LapackInt): LapackInt {.importc: "LAPACKE_ssytrs_3",
    dynlib: libName.}
proc ssytrs3Work*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr cfloat; lda: LapackInt; e: ptr cfloat; ipiv: ptr LapackInt;
                 b: ptr cfloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ssytrs_3_work", dynlib: libName.}
proc zhecon3*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexDouble;
             lda: LapackInt; e: ptr LapackComplexDouble; ipiv: ptr LapackInt;
             anorm: cdouble; rcond: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zhecon_3", dynlib: libName.}
proc zhecon3Work*(matrixLayout: cint; uplo: char; n: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt;
                 e: ptr LapackComplexDouble; ipiv: ptr LapackInt; anorm: cdouble;
                 rcond: ptr cdouble; work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zhecon_3_work", dynlib: libName.}
proc zhesvRk*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; e: ptr LapackComplexDouble;
             ipiv: ptr LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zhesv_rk", dynlib: libName.}
proc zhesvRkWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt;
                 e: ptr LapackComplexDouble; ipiv: ptr LapackInt;
                 b: ptr LapackComplexDouble; ldb: LapackInt;
                 work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zhesv_rk_work", dynlib: libName.}
proc zhetrfRk*(matrixLayout: cint; uplo: char; n: LapackInt;
              a: ptr LapackComplexDouble; lda: LapackInt; e: ptr LapackComplexDouble;
              ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_zhetrf_rk",
    dynlib: libName.}
proc zhetrfRkWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                  a: ptr LapackComplexDouble; lda: LapackInt;
                  e: ptr LapackComplexDouble; ipiv: ptr LapackInt;
                  work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zhetrf_rk_work", dynlib: libName.}
proc zhetri3*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexDouble;
             lda: LapackInt; e: ptr LapackComplexDouble; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zhetri_3", dynlib: libName.}
proc zhetri3Work*(matrixLayout: cint; uplo: char; n: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt;
                 e: ptr LapackComplexDouble; ipiv: ptr LapackInt;
                 work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zhetri_3_work", dynlib: libName.}
proc zhetrs3*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; e: ptr LapackComplexDouble;
             ipiv: ptr LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zhetrs_3", dynlib: libName.}
proc zhetrs3Work*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt;
                 e: ptr LapackComplexDouble; ipiv: ptr LapackInt;
                 b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zhetrs_3_work", dynlib: libName.}
proc zsycon3*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexDouble;
             lda: LapackInt; e: ptr LapackComplexDouble; ipiv: ptr LapackInt;
             anorm: cdouble; rcond: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zsycon_3", dynlib: libName.}
proc zsycon3Work*(matrixLayout: cint; uplo: char; n: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt;
                 e: ptr LapackComplexDouble; ipiv: ptr LapackInt; anorm: cdouble;
                 rcond: ptr cdouble; work: ptr LapackComplexDouble): LapackInt {.
    importc: "LAPACKE_zsycon_3_work", dynlib: libName.}
proc zsysvAa*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
             b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zsysv_aa", dynlib: libName.}
proc zsysvAaWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                 b: ptr LapackComplexDouble; ldb: LapackInt;
                 work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zsysv_aa_work", dynlib: libName.}
proc zsysvRk*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; e: ptr LapackComplexDouble;
             ipiv: ptr LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zsysv_rk", dynlib: libName.}
proc zsysvRkWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt;
                 e: ptr LapackComplexDouble; ipiv: ptr LapackInt;
                 b: ptr LapackComplexDouble; ldb: LapackInt;
                 work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zsysv_rk_work", dynlib: libName.}
proc zsytrfAa*(matrixLayout: cint; uplo: char; n: LapackInt;
              a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zsytrf_aa", dynlib: libName.}
proc zsytrfAaWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                  a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                  work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zsytrf_aa_work", dynlib: libName.}
proc zsytrfRk*(matrixLayout: cint; uplo: char; n: LapackInt;
              a: ptr LapackComplexDouble; lda: LapackInt; e: ptr LapackComplexDouble;
              ipiv: ptr LapackInt): LapackInt {.importc: "LAPACKE_zsytrf_rk",
    dynlib: libName.}
proc zsytrfRkWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                  a: ptr LapackComplexDouble; lda: LapackInt;
                  e: ptr LapackComplexDouble; ipiv: ptr LapackInt;
                  work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zsytrf_rk_work", dynlib: libName.}
proc zsytri3*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr LapackComplexDouble;
             lda: LapackInt; e: ptr LapackComplexDouble; ipiv: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zsytri_3", dynlib: libName.}
proc zsytri3Work*(matrixLayout: cint; uplo: char; n: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt;
                 e: ptr LapackComplexDouble; ipiv: ptr LapackInt;
                 work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zsytri_3_work", dynlib: libName.}
proc zsytrs3*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
             a: ptr LapackComplexDouble; lda: LapackInt; e: ptr LapackComplexDouble;
             ipiv: ptr LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zsytrs_3", dynlib: libName.}
proc zsytrs3Work*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt;
                 e: ptr LapackComplexDouble; ipiv: ptr LapackInt;
                 b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zsytrs_3_work", dynlib: libName.}
proc zsytrsAa*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
              a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
              b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zsytrs_aa", dynlib: libName.}
proc zsytrsAaWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                  a: ptr LapackComplexDouble; lda: LapackInt; ipiv: ptr LapackInt;
                  b: ptr LapackComplexDouble; ldb: LapackInt;
                  work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zsytrs_aa_work", dynlib: libName.}
proc chbev2stage*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                 kd: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
                 w: ptr cfloat; z: ptr LapackComplexFloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_chbev_2stage", dynlib: libName.}
proc chbev2stageWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                     kd: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
                     w: ptr cfloat; z: ptr LapackComplexFloat; ldz: LapackInt;
                     work: ptr LapackComplexFloat; lwork: LapackInt;
                     rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_chbev_2stage_work", dynlib: libName.}
proc chbevd2stage*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                  kd: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
                  w: ptr cfloat; z: ptr LapackComplexFloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_chbevd_2stage", dynlib: libName.}
proc chbevd2stageWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                      kd: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
                      w: ptr cfloat; z: ptr LapackComplexFloat; ldz: LapackInt;
                      work: ptr LapackComplexFloat; lwork: LapackInt;
                      rwork: ptr cfloat; lrwork: LapackInt; iwork: ptr LapackInt;
                      liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_chbevd_2stage_work", dynlib: libName.}
proc chbevx2stage*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                  kd: LapackInt; ab: ptr LapackComplexFloat; ldab: LapackInt;
                  q: ptr LapackComplexFloat; ldq: LapackInt; vl: cfloat; vu: cfloat;
                  il: LapackInt; iu: LapackInt; abstol: cfloat; m: ptr LapackInt;
                  w: ptr cfloat; z: ptr LapackComplexFloat; ldz: LapackInt;
                  ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_chbevx_2stage", dynlib: libName.}
proc chbevx2stageWork*(matrixLayout: cint; jobz: char; range: char; uplo: char;
                      n: LapackInt; kd: LapackInt; ab: ptr LapackComplexFloat;
                      ldab: LapackInt; q: ptr LapackComplexFloat; ldq: LapackInt;
                      vl: cfloat; vu: cfloat; il: LapackInt; iu: LapackInt;
                      abstol: cfloat; m: ptr LapackInt; w: ptr cfloat;
                      z: ptr LapackComplexFloat; ldz: LapackInt;
                      work: ptr LapackComplexFloat; lwork: LapackInt;
                      rwork: ptr cfloat; iwork: ptr LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_chbevx_2stage_work", dynlib: libName.}
proc cheev2stage*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt; w: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cheev_2stage", dynlib: libName.}
proc cheev2stageWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                     a: ptr LapackComplexFloat; lda: LapackInt; w: ptr cfloat;
                     work: ptr LapackComplexFloat; lwork: LapackInt;
                     rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cheev_2stage_work", dynlib: libName.}
proc cheevd2stage*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                  a: ptr LapackComplexFloat; lda: LapackInt; w: ptr cfloat): LapackInt {.
    importc: "LAPACKE_cheevd_2stage", dynlib: libName.}
proc cheevd2stageWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                      a: ptr LapackComplexFloat; lda: LapackInt; w: ptr cfloat;
                      work: ptr LapackComplexFloat; lwork: LapackInt;
                      rwork: ptr cfloat; lrwork: LapackInt; iwork: ptr LapackInt;
                      liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cheevd_2stage_work", dynlib: libName.}
proc cheevr2stage*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                  a: ptr LapackComplexFloat; lda: LapackInt; vl: cfloat; vu: cfloat;
                  il: LapackInt; iu: LapackInt; abstol: cfloat; m: ptr LapackInt;
                  w: ptr cfloat; z: ptr LapackComplexFloat; ldz: LapackInt;
                  isuppz: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_cheevr_2stage", dynlib: libName.}
proc cheevr2stageWork*(matrixLayout: cint; jobz: char; range: char; uplo: char;
                      n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                      vl: cfloat; vu: cfloat; il: LapackInt; iu: LapackInt;
                      abstol: cfloat; m: ptr LapackInt; w: ptr cfloat;
                      z: ptr LapackComplexFloat; ldz: LapackInt;
                      isuppz: ptr LapackInt; work: ptr LapackComplexFloat;
                      lwork: LapackInt; rwork: ptr cfloat; lrwork: LapackInt;
                      iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_cheevr_2stage_work", dynlib: libName.}
proc cheevx2stage*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                  a: ptr LapackComplexFloat; lda: LapackInt; vl: cfloat; vu: cfloat;
                  il: LapackInt; iu: LapackInt; abstol: cfloat; m: ptr LapackInt;
                  w: ptr cfloat; z: ptr LapackComplexFloat; ldz: LapackInt;
                  ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_cheevx_2stage", dynlib: libName.}
proc cheevx2stageWork*(matrixLayout: cint; jobz: char; range: char; uplo: char;
                      n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                      vl: cfloat; vu: cfloat; il: LapackInt; iu: LapackInt;
                      abstol: cfloat; m: ptr LapackInt; w: ptr cfloat;
                      z: ptr LapackComplexFloat; ldz: LapackInt;
                      work: ptr LapackComplexFloat; lwork: LapackInt;
                      rwork: ptr cfloat; iwork: ptr LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_cheevx_2stage_work", dynlib: libName.}
proc chegv2stage*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
                 n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                 b: ptr LapackComplexFloat; ldb: LapackInt; w: ptr cfloat): LapackInt {.
    importc: "LAPACKE_chegv_2stage", dynlib: libName.}
proc chegv2stageWork*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
                     n: LapackInt; a: ptr LapackComplexFloat; lda: LapackInt;
                     b: ptr LapackComplexFloat; ldb: LapackInt; w: ptr cfloat;
                     work: ptr LapackComplexFloat; lwork: LapackInt;
                     rwork: ptr cfloat): LapackInt {.
    importc: "LAPACKE_chegv_2stage_work", dynlib: libName.}
proc dsbev2stage*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                 kd: LapackInt; ab: ptr cdouble; ldab: LapackInt; w: ptr cdouble;
                 z: ptr cdouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_dsbev_2stage", dynlib: libName.}
proc dsbev2stageWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                     kd: LapackInt; ab: ptr cdouble; ldab: LapackInt; w: ptr cdouble;
                     z: ptr cdouble; ldz: LapackInt; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dsbev_2stage_work", dynlib: libName.}
proc dsbevd2stage*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                  kd: LapackInt; ab: ptr cdouble; ldab: LapackInt; w: ptr cdouble;
                  z: ptr cdouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_dsbevd_2stage", dynlib: libName.}
proc dsbevd2stageWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                      kd: LapackInt; ab: ptr cdouble; ldab: LapackInt; w: ptr cdouble;
                      z: ptr cdouble; ldz: LapackInt; work: ptr cdouble;
                      lwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dsbevd_2stage_work", dynlib: libName.}
proc dsbevx2stage*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                  kd: LapackInt; ab: ptr cdouble; ldab: LapackInt; q: ptr cdouble;
                  ldq: LapackInt; vl: cdouble; vu: cdouble; il: LapackInt;
                  iu: LapackInt; abstol: cdouble; m: ptr LapackInt; w: ptr cdouble;
                  z: ptr cdouble; ldz: LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsbevx_2stage", dynlib: libName.}
proc dsbevx2stageWork*(matrixLayout: cint; jobz: char; range: char; uplo: char;
                      n: LapackInt; kd: LapackInt; ab: ptr cdouble; ldab: LapackInt;
                      q: ptr cdouble; ldq: LapackInt; vl: cdouble; vu: cdouble;
                      il: LapackInt; iu: LapackInt; abstol: cdouble; m: ptr LapackInt;
                      w: ptr cdouble; z: ptr cdouble; ldz: LapackInt; work: ptr cdouble;
                      lwork: LapackInt; iwork: ptr LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsbevx_2stage_work", dynlib: libName.}
proc dsyev2stage*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                 a: ptr cdouble; lda: LapackInt; w: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsyev_2stage", dynlib: libName.}
proc dsyev2stageWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                     a: ptr cdouble; lda: LapackInt; w: ptr cdouble; work: ptr cdouble;
                     lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dsyev_2stage_work", dynlib: libName.}
proc dsyevd2stage*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                  a: ptr cdouble; lda: LapackInt; w: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsyevd_2stage", dynlib: libName.}
proc dsyevd2stageWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                      a: ptr cdouble; lda: LapackInt; w: ptr cdouble; work: ptr cdouble;
                      lwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dsyevd_2stage_work", dynlib: libName.}
proc dsyevr2stage*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                  a: ptr cdouble; lda: LapackInt; vl: cdouble; vu: cdouble;
                  il: LapackInt; iu: LapackInt; abstol: cdouble; m: ptr LapackInt;
                  w: ptr cdouble; z: ptr cdouble; ldz: LapackInt; isuppz: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsyevr_2stage", dynlib: libName.}
proc dsyevr2stageWork*(matrixLayout: cint; jobz: char; range: char; uplo: char;
                      n: LapackInt; a: ptr cdouble; lda: LapackInt; vl: cdouble;
                      vu: cdouble; il: LapackInt; iu: LapackInt; abstol: cdouble;
                      m: ptr LapackInt; w: ptr cdouble; z: ptr cdouble; ldz: LapackInt;
                      isuppz: ptr LapackInt; work: ptr cdouble; lwork: LapackInt;
                      iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dsyevr_2stage_work", dynlib: libName.}
proc dsyevx2stage*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                  a: ptr cdouble; lda: LapackInt; vl: cdouble; vu: cdouble;
                  il: LapackInt; iu: LapackInt; abstol: cdouble; m: ptr LapackInt;
                  w: ptr cdouble; z: ptr cdouble; ldz: LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsyevx_2stage", dynlib: libName.}
proc dsyevx2stageWork*(matrixLayout: cint; jobz: char; range: char; uplo: char;
                      n: LapackInt; a: ptr cdouble; lda: LapackInt; vl: cdouble;
                      vu: cdouble; il: LapackInt; iu: LapackInt; abstol: cdouble;
                      m: ptr LapackInt; w: ptr cdouble; z: ptr cdouble; ldz: LapackInt;
                      work: ptr cdouble; lwork: LapackInt; iwork: ptr LapackInt;
                      ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsyevx_2stage_work", dynlib: libName.}
proc dsygv2stage*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
                 n: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble;
                 ldb: LapackInt; w: ptr cdouble): LapackInt {.
    importc: "LAPACKE_dsygv_2stage", dynlib: libName.}
proc dsygv2stageWork*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
                     n: LapackInt; a: ptr cdouble; lda: LapackInt; b: ptr cdouble;
                     ldb: LapackInt; w: ptr cdouble; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dsygv_2stage_work", dynlib: libName.}
proc ssbev2stage*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                 kd: LapackInt; ab: ptr cfloat; ldab: LapackInt; w: ptr cfloat;
                 z: ptr cfloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_ssbev_2stage", dynlib: libName.}
proc ssbev2stageWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                     kd: LapackInt; ab: ptr cfloat; ldab: LapackInt; w: ptr cfloat;
                     z: ptr cfloat; ldz: LapackInt; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ssbev_2stage_work", dynlib: libName.}
proc ssbevd2stage*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                  kd: LapackInt; ab: ptr cfloat; ldab: LapackInt; w: ptr cfloat;
                  z: ptr cfloat; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_ssbevd_2stage", dynlib: libName.}
proc ssbevd2stageWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                      kd: LapackInt; ab: ptr cfloat; ldab: LapackInt; w: ptr cfloat;
                      z: ptr cfloat; ldz: LapackInt; work: ptr cfloat; lwork: LapackInt;
                      iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ssbevd_2stage_work", dynlib: libName.}
proc ssbevx2stage*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                  kd: LapackInt; ab: ptr cfloat; ldab: LapackInt; q: ptr cfloat;
                  ldq: LapackInt; vl: cfloat; vu: cfloat; il: LapackInt; iu: LapackInt;
                  abstol: cfloat; m: ptr LapackInt; w: ptr cfloat; z: ptr cfloat;
                  ldz: LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ssbevx_2stage", dynlib: libName.}
proc ssbevx2stageWork*(matrixLayout: cint; jobz: char; range: char; uplo: char;
                      n: LapackInt; kd: LapackInt; ab: ptr cfloat; ldab: LapackInt;
                      q: ptr cfloat; ldq: LapackInt; vl: cfloat; vu: cfloat;
                      il: LapackInt; iu: LapackInt; abstol: cfloat; m: ptr LapackInt;
                      w: ptr cfloat; z: ptr cfloat; ldz: LapackInt; work: ptr cfloat;
                      lwork: LapackInt; iwork: ptr LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ssbevx_2stage_work", dynlib: libName.}
proc ssyev2stage*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                 a: ptr cfloat; lda: LapackInt; w: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssyev_2stage", dynlib: libName.}
proc ssyev2stageWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                     a: ptr cfloat; lda: LapackInt; w: ptr cfloat; work: ptr cfloat;
                     lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ssyev_2stage_work", dynlib: libName.}
proc ssyevd2stage*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                  a: ptr cfloat; lda: LapackInt; w: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssyevd_2stage", dynlib: libName.}
proc ssyevd2stageWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                      a: ptr cfloat; lda: LapackInt; w: ptr cfloat; work: ptr cfloat;
                      lwork: LapackInt; iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ssyevd_2stage_work", dynlib: libName.}
proc ssyevr2stage*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                  a: ptr cfloat; lda: LapackInt; vl: cfloat; vu: cfloat; il: LapackInt;
                  iu: LapackInt; abstol: cfloat; m: ptr LapackInt; w: ptr cfloat;
                  z: ptr cfloat; ldz: LapackInt; isuppz: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ssyevr_2stage", dynlib: libName.}
proc ssyevr2stageWork*(matrixLayout: cint; jobz: char; range: char; uplo: char;
                      n: LapackInt; a: ptr cfloat; lda: LapackInt; vl: cfloat;
                      vu: cfloat; il: LapackInt; iu: LapackInt; abstol: cfloat;
                      m: ptr LapackInt; w: ptr cfloat; z: ptr cfloat; ldz: LapackInt;
                      isuppz: ptr LapackInt; work: ptr cfloat; lwork: LapackInt;
                      iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ssyevr_2stage_work", dynlib: libName.}
proc ssyevx2stage*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                  a: ptr cfloat; lda: LapackInt; vl: cfloat; vu: cfloat; il: LapackInt;
                  iu: LapackInt; abstol: cfloat; m: ptr LapackInt; w: ptr cfloat;
                  z: ptr cfloat; ldz: LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ssyevx_2stage", dynlib: libName.}
proc ssyevx2stageWork*(matrixLayout: cint; jobz: char; range: char; uplo: char;
                      n: LapackInt; a: ptr cfloat; lda: LapackInt; vl: cfloat;
                      vu: cfloat; il: LapackInt; iu: LapackInt; abstol: cfloat;
                      m: ptr LapackInt; w: ptr cfloat; z: ptr cfloat; ldz: LapackInt;
                      work: ptr cfloat; lwork: LapackInt; iwork: ptr LapackInt;
                      ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ssyevx_2stage_work", dynlib: libName.}
proc ssygv2stage*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
                 n: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat;
                 ldb: LapackInt; w: ptr cfloat): LapackInt {.
    importc: "LAPACKE_ssygv_2stage", dynlib: libName.}
proc ssygv2stageWork*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
                     n: LapackInt; a: ptr cfloat; lda: LapackInt; b: ptr cfloat;
                     ldb: LapackInt; w: ptr cfloat; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ssygv_2stage_work", dynlib: libName.}
proc zhbev2stage*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                 kd: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
                 w: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_zhbev_2stage", dynlib: libName.}
proc zhbev2stageWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                     kd: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
                     w: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt;
                     work: ptr LapackComplexDouble; lwork: LapackInt;
                     rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zhbev_2stage_work", dynlib: libName.}
proc zhbevd2stage*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                  kd: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
                  w: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt): LapackInt {.
    importc: "LAPACKE_zhbevd_2stage", dynlib: libName.}
proc zhbevd2stageWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                      kd: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
                      w: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt;
                      work: ptr LapackComplexDouble; lwork: LapackInt;
                      rwork: ptr cdouble; lrwork: LapackInt; iwork: ptr LapackInt;
                      liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zhbevd_2stage_work", dynlib: libName.}
proc zhbevx2stage*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                  kd: LapackInt; ab: ptr LapackComplexDouble; ldab: LapackInt;
                  q: ptr LapackComplexDouble; ldq: LapackInt; vl: cdouble; vu: cdouble;
                  il: LapackInt; iu: LapackInt; abstol: cdouble; m: ptr LapackInt;
                  w: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt;
                  ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zhbevx_2stage", dynlib: libName.}
proc zhbevx2stageWork*(matrixLayout: cint; jobz: char; range: char; uplo: char;
                      n: LapackInt; kd: LapackInt; ab: ptr LapackComplexDouble;
                      ldab: LapackInt; q: ptr LapackComplexDouble; ldq: LapackInt;
                      vl: cdouble; vu: cdouble; il: LapackInt; iu: LapackInt;
                      abstol: cdouble; m: ptr LapackInt; w: ptr cdouble;
                      z: ptr LapackComplexDouble; ldz: LapackInt;
                      work: ptr LapackComplexDouble; lwork: LapackInt;
                      rwork: ptr cdouble; iwork: ptr LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zhbevx_2stage_work", dynlib: libName.}
proc zheev2stage*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt; w: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zheev_2stage", dynlib: libName.}
proc zheev2stageWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                     a: ptr LapackComplexDouble; lda: LapackInt; w: ptr cdouble;
                     work: ptr LapackComplexDouble; lwork: LapackInt;
                     rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zheev_2stage_work", dynlib: libName.}
proc zheevd2stage*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                  a: ptr LapackComplexDouble; lda: LapackInt; w: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zheevd_2stage", dynlib: libName.}
proc zheevd2stageWork*(matrixLayout: cint; jobz: char; uplo: char; n: LapackInt;
                      a: ptr LapackComplexDouble; lda: LapackInt; w: ptr cdouble;
                      work: ptr LapackComplexDouble; lwork: LapackInt;
                      rwork: ptr cdouble; lrwork: LapackInt; iwork: ptr LapackInt;
                      liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zheevd_2stage_work", dynlib: libName.}
proc zheevr2stage*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                  a: ptr LapackComplexDouble; lda: LapackInt; vl: cdouble; vu: cdouble;
                  il: LapackInt; iu: LapackInt; abstol: cdouble; m: ptr LapackInt;
                  w: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt;
                  isuppz: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zheevr_2stage", dynlib: libName.}
proc zheevr2stageWork*(matrixLayout: cint; jobz: char; range: char; uplo: char;
                      n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                      vl: cdouble; vu: cdouble; il: LapackInt; iu: LapackInt;
                      abstol: cdouble; m: ptr LapackInt; w: ptr cdouble;
                      z: ptr LapackComplexDouble; ldz: LapackInt;
                      isuppz: ptr LapackInt; work: ptr LapackComplexDouble;
                      lwork: LapackInt; rwork: ptr cdouble; lrwork: LapackInt;
                      iwork: ptr LapackInt; liwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zheevr_2stage_work", dynlib: libName.}
proc zheevx2stage*(matrixLayout: cint; jobz: char; range: char; uplo: char; n: LapackInt;
                  a: ptr LapackComplexDouble; lda: LapackInt; vl: cdouble; vu: cdouble;
                  il: LapackInt; iu: LapackInt; abstol: cdouble; m: ptr LapackInt;
                  w: ptr cdouble; z: ptr LapackComplexDouble; ldz: LapackInt;
                  ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zheevx_2stage", dynlib: libName.}
proc zheevx2stageWork*(matrixLayout: cint; jobz: char; range: char; uplo: char;
                      n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                      vl: cdouble; vu: cdouble; il: LapackInt; iu: LapackInt;
                      abstol: cdouble; m: ptr LapackInt; w: ptr cdouble;
                      z: ptr LapackComplexDouble; ldz: LapackInt;
                      work: ptr LapackComplexDouble; lwork: LapackInt;
                      rwork: ptr cdouble; iwork: ptr LapackInt; ifail: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zheevx_2stage_work", dynlib: libName.}
proc zhegv2stage*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
                 n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                 b: ptr LapackComplexDouble; ldb: LapackInt; w: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zhegv_2stage", dynlib: libName.}
proc zhegv2stageWork*(matrixLayout: cint; itype: LapackInt; jobz: char; uplo: char;
                     n: LapackInt; a: ptr LapackComplexDouble; lda: LapackInt;
                     b: ptr LapackComplexDouble; ldb: LapackInt; w: ptr cdouble;
                     work: ptr LapackComplexDouble; lwork: LapackInt;
                     rwork: ptr cdouble): LapackInt {.
    importc: "LAPACKE_zhegv_2stage_work", dynlib: libName.}
proc mklCgetrfnp*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                 a: ptr LapackComplexFloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_cgetrfnp", dynlib: libName.}
proc mklCgetrfnpWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                     a: ptr LapackComplexFloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_cgetrfnp_work", dynlib: libName.}
proc mklDgetrfnp*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
                 lda: LapackInt): LapackInt {.importc: "LAPACKE_mkl_dgetrfnp",
    dynlib: libName.}
proc mklDgetrfnpWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cdouble;
                     lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_dgetrfnp_work", dynlib: libName.}
proc mklSgetrfnp*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
                 lda: LapackInt): LapackInt {.importc: "LAPACKE_mkl_sgetrfnp",
    dynlib: libName.}
proc mklSgetrfnpWork*(matrixLayout: cint; m: LapackInt; n: LapackInt; a: ptr cfloat;
                     lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_sgetrfnp_work", dynlib: libName.}
proc mklZgetrfnp*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                 a: ptr LapackComplexDouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_zgetrfnp", dynlib: libName.}
proc mklZgetrfnpWork*(matrixLayout: cint; m: LapackInt; n: LapackInt;
                     a: ptr LapackComplexDouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_zgetrfnp_work", dynlib: libName.}
proc mklCgetrinp*(matrixLayout: cint; n: LapackInt; a: ptr LapackComplexFloat;
                 lda: LapackInt): LapackInt {.importc: "LAPACKE_mkl_cgetrinp",
    dynlib: libName.}
proc mklCgetrinpWork*(matrixLayout: cint; n: LapackInt; a: ptr LapackComplexFloat;
                     lda: LapackInt; work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_cgetrinp_work", dynlib: libName.}
proc mklDgetrinp*(matrixLayout: cint; n: LapackInt; a: ptr cdouble; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_dgetrinp", dynlib: libName.}
proc mklDgetrinpWork*(matrixLayout: cint; n: LapackInt; a: ptr cdouble; lda: LapackInt;
                     work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_dgetrinp_work", dynlib: libName.}
proc mklSgetrinp*(matrixLayout: cint; n: LapackInt; a: ptr cfloat; lda: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_sgetrinp", dynlib: libName.}
proc mklSgetrinpWork*(matrixLayout: cint; n: LapackInt; a: ptr cfloat; lda: LapackInt;
                     work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_sgetrinp_work", dynlib: libName.}
proc mklZgetrinp*(matrixLayout: cint; n: LapackInt; a: ptr LapackComplexDouble;
                 lda: LapackInt): LapackInt {.importc: "LAPACKE_mkl_zgetrinp",
    dynlib: libName.}
proc mklZgetrinpWork*(matrixLayout: cint; n: LapackInt; a: ptr LapackComplexDouble;
                     lda: LapackInt; work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_mkl_zgetrinp_work", dynlib: libName.}
proc setNancheck*(flag: cint) {.importc: "LAPACKE_set_nancheck", dynlib: libName.}
proc getNancheck*(): cint {.importc: "LAPACKE_get_nancheck", dynlib: libName.}
proc chesvAa2stage*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                   a: ptr LapackComplexFloat; lda: LapackInt;
                   tb: ptr LapackComplexFloat; ltb: LapackInt; ipiv: ptr LapackInt;
                   ipiv2: ptr LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_chesv_aa_2stage", dynlib: libName.}
proc chetrfAa2stage*(matrixLayout: cint; uplo: char; n: LapackInt;
                    a: ptr LapackComplexFloat; lda: LapackInt;
                    tb: ptr LapackComplexFloat; ltb: LapackInt; ipiv: ptr LapackInt;
                    ipiv2: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_chetrf_aa_2stage", dynlib: libName.}
proc chesvAa2stageWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                       a: ptr LapackComplexFloat; lda: LapackInt;
                       tb: ptr LapackComplexFloat; ltb: LapackInt;
                       ipiv: ptr LapackInt; ipiv2: ptr LapackInt;
                       b: ptr LapackComplexFloat; ldb: LapackInt;
                       work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_chesv_aa_2stage_work", dynlib: libName.}
proc chetrfAa2stageWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                        a: ptr LapackComplexFloat; lda: LapackInt;
                        tb: ptr LapackComplexFloat; ltb: LapackInt;
                        ipiv: ptr LapackInt; ipiv2: ptr LapackInt;
                        work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_chetrf_aa_2stage_work", dynlib: libName.}
proc chetrsAa2stage*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                    a: ptr LapackComplexFloat; lda: LapackInt;
                    tb: ptr LapackComplexFloat; ltb: LapackInt; ipiv: ptr LapackInt;
                    ipiv2: ptr LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_chetrs_aa_2stage", dynlib: libName.}
proc chetrsAa2stageWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                        a: ptr LapackComplexFloat; lda: LapackInt;
                        tb: ptr LapackComplexFloat; ltb: LapackInt;
                        ipiv: ptr LapackInt; ipiv2: ptr LapackInt;
                        b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_chetrs_aa_2stage_work", dynlib: libName.}
proc csysvAa2stage*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                   a: ptr LapackComplexFloat; lda: LapackInt;
                   tb: ptr LapackComplexFloat; ltb: LapackInt; ipiv: ptr LapackInt;
                   ipiv2: ptr LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_csysv_aa_2stage", dynlib: libName.}
proc csysvAa2stageWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                       a: ptr LapackComplexFloat; lda: LapackInt;
                       tb: ptr LapackComplexFloat; ltb: LapackInt;
                       ipiv: ptr LapackInt; ipiv2: ptr LapackInt;
                       b: ptr LapackComplexFloat; ldb: LapackInt;
                       work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_csysv_aa_2stage_work", dynlib: libName.}
proc csytrfAa2stage*(matrixLayout: cint; uplo: char; n: LapackInt;
                    a: ptr LapackComplexFloat; lda: LapackInt;
                    tb: ptr LapackComplexFloat; ltb: LapackInt; ipiv: ptr LapackInt;
                    ipiv2: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_csytrf_aa_2stage", dynlib: libName.}
proc csytrfAa2stageWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                        a: ptr LapackComplexFloat; lda: LapackInt;
                        tb: ptr LapackComplexFloat; ltb: LapackInt;
                        ipiv: ptr LapackInt; ipiv2: ptr LapackInt;
                        work: ptr LapackComplexFloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_csytrf_aa_2stage_work", dynlib: libName.}
proc csytrsAa2stage*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                    a: ptr LapackComplexFloat; lda: LapackInt;
                    tb: ptr LapackComplexFloat; ltb: LapackInt; ipiv: ptr LapackInt;
                    ipiv2: ptr LapackInt; b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_csytrs_aa_2stage", dynlib: libName.}
proc csytrsAa2stageWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                        a: ptr LapackComplexFloat; lda: LapackInt;
                        tb: ptr LapackComplexFloat; ltb: LapackInt;
                        ipiv: ptr LapackInt; ipiv2: ptr LapackInt;
                        b: ptr LapackComplexFloat; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_csytrs_aa_2stage_work", dynlib: libName.}
proc dsysvAa2stage*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                   a: ptr cdouble; lda: LapackInt; tb: ptr cdouble; ltb: LapackInt;
                   ipiv: ptr LapackInt; ipiv2: ptr LapackInt; b: ptr cdouble;
                   ldb: LapackInt): LapackInt {.importc: "LAPACKE_dsysv_aa_2stage",
    dynlib: libName.}
proc dsysvAa2stageWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                       a: ptr cdouble; lda: LapackInt; tb: ptr cdouble; ltb: LapackInt;
                       ipiv: ptr LapackInt; ipiv2: ptr LapackInt; b: ptr cdouble;
                       ldb: LapackInt; work: ptr cdouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dsysv_aa_2stage_work", dynlib: libName.}
proc dsytrfAa2stage*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                    lda: LapackInt; tb: ptr cdouble; ltb: LapackInt;
                    ipiv: ptr LapackInt; ipiv2: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_dsytrf_aa_2stage", dynlib: libName.}
proc dsytrfAa2stageWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cdouble;
                        lda: LapackInt; tb: ptr cdouble; ltb: LapackInt;
                        ipiv: ptr LapackInt; ipiv2: ptr LapackInt; work: ptr cdouble;
                        lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_dsytrf_aa_2stage_work", dynlib: libName.}
proc dsytrsAa2stage*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                    a: ptr cdouble; lda: LapackInt; tb: ptr cdouble; ltb: LapackInt;
                    ipiv: ptr LapackInt; ipiv2: ptr LapackInt; b: ptr cdouble;
                    ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dsytrs_aa_2stage", dynlib: libName.}
proc dsytrsAa2stageWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                        a: ptr cdouble; lda: LapackInt; tb: ptr cdouble; ltb: LapackInt;
                        ipiv: ptr LapackInt; ipiv2: ptr LapackInt; b: ptr cdouble;
                        ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_dsytrs_aa_2stage_work", dynlib: libName.}
proc zsysvAa2stage*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                   a: ptr LapackComplexDouble; lda: LapackInt;
                   tb: ptr LapackComplexDouble; ltb: LapackInt; ipiv: ptr LapackInt;
                   ipiv2: ptr LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zsysv_aa_2stage", dynlib: libName.}
proc zsysvAa2stageWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                       a: ptr LapackComplexDouble; lda: LapackInt;
                       tb: ptr LapackComplexDouble; ltb: LapackInt;
                       ipiv: ptr LapackInt; ipiv2: ptr LapackInt;
                       b: ptr LapackComplexDouble; ldb: LapackInt;
                       work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zsysv_aa_2stage_work", dynlib: libName.}
proc zsytrfAa2stage*(matrixLayout: cint; uplo: char; n: LapackInt;
                    a: ptr LapackComplexDouble; lda: LapackInt;
                    tb: ptr LapackComplexDouble; ltb: LapackInt; ipiv: ptr LapackInt;
                    ipiv2: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zsytrf_aa_2stage", dynlib: libName.}
proc zsytrfAa2stageWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                        a: ptr LapackComplexDouble; lda: LapackInt;
                        tb: ptr LapackComplexDouble; ltb: LapackInt;
                        ipiv: ptr LapackInt; ipiv2: ptr LapackInt;
                        work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zsytrf_aa_2stage_work", dynlib: libName.}
proc zsytrsAa2stage*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                    a: ptr LapackComplexDouble; lda: LapackInt;
                    tb: ptr LapackComplexDouble; ltb: LapackInt; ipiv: ptr LapackInt;
                    ipiv2: ptr LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zsytrs_aa_2stage", dynlib: libName.}
proc zsytrsAa2stageWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                        a: ptr LapackComplexDouble; lda: LapackInt;
                        tb: ptr LapackComplexDouble; ltb: LapackInt;
                        ipiv: ptr LapackInt; ipiv2: ptr LapackInt;
                        b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zsytrs_aa_2stage_work", dynlib: libName.}
proc ssysvAa2stage*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                   a: ptr cfloat; lda: LapackInt; tb: ptr cfloat; ltb: LapackInt;
                   ipiv: ptr LapackInt; ipiv2: ptr LapackInt; b: ptr cfloat;
                   ldb: LapackInt): LapackInt {.importc: "LAPACKE_ssysv_aa_2stage",
    dynlib: libName.}
proc ssysvAa2stageWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                       a: ptr cfloat; lda: LapackInt; tb: ptr cfloat; ltb: LapackInt;
                       ipiv: ptr LapackInt; ipiv2: ptr LapackInt; b: ptr cfloat;
                       ldb: LapackInt; work: ptr cfloat; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ssysv_aa_2stage_work", dynlib: libName.}
proc ssytrfAa2stage*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                    lda: LapackInt; tb: ptr cfloat; ltb: LapackInt;
                    ipiv: ptr LapackInt; ipiv2: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_ssytrf_aa_2stage", dynlib: libName.}
proc ssytrfAa2stageWork*(matrixLayout: cint; uplo: char; n: LapackInt; a: ptr cfloat;
                        lda: LapackInt; tb: ptr cfloat; ltb: LapackInt;
                        ipiv: ptr LapackInt; ipiv2: ptr LapackInt; work: ptr cfloat;
                        lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_ssytrf_aa_2stage_work", dynlib: libName.}
proc ssytrsAa2stage*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                    a: ptr cfloat; lda: LapackInt; tb: ptr cfloat; ltb: LapackInt;
                    ipiv: ptr LapackInt; ipiv2: ptr LapackInt; b: ptr cfloat;
                    ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ssytrs_aa_2stage", dynlib: libName.}
proc ssytrsAa2stageWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                        a: ptr cfloat; lda: LapackInt; tb: ptr cfloat; ltb: LapackInt;
                        ipiv: ptr LapackInt; ipiv2: ptr LapackInt; b: ptr cfloat;
                        ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_ssytrs_aa_2stage_work", dynlib: libName.}
proc zhesvAa2stage*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                   a: ptr LapackComplexDouble; lda: LapackInt;
                   tb: ptr LapackComplexDouble; ltb: LapackInt; ipiv: ptr LapackInt;
                   ipiv2: ptr LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zhesv_aa_2stage", dynlib: libName.}
proc zhesvAa2stageWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                       a: ptr LapackComplexDouble; lda: LapackInt;
                       tb: ptr LapackComplexDouble; ltb: LapackInt;
                       ipiv: ptr LapackInt; ipiv2: ptr LapackInt;
                       b: ptr LapackComplexDouble; ldb: LapackInt;
                       work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zhesv_aa_2stage_work", dynlib: libName.}
proc zhetrfAa2stageWork*(matrixLayout: cint; uplo: char; n: LapackInt;
                        a: ptr LapackComplexDouble; lda: LapackInt;
                        tb: ptr LapackComplexDouble; ltb: LapackInt;
                        ipiv: ptr LapackInt; ipiv2: ptr LapackInt;
                        work: ptr LapackComplexDouble; lwork: LapackInt): LapackInt {.
    importc: "LAPACKE_zhetrf_aa_2stage_work", dynlib: libName.}
proc zhetrfAa2stage*(matrixLayout: cint; uplo: char; n: LapackInt;
                    a: ptr LapackComplexDouble; lda: LapackInt;
                    tb: ptr LapackComplexDouble; ltb: LapackInt; ipiv: ptr LapackInt;
                    ipiv2: ptr LapackInt): LapackInt {.
    importc: "LAPACKE_zhetrf_aa_2stage", dynlib: libName.}
proc zhetrsAa2stage*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                    a: ptr LapackComplexDouble; lda: LapackInt;
                    tb: ptr LapackComplexDouble; ltb: LapackInt; ipiv: ptr LapackInt;
                    ipiv2: ptr LapackInt; b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zhetrs_aa_2stage", dynlib: libName.}
proc zhetrsAa2stageWork*(matrixLayout: cint; uplo: char; n: LapackInt; nrhs: LapackInt;
                        a: ptr LapackComplexDouble; lda: LapackInt;
                        tb: ptr LapackComplexDouble; ltb: LapackInt;
                        ipiv: ptr LapackInt; ipiv2: ptr LapackInt;
                        b: ptr LapackComplexDouble; ldb: LapackInt): LapackInt {.
    importc: "LAPACKE_zhetrs_aa_2stage_work", dynlib: libName.}

