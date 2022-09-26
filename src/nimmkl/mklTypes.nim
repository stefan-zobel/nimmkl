when defined(windows):
  const
    libSuffix = ".dll"
    libPrefix = ""
elif defined(macosx):
  const
    libSuffix = ".dylib"
    libPrefix = "lib"
else:
  const
    libSuffix = ".so(||.3|.2|.1|.0)"
    libPrefix = "lib"

const
  mkl {.strdefine.} = "mkl_rt"
  libName* = libPrefix & mkl & libSuffix

{.hint: "Using MKL library with name: " & libName .}


type
  Mkl_Int* = cint
  Mkl_Int8* = cchar
  Mkl_Int16* = cshort
  Mkl_Int32* = cint


type
  MKL_Complex8* {.importc: "MKL_Complex8", header: "mkl_types.h", bycopy.} = object
    real* {.importc: "real".}: cfloat
    imag* {.importc: "imag".}: cfloat

type
  MKL_Complex16* {.importc: "MKL_Complex16", header: "mkl_types.h", bycopy.} = object
    real* {.importc: "real".}: cdouble
    imag* {.importc: "imag".}: cdouble


type
  Cblas_Layout* {.size: sizeof(cint).} = enum
    CblasRowMajor = 101,
    CblasColMajor = 102
  Cblas_Transpose* {.size: sizeof(cint).} = enum
    CblasNoTrans = 111,
    CblasTrans = 112,
    CblasConjTrans = 113
  Cblas_Uplo* {.size: sizeof(cint).} = enum
    CblasUpper = 121,
    CblasLower = 122
  Cblas_Diag* {.size: sizeof(cint).} = enum
    CblasNonUnit = 131,
    CblasUnit = 132
  Cblas_Side* {.size: sizeof(cint).} = enum
    CblasLeft = 141,
    CblasRight = 142
  Cblas_Identifier* {.size: sizeof(cint).} = enum
    CblasAMatrix = 161,
    CblasBMatrix = 162
  Cblas_Offset* {.size: sizeof(cint).} = enum
    CblasRowOffset = 171,
    CblasColOffset = 172,
    CblasFixOffset = 173

type
  Cblas_Index* = csize_t


type
  LapackInt* = Mkl_Int
  LapackLogical* = LapackInt
  LapackComplexFloat* = MKL_Complex8
  LapackComplexDouble* = MKL_Complex16

