import mklTypes

type
  XerblaEntry* = proc (name: cstring; num: ptr cint; len: cint) {.cdecl.}

proc mklSetXerbla*(xerbla: XerblaEntry): XerblaEntry {.cdecl,
    importc: "mkl_set_xerbla", dynlib: libName.}
