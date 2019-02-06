# timstof
## Set up
Update the library paths in timsdata.py so that it points to the binary files in lib

    if sys.platform[:5] == "win32":
      libname = "[put path here]\timsdata.dll"
    elif sys.platform[:5] == "linux":
      libname = "[put path here]/libtimsdata.so"
