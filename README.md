# timstof
## Set up
Update the library paths in src/timsdata.py so that it points to the binary files in lib

    if sys.platform[:5] == "win32":
      libname = "[put path here]\timsdata.dll"
    elif sys.platform[:5] == "linux":
      libname = "[put path here]/libtimsdata.so"
Packages required:
   * numpy
   * matplotlib
Can be installed with
    pip install numpy matplotlib
Python Version = 2.7

## To Run:
    python extract_msn_nopd.py input/ output/
