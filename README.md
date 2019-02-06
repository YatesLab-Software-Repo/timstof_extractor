# timstof
## Set up
Python Version = 2.7

Packages required:
   * numpy
   * matplotlib
   
To Install Packages Run:

    pip install numpy matplotlib
Update the library paths in src/timsdata.py so that it points to the binary files in lib

    if sys.platform[:5] == "win32":
      libname = "[put path here]\timsdata.dll"
    elif sys.platform[:5] == "linux":
      libname = "[put path here]/libtimsdata.so"



## To Run:
    python extract_msn_nopd.py input/ output/
