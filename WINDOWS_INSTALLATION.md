# Installation For Windows Computers

### 1. Install Visual Studio C++ 2017 :
 1. Download the correct version of Microsoft Visual C++ Redistributable for Visual Studio 2017 https://support.microsoft.com/en-us/help/2977003/the-latest-supported-visual-c-downloads
 * If Windows is 64 bit, download and run  "vc_redist.x64.exe"
 * If Windows is 32 bit, download and run "vc_redist.x86.exe"
2. Restart computer

### 2. Install Python 3.6 and required packages:
1. Download and Run https://repo.anaconda.com/archive/Anaconda3-2019.03-Windows-x86_64.exe
2. During Installation, choose the option to add CONDA to PATH
3. After Installation, open command prompt
4. Go to directory of the timstof package
5. Run the following command
    ````
        conda env create -f environment.yml
    ````

### 3. Configure Python Scripts
* Update the library paths in src/timsdata.py so that it points to the correct binary files in lib
* If Windows is 64 bit: use "lib/win64/timsdata.dll"
* If Windows is 32 bit: use "lib/win32/timsdata.dll"
* Example:
    ```
    if sys.platform[:5] == "win32":
        libname = "[put path here]/lib/win64/timsdata.dll"
    ```

# To Run
```
conda activate timstof3
python extract_msn_nopd.py input/ output/
```
