# elegant - SDDS - python

2019.02.06

This is a readme file that documents my installation of elegant on my personal machine (windwos 10).

---

## elegant

1. Download the elegant.msi file from the [APS website](https://www.aps.anl.gov/Accelerator-Operations-Physics/Software#elegant).
2. Run it and install elegant to desired location (e.g. C:\Program Files\APS\Elegant x64).
3. Add this location to the Windows path.
    - `Control Panel > System > Advanced system settings > Environment variables...`
    - Inside 'System variables' find 'Path'. Select it and click 'Edit'
    - Click 'New' (or 'Edit' if you are changing the location of elegant)
    - Write the path to the location of the `elegant.exe` file.

After performing these steps, it should be possible to call elegant from the command prompt (or PowerShell).

```PowerShell
PS D:\Dropbox\RBT\beam_condi\elegant_sims>elegant
```

---

## SDDS and python

SDDS is a file type for data that was developed by APS (advanced photon source at ANL). It is the main output of elegant. There is an SDDS toolkit, which handles these files.

The general idea of any of the functions of the SDDS toolkit is that it should perform an operation on the data and return the same structure, so that another function can take the output and perform its modifications. In this way, a sequence of toolkit functions can be strung together to produce a final result ( R = Ofinal...O3.O2.O1.data ).

The Python SDDS binaries aim to translate the SDDS toolkit into a pythonic mudule.

!!! At the time of writing (2020.11.19) Python 3.8 is supported. !!!

1. The python library can be downloaded here: [SDDS Python3.8.msi](https://www.aps.anl.gov/Accelerator-Operations-Physics/Software)
2. Run the `.msi` file and install.
   1. This will make a directory (default is `c:\Pythonxx\`) that will contain the following:

        ```filestructure
            ..\
            - DLLs
                - .dll files
                - sddsdata.pyd (DLL called by sdds.py)
            - Lib
                - sdds.py (SDDS class definition)
            - SDDS_demo (demo on how to use SDDS for python)
        ```

   2. If Python is indeed installed as a standalone, then the .dll, sddsdata.pyd (which is also a DLL), and the sdds.py files would be placed in the DLLs and Lib folders of the installation. However, I am using Anaconda3 as my python environment manager. This means that the DLLs and Lib folders are not found in the default path `c:\Pythonxx`.
       1. Copy and paste the .dll and .pyd files into the directory: `c:\Anaconda3\envs\accel_phys\DLLs\`. The generic form is `c:\Anaconda3\envs\name_of_venv\DLLs\`.
       2. Copy and paste the `sdds.py` file in the directory `c:\Anaconda3\envs\accel_phys\Lib\`. The generic form is `c:\Anaconda3\envs\name_of_venv\Lib\`.
       3. Now, it will be possible to simply `import sdds` for the venv 'accel_phys' or 'name_of_venv'.
       4. No need to add path variables to `sys.path` and install packages with setuptools.

---
useful links

[python packages - section 6.4](https://docs.python.org/3/tutorial/modules.html#the-module-search-path)

[packaging python projects](https://packaging.python.org/tutorials/packaging-projects/)

[Stack overflow about SDDS and python by a physicist at Stanford](https://stackoverflow.com/questions/28222655/importing-class-in-python-subpackage-imports-more-than-requested)

---

Created by Ivan Gadjev. 2019.02.06