<!-- title: python-anaconda-vscode-radia -->

# Setting up Anaconda and VS Code
---

1. Download Anaconda Navigator installer
    - Run it. 
    - As recommmended, install for User.
2. In order to reproduce the python settings used for the development of custom libraries related to accelerator physics and such.
    - You can use a terminal to create the environment:
        
        `$ conda create --name accel_phys_38 --file accel_phys_38.yml`
      - Refer to this [cheat sheet](https://docs.conda.io/projects/conda/en/latest/_downloads/843d9e0198f2a193a3484886fa28163c/conda-cheatsheet.pdf)


    - The same thing can be done inside Anaconda Navigator's graphical interface. Use the "Import" feature.
    - Inside Anaconda Navigator, import the python environment by using the `accel_phys_38.yml` file. The file is output by Anaconda looks like this:
    ```
        # This file may be used to create an environment using:
        # $ conda create --name <env> --file <this file>
        # platform: win-64
        @EXPLICIT
        https://repo.anaconda.com/pkgs/main/win-64/blas-1.0-mkl.tar.bz2
        https://conda.anaconda.org/conda-forge/win-64/ca-certificates-2020.6.20-hecda079_0.tar.bz2
        https://repo.anaconda.com/pkgs/main/win-64/icc_rt-2019.0.0-h0cc432a_1.tar.bz2
        ...
        ...
        ...
    ```
3. Download and install VS Code. 
    - Install the python and pylance extensions. 
    - Set the default terminal inside VS Code to use command prompt (CMD). Using PowerShell fails to _activate_ Anaconda environments. Probably has to do with where PowerShell looks for path variables and such.
        - Open a terminal and "select Default Shell"
4. Select python interpreter. Once VS Code is open, select a python interpreter to handle jupyter notebooks and scripts. 
    - `ctrl + shift + P`
    - Begin typing "select python interpreter" and click on the command. The available Anaconda environments should show up. Select one. (Yay! VS Code is smart)

At this point, the VS Code + Anaconda is set up. Setting up jupyter server is omitted, because it is assumed that the native VS Code jupyter notebook viewer is used.

# Setting up Radia for python
---

You must install `git` in order to clone repos from github. However, you may simply go to the repo website, download the repo, and extract it to the directory you want. 

The python version of Radia is available from Prof. Oleg Chubar's github repository: https://github.com/ochubar/Radia.git

 - Open a terminal to the directory where Radia will be housed.
 - Make sure you are connected to the internet. 
 - `$ git clone https://github.com/ochubar/Radia.git`

Now, Radia repo is cloned to the specified directory.

To use pyRadia in a jupyter notebook, the location of the compiled python dll (`radia.pyd`) must be added to the path. Suppose you cloned Radia to the directory `D:/some_folder/`, here is the sample code to import the Radia module:
```
import sys
sys.path.append(r'D:/some_folder/Radia/env/radia_python/')
import radia
```
Now, you may get an error that says the specified DLL cannot be loaded. One reason this may happen is that the compiled `pyd` file may be for a different version of python than that of the loaded environment. What I have found to work is to rename the file that matches my python version to `radia.pyd`. e.g. the env uses python 3.8.x and there is a file named `../Radia/env/radia_python/radia_py3_8_x64.pyd`. Rename it to `radia.pyd` and then use `import radia`. You can also, skip renaming the file, but then you have write out the name in the import line: `import radia_py3_8_x64`, which is annoying.

