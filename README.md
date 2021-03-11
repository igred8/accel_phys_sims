Accelerator Physics Simulations API
==============================

Python based functionality that handles the control and interaction of some commonly used numerical simulation tools for particle accelerator physics.

The library has the ability to:
 - run elegant, genesis1.3, ICS code from LLNL, and GPT (future) from inside jupyter notebooks or python scripts
 - load SDDS files from elegant and ones produced by GPT (using `gdf2sdds`)
 - manipulate 6D particle distributions and handle I/O

See other readme files for details on setup and use.

Project Organization
------------

    ├── LICENSE
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   └── sim_out        <- Output data from simulations. Subfolders by date sim was run.
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    ├── src                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   │
    │   ├── features       <- Python libraries to turn simulation data into features for modeling
    │   │   └── build_features.py
    

--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
