README
======

FEV\@KEGG
---------
FEV\@KEGG allows for easy analysis of metabolic networks of organisms in KEGG (Kyoto Encyclopedia of Genes and Genomes).
Read the API documentation here: https://fev-kegg.readthedocs.io


Restrictions
____________
- You **MUST** make absolutely sure to comply with the conditions of using KEGG and its API: http://www.kegg.jp/kegg/legal.html and http://www.kegg.jp/kegg/rest/.
- If you have access to an offline copy of KEGG, you **MUST NOT** use the default Database and Download modules, since they cause a lot of load on KEGG servers. Instead, contact me, so we can integrate your offline copy to be used before anything is downloaded.


Features
________
- convert data from KEGG PATHWAY and KEGG GENE to organism-specific graphs
- graphs link substrates/products with reactions, genes, EC numbers, or abstract 'enzymes'
- cache downloads from KEGG, graphs, and any other computational result
- build groups of organisms, allows for fusing their graphs into a common metabolism
- gather groups from NCBI or KEGG taxonomy, using KEGG BRITE
- gather clades from NCBI taxonomy and compare their 'core' metabolism
- find paralogs/orthologs using KEGG SSDB
- find possible gene duplications or neofunctionalisations
- calculate robustness metrics between graphs, organisms, groups of organisms, or clades
- ... anything you can think of using graphs derived from KEGG


Install
-------
Use pip to install FEV\@KEGG and to automagically install all dependencies:
``pip install FEV_KEGG``

If you are on Python 3.4, you will have to use ``pip install FEV_KEGG[python34]`` to pull in the backported *typing* package.


Where to start?
---------------
After successful installation, you might want to take a look at the "experiment" scripts in *FEV_KEGG/Experiments*.
These scripts consecutively involve more and more functionality of this library. They were used during development, step by step adding and testing another layer of functionality or abstraction.
Therefore, they might be useful to you in learning how to use this very functionality.

Also, take a look at the API documentation: https://fev-kegg.readthedocs.io

If any questions remain, feel free to report an issue: https://github.com/ryhaberecht/FEV-KEGG/issues


Dependencies
------------
These are automatically installed by pip.

- Python 3.4+
- NetworkX
- anytree
- jsonpickle
- tqdm
- BeautifulSoup
- retrying
- appdirs
- typing (for Python 3.4 only)


Optional Dependencies
---------------------
If you want to draw a graph to an image file:

- PyGraphviz
- Graphviz (non-python software you will have to install manually!)

Use ``pip install FEV_KEGG[draw_image]``.

|

If you want to draw a graph in a pop-up window:

- Matplotlib

Use ``pip install FEV_KEGG[draw_window]``

|

Exporting to GraphML or GML works without any optional dependencies.


Included Dependencies
---------------------
These have been partially copied into this project to avoid unnecessarily big dependencies and allow for minor changes.

- Bio.KEGG from Biopython in lib.Biopython.KEGG


Recommendations
---------------
- SSD

When handling 500 organisms from KEGG at once:

- 64 bit operating system
- 4 GB RAM
- 20 GB disk space for cache

When handling all ~5000 organisms in KEGG at once:

- 64 bit operating system
- 12 GB RAM
- 100 GB disk space for cache


Developer's System
------------------
- cPython 3.4.6
- x86-64 Linux (OpenSUSE Leap 42.3)
- 16 GB RAM
- 1 CPU, 2 Cores, 4 Threads
- SSD

.. _readme-cache-reference:

Caching
-------
- The cache directory path is set up in the 'settings.py' file on the top level of the project. Per default, it points to your user's cache directory as defined by your OS.

  - Linux/Unix: ~/.cache/FEV-KEGG
  
  - OS X: ~/Library/Caches/FEV-KEGG
  
  - Windows: C:\\Users\\username\\AppData\\Local\\ryh\\FEV-KEGG\\Cache
  
- All downloads from KEGG are cached automatically. Also, basic graphs are cached by organism. These default cachings alone can grow the cache directory to 100 GB size!
- You can cache any function's result using the @cache decorator, see :func:`FEV_KEGG.KEGG.File.cache`. Watch out to remember the path and file name and not to overwrite any other cached files.
- To cause a download of the newest version of data from KEGG, you have to delete the cached file manually. Have a look inside the 'cache' folder, file paths and names should be self-explanatory.
- On Linux with supporting file systems, disabling atime (file access time) for the cache directory and all its contents might improve performance: sudo chattr -R +A ~/.cache/FEV-KEGG
