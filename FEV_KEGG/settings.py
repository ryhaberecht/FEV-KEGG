"""
Settings used across the application.
"""





### biological ###

defaultEvalue = 1e-15
"""
Default threshold for the statistical expectation value (E-value), below which a sequence alignment is considered significant.
This can usually be overwritten for a whole module, or in a function's parameters.
"""

defaultOneOrganismPerSpecies = True
"""
Whether to use only one organism per species in a taxonomy. This is useful to reduce overrepresentation of species with many sequenced genomes in a group. 
"""

defaultNoMultifunctional = True
"""
Whether to exclude multifunctional enzymes from analysis. This is useful to exclude complex cases of multi-domain enzymes and known promiscuous enzymes.
"""




### technical ###

import appdirs
cachePath = appdirs.user_cache_dir(appname = 'FEV-KEGG', appauthor = 'ryh') # according to your OS
"""
File system path for the cache directory.

This directory is used for caching downloads and calculations. It is automatically set depending on your OS's default cache folder.

See Also
--------
:ref:`readme-cache-reference` : Readmy entry explaining caching in general.
"""

verbosity = 1
"""
Verbosity of the output on the console.

0 = no unecessary output, only errors.
1 = informative output, including progess for long-running tasks.
2 = including short description of the process shown as a progress bar, e.g. 'Calculating redundancy...'
3 = debug output, including whether a file is downloaded or already in cache.
"""


import multiprocessing
try:
    processes = int( multiprocessing.cpu_count() )
except NotImplementedError:
    processes = 1
    """
    Number of processes to calculate with.
    
    Depends on your computer. Defaults to the number of logical CPU cores, or if unretrievable, to 1.
    
    Warnings
    --------
    Processes have much more overhead than threads.
    """

automaticallyStartProcessPool = True
"""
Whether a pool of processes should automatically spin up.

You will want this *False* if:

- you ``import FEV_KEGG`` inside a background thread or a child process *only*, because starting a process pool only works from within the main thread of the main process **OR**
- you already have other means of providing parallelism **OR**
- you want to save on the rather big overhead of spinning up, using, and tearing down processes. Especially when you ``import FEV_KEGG`` without actually using the process pool, e.g. when using something like Sphinx with autodoc.

If you do have this *False*, you will have to either:

- run :func:`FEV_KEGG.startProcessPool()`. You should do this as early as possible, because your whole main process will be copied, including all the RAM it hitherto occupied. Will only work within the main thread of the main process **OR**
- replace :attr:`FEV_KEGG.Util.Parallelism.processPool` with a :class:`concurrent.futures.Executor` of your choice, or else some functions will fail and raise an error, especially in :class:`FEV_KEGG.KEGG.Organism.Group`. The object may be a dummy, as long as it adheres to the interface and successfully works on tasks.

You will want this *True* if:

- you ``import FEV_KEGG`` in your main thread of the main process only **AND**
- you want to use functions utilising parallelism in the process pool, especially in :class:`FEV_KEGG.KEGG.Organism.Group` **AND**
- you rarely import without actually using the process pool, because spinning up and tearing down processes can take several seconds.

See Also
--------
FEV_KEGG.Util.Parallelism.startProcessPool : Spins up a pool of processes. No further setup required.
"""

import math
downloadThreads = 32 # in main process, default: 64
"""
Number of threads to download with at once.

Depends on the speed of your internet connection.

Warnings
--------
Should only be used inside the main process. Else, the total number of download threads will multiply by the number of available processors!
"""

downloadThreadsSSDB = 12 # in main process, default: 16
"""
Number of threads to download with at once from KEGG SSDB.

Depends on the speed of your internet connection.
Also depends on SSDB, because this part of KEGG can obviously withstand far fewer parallel connections.
If you get frequent HTTP 403 Forbidden errors, despite of retrying upon such errors, lower this value until the progress bar moves smoothly.

Warnings
--------
Should only be used inside the main process. Else, the total number of download threads will multiply by the number of available processors!
"""

downloadThreadsPerProcess = math.floor( downloadThreads / processes ) # in sub-processes
"""
Number of threads to download with at once, divided by the number of available processors.

Should be used inside background processes. Depends on your computer.
"""

downloadThreadsPerProcessSSDB = math.floor( downloadThreadsSSDB / processes ) # in sub-processes
"""
Number of threads to download with at once from KEGG SSDB, divided by the number of available processors.

Should be used inside background processes. Depends on your computer.
"""

fileThreads = 8 # in main process, default: 8
"""
Number of threads to access files with at once.

Depends on the speed of your data storage.

Warnings
--------
Should only be used inside the main process. Else, the total number of file threads will multiply by the number of available processors!
"""
fileThreadsPerProcess = math.floor( fileThreads / processes ) # in sub-processes
"""
Number of threads to acces files with at once, divided by the number of available processors.

Should be used inside background processes. Depends on your computer.
"""


# parameters for retrying failed downloads
downloadTimeout = 10
"""
Time in seconds after which a blocking socket is considered failed.

Should be well below `retryDownloadMax`, or else you will not retry at all.
"""

downloadTimeoutSocket = math.floor(downloadTimeout / 2)
"""
Quirks for `downloadTimeout` in Python 3.4 (and maybe above).

For some reason, in Python 3.4, sockets double the timeout value passed to them.
"""

retryDownloadMax = 60000 # default: 60000
"""
Maximum total time in milliseconds a download should be retried before giving up.

Warnings
--------
Giving up usually throws an error and halts the whole program!
"""

retryDownloadBackoffFactor = 1000 # default: 1000
"""
Factor in milliseconds for exponential backoff.

The backoff function waits ``2^x * retryDownloadBackoffFactor`` milliseconds between retries, where *x* is the count of tries already failed.
"""

retryDownloadBackoffMax = 10000 # default: 10000
"""
Maximum time between two retries, which the exponential backoff function can not exceed.
"""

organismInfoExpiration = 86400 # 24h
"""
Seconds after which an organism info file is considered outdated. This can be useful when relying upon a current database size for calculating E-values for a :class:`FEV_KEGG.KEGG.SSDB.Match`.
"""