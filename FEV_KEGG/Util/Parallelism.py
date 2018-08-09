def printBelowProgress(message):
    """
    Does a print() beneath any progress bars.
    
    Simply prints as many newlines as there are processes in use by the process pool.
    Finally, prints the `message`.
    
    Parameters
    ----------
    message : str
        Message to be printed below any progress bars.
    
    Note
    ----
    This is an indirect measure of currently possibly rendered progress bars, at best!
    There is no check if the possible progress bars are actually rendered.
    Therefore, a surplus of empty lines may appear, or there may not be enough.
    """
    import FEV_KEGG.settings
    for _ in range(0, FEV_KEGG.settings.processes + 1):
        print('')
    if message is not None:
        print( message )

def isMainProcess():
    """
    Check if this is running in the main process.
    
    Returns
    -------
    bool
        *True* if this is the main process (any thread), *False* if not.
    """
    import multiprocessing
    return not type(multiprocessing.current_process()) == multiprocessing.Process

def isMainThread():
    """
    Check if this is running in the main thread.
    
    Returns
    -------
    bool
        *True* if this is the main thread (of any process), *False* if not.
    """
    import threading
    return threading.current_thread() is threading.main_thread()

def isMainThreadInMainProcess():
    """
    Check if this is running in the main thread of the main process.
    
    Returns
    -------
    bool
        *True* if this is the main thread of the main process, *False* if not.
    """
    return isMainProcess() and isMainThread()

def keyboardInterruptHandler(processPoolFutures = None, threadPool = None, threadPoolFutures = None, silent = False, terminateProcess = False):
    """
    Handle :class:`KeyboardInterrupt`.
    
    Tries to consistently handle :class:`KeyboardInterrupt` in all current implementations of parallel computation.
    If this is inside a background process, cancel `threadPool` tasks, running and scheduled (if listed inside `threadPoolFutures`).
    If this is inside the main process, also cancel `processPoolFutures` and shut down the global process pool itself.
    May also exit this very process, if `terminateProcess` == *True*.
    
    Parameters
    ----------
    processPoolFutures : Iterable, optional
        An :class:`Iterable` of all :class:`Future` currently scheduled in the process pool.
    threadPool : ThreadPoolExecutor, optional
        The thread pool currently running within the main process. This pool usually handles the process pool's results in a parallel manner.
    threadPoolFutures : Iterable, optional
        The :class:`Iterable` of all :class:`Future` currently scheduled in the `threadPool`.
    silent : bool, optional
        If *True*, prints messages of the current state of handling the interrupt.
    terminateProcess : bool, optional
        If *True*, calls :func:`sys.exit()` on this very process, at the end of this function.
    
    Warnings
    --------
    Only ever call this from main thread (in any process)! Else, you might get unexpected and dangerous behaviour!
    
    Note
    ----
    `terminateProcess` is useful to cancel tasks already scheduled inside the *call_queue* of the pool, because there is no method to empty this queue and cancel already scheduled tasks.
    Keep in mind that this queue is **always** populated with *p* + 1 tasks, while *p* (different) tasks are being executed, with *p* being the number of processes in the pool.
    However, this only works as long as :class:`ProcessPoolExecutor` does **not** re-create broken processes!
    """
    import sys
    
    # tell running threads to cancel. Including the ones in a thread pool's call_queue.
    if terminateProcess is True or threadPool is not None:
        enableShallCancelThreads()
    
    # tell futures scheduled in sub-processes to cancel
    if isMainProcess():
        if silent is False:
            printBelowProgress( 'Program interrupted by keyboard. Exiting, please wait...' )
            print( 'Grinding processes to a halt...' )
        
        # cancel futures. This does NOT include futures already scheduled in the process' call_queue! Which are always n+1 (n = number of processes) futures.
        if processPoolFutures is not None:
            for future in processPoolFutures:
                future.cancel()
        
        # tell the pool to shut down. This does NOT stop already scheduled futures to cancel or even abort! Pending futures are canceled above. Running futures (and futures already inside the call_queue) can NOT be aborted through the current API.
        global processPool
        if processPool is not None:
            processPool.shutdown(wait = not silent)
        
        # Running futures (and futures already inside the call_queue) can NOT be cancelled. The only way is to terminate all processes of the pool, see sys.exit() below.
    
    # tell scheduled threads to cancel
    if threadPool is not None:
        if silent is False: 
            print( 'Knotting together loose threads...' )
        
        if threadPoolFutures is not None:
            for future in threadPoolFutures:
                future.cancel()
                
        threadPool.shutdown(wait = True)
        
        # reset shallCancel signal, so next work item can function correctly
        if terminateProcess is False:
            resetShallCancelThreads()
    
    # terminate process. This allows to cancel futures already inside the call_queue of the pool. Only works because, currently, the pool does not re-create broken processes.
    if terminateProcess is True:
        if silent is True:
            sys.exit()
        else:
            sys.exit('Done. SNAFU')
    
def getNumberOfThreadsDownload(isSSDB = False):
    """
    Number of threads allowed for downloading.
    
    Depends on whether this function is executed inside a process pool, because all processes will have to share the total number of allowed download threads, in order not to exceed it.
    
    Parameters
    ----------
    isSSDB : bool. optional
        IF *True*, uses a lower number of threads, because KEGG SSDB can not handle too many parallel connections.
    
    Returns
    -------
    int
        The number of allowed threads for downloading.
    """
    import FEV_KEGG.settings
    if isMainProcess():
        if isSSDB:
            return FEV_KEGG.settings.downloadThreadsSSDB
        else:
            return FEV_KEGG.settings.downloadThreads
    else:
        if isSSDB:
            return FEV_KEGG.settings.downloadThreadsPerProcessSSDB
        else:
            return FEV_KEGG.settings.downloadThreadsPerProcess
    
def getNumberOfThreadsFile():
    """
    Number of thread allowed for file access.
    
    Depends on whether this function is executed inside a process pool, because all processes will have to share the total number of allowed file access threads, in order not to exceed it.
    
    Returns
    -------
    int
        The number of allowed threads for file access.
    """
    import FEV_KEGG.settings
    if isMainProcess():
        return FEV_KEGG.settings.fileThreads
    else:
        return FEV_KEGG.settings.fileThreadsPerProcess

def getTqdmPosition():
    """
    Get position where this process' progress bar shall be printed.
    
    Interprets current process' name as the index for positioning its :mod:`tqdm` progress bar alongside the progress bar of all other processes.
    
    Returns
    -------
    int
        Position index of this process' progress bar, to be used by :mod:`tqdm`.
    
    Note
    ----
    This function depends on the naming of processes of your OS!
    You might have to replace this function to get clean positioning of progress bars.
    """
    if not isMainProcess():
        import multiprocessing
        import re
        try:
            position = int( re.sub('^.+-(?=[0-9]+$)', '', multiprocessing.current_process().name) )
            return position
        except Exception:
            return 1 
    else:
        return 0

processPool = None
"""
The global process pool.

If you want to run FEV\@KEGG inside another thread or process, your will have to replace this with a :class:`concurrent.futures.Executor` of your choice.

See Also
--------
startProcessPool : Default creator of a process pool.
"""

def startProcessPool():
    """
    Creates the process pool and stores it in :attr:`processPool`.
    
    You might not want this when you already run your own :class:`concurrent.futures.Executor`. If you run your own, let :attr:`processPool` point to it.
    
    See Also
    --------
    FEV_KEGG.lib.Python.concurrent.futures.InterruptibleProcessPoolExecutor : Our own executor used as a process pool, because it handles :class:`KeyboardInterrupt` nicer than the default implementation.
    
    Note
    ----
    A process pool will only work from within the main thread in the main process of your application.
    If you want to run FEV\@KEGG inside another thread or process, your will have to replace :attr:`processPool` with a :class:`concurrent.futures.Executor` that allows running outside the main thread. This may be a dummy class.
    """
    global processPool
    if isMainThreadInMainProcess() is True and processPool is None:
    
        # create process pool at start. Prevents forking a heavy process.
        from FEV_KEGG.settings import processes
        import FEV_KEGG.lib.Python.concurrent.futures
        # InterruptibleProcessPoolExecutor ovverides the process worker function of concurrent.futures.ProcessPoolExecutor, keeping a KeyboardInterrupt during queue-pull from causing a stack trace print
        processPool = FEV_KEGG.lib.Python.concurrent.futures.InterruptibleProcessPoolExecutor( processes )
        
        # terminate process pool on exit. Frees resources.
        import atexit
        atexit.register(lambda x: x.shutdown(wait = False), processPool)    


shallCancelThreads = False
"""
Allows threads to be canceled in a timely manner, across the whole process.

Specifically, any task currently running in a thread of this process can check this variable to see whether it is supposed to abort and return early.

Note
----
This is necessary because :class:`concurrent.futures.ThreadPoolExecutor` does not have support for neither pre-emptively aborting tasks (*call_queue*) nor emptying the pending task queue (*pending_work_items*).
"""

def resetShallCancelThreads():
    """
    Reset the global flag for cancelling all threads to "not cancelling".
    """
    global shallCancelThreads
    shallCancelThreads = False

def enableShallCancelThreads():
    """
    Set the global flag for cancelling all threads to "cancelling".
    """
    global shallCancelThreads
    shallCancelThreads = True

def getShallCancelThreads():
    """
    Get the global flag for cancelling all threads.
    
    Returns
    -------
    bool
        Whether all tasks currently running in threads of this process shall abort and return early.
    """
    global shallCancelThreads
    return shallCancelThreads