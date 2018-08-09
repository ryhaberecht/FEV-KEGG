from FEV_KEGG import settings
    
def startProcessPool():
    """
    Shortcut for starting the process pool for parallel computation.
    
    See Also
    --------
    FEV_KEGG.Util.Parallelism.startProcessPool : The actual function doing the work.
    """
    from FEV_KEGG.Util import Parallelism
    Parallelism.startProcessPool()
    
if settings.automaticallyStartProcessPool:
    startProcessPool()