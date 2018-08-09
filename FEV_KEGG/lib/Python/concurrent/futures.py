import concurrent.futures
from concurrent.futures.process import _ResultItem, os
import sys
import traceback

def process_worker(call_queue, result_queue):
    """
    A copy of Python's process_worker function in :class:`concurrent.futures.ProcessPoolExecutor`.
    
    This copy was changed to not die on KeyboardInterrupt, but to exit gracefully.
    Also, no traceback is printed upon :class:`NoKnownPathwaysError` or :class:`CancelledError`.
    
    Note
    ----
    Copyright © 2001-2018 Python Software Foundation; All Rights Reserved
    """
    while True:
        try:
            call_item = call_queue.get(block=True)
            if call_item is None:
                # Wake up queue management thread
                result_queue.put(os.getpid())
                return
        except KeyboardInterrupt as e:
            sys.exit()
            
        try:
            r = call_item.fn(*call_item.args, **call_item.kwargs)
        except BaseException as e:
            from FEV_KEGG.KEGG.Database import NoKnownPathwaysError
            if isinstance(e, NoKnownPathwaysError) or isinstance(e, concurrent.futures._base.CancelledError):
                pass
            elif isinstance(e, Exception):
                traceback.print_exc()
            
            result_queue.put(_ResultItem(call_item.work_id,
                                         exception=e))
        else:
            result_queue.put(_ResultItem(call_item.work_id,
                                         result=r))
        
            


class InterruptibleProcessPoolExecutor(concurrent.futures.ProcessPoolExecutor):
    """
    This class inherits from :class:`concurrent.futures.ProcessPoolExecutor` and only replaces the worker function which pulls work items from the call_queue and reports their results back on the result_queue.
    A KeyboardInterrupt can be caught during execution of the work item, but not when the next work item is being pulled from the call_queue! A KeyboardInterrupt stack trace is printed and the process dies.
    This class catches a KeyboardInterrupt during call_queue-pulling and causes the current process (one within the pool) to sys.exit() immediately, without printing a KeyboardInterrupt stack trace.
    Per default, sys.exit() does not cause a stack trace print.
    """
    concurrent.futures.process._process_worker = process_worker