import os
from typing import List, Generator, Set
from FEV_KEGG.settings import cachePath

import pickle
from FEV_KEGG.Util import Parallelism
import contextlib
import tempfile
import functools

from io import IOBase


def createPath(fileName):
    """
    Extracts the path from fileName and ensures all folders exist.
    
    `fileName` is relative to your cache folder! See :attr:`FEV_KEGG.settings.cachePath`.
    
    Parameters
    ----------
    fileName : str
        Path and file name in a format your OS understands. Something like 'subfolder/another_folder/myFile.txt' should most likely work.
    """
    os.makedirs(os.path.dirname(os.path.join(cachePath, fileName)), exist_ok=True)

def writeToFile(content: 'will be encoded in UTF-8', fileName: 'will be overwritten, if already present', atomic = True):
    """
    Writes `content` to a text file.
    
    `fileName` is relative to your cache folder! See :attr:`FEV_KEGG.settings.cachePath`.
    File will be overwritten completely, if already present.
    
    Parameters
    ----------
    content : str
        Content of the file. Will be encoded into UTF-8.
    fileName : str
        Path and name of the file, in a format your OS understands. Something like 'subfolder/another_folder/myFile.txt' should most likely work.
    atomic : bool, optional
        If *True*, write file atomically.
    
    Raises
    ------
    ValueError
        Encoding into UTF-8 failed.
    OSError
        File could not be opened.
    """
    createPath(fileName)
    path = os.path.join(cachePath, fileName)
    if atomic:
        with atomic_write(path, text = True) as file:
            file.write(content)
    else:
        with open(path, 'w', encoding = 'utf_8', errors = 'strict') as file:
            file.write(content)

def writeToFileBytes(data, fileName: 'will be overwritten, if already present', atomic = True):
    """
    Writes `content` to a binary file.
    
    `fileName` is relative to your cache folder! See :attr:`FEV_KEGG.settings.cachePath`.
    File will be overwritten completely, if already present.
    
    Parameters
    ----------
    data : bytes or bytearray
        Content of the file. Will not be encoded.
    fileName : str
        Path and name of the file, in a format your OS understands. Something like 'subfolder/another_folder/myFile.txt' should most likely work.
    atomic : bool, optional
        If *True*, write file atomically.
        
    Raises
    ------
    OSError
        File could not be opened.
    """
    createPath(fileName)
    path = os.path.join(cachePath, fileName)
    if atomic:
        with atomic_write(path, text = False) as file:
            file.write(data)
    else:
        with open(path, 'wb') as file:
            file.write(data)

def readStringFromFileAtOnce(fileName) -> str:
    """
    Reads from a text file.
    
    `fileName` is relative to your cache folder! See :attr:`FEV_KEGG.settings.cachePath`. Reads the file fully at once (more memory intensive).
    
    Parameters
    ----------
    fileName : str
        Path and name of the file, in a format your OS understands. Something like 'subfolder/another_folder/myFile.txt' should most likely work.
    Returns
    -------
    str
        Content of the file. Will be decoded from UTF-8.
    
    Raises
    ------
    ValueError
        Decoding from UTF-8 failed.
    OSError
        File could not be opened.
    """
    with open(os.path.join(cachePath, fileName), 'r', encoding = 'utf_8', errors = 'strict') as file:
        
        string = ''.join([line for line in file])
        return string

def readBytesFromFileAtOnce(fileName) -> bytes:
    """
    Reads from a binary file.
    
    `fileName` is relative to your cache folder! See :attr:`FEV_KEGG.settings.cachePath`. Reads the file fully at once (more memory intensive).
    
    Parameters
    ----------
    data : bytes or bytearray
        Content of the file. Will not be encoded.
    fileName : str
        Path and name of the file, in a format your OS understands. Something like 'subfolder/another_folder/myFile.txt' should most likely work.
    atomic : bool, optional
        If *True*, write file atomically.
    
    Returns
    -------
    bytes
        Content of the file. Not decoded.
    
    Raises
    ------
    OSError
        File could not be opened.
    """
    with open(os.path.join(cachePath, fileName), 'rb') as file:
        return file.read()

def readListFromFileAtOnce(fileName) -> List[str]:
    """
    Reads list from a text file.
    
    `fileName` is relative to your cache folder! See :attr:`FEV_KEGG.settings.cachePath`. Reads the file fully at once (more memory intensive). Does **NOT** return the newline characters.
    
    Parameters
    ----------
    fileName : str
        Path and name of the file, in a format your OS understands. Something like 'subfolder/another_folder/myFile.txt' should most likely work.
    Returns
    -------
    List[str]
        Content of the file as a list. Will be decoded from UTF-8.
    
    Raises
    ------
    ValueError
        Decoding from UTF-8 failed.
    OSError
        File could not be opened.
    """
    with open(os.path.join(cachePath, fileName), 'r', encoding = 'utf_8', errors = 'strict') as file:
        outputList = []
        for line in file:
            outputList.append(line.rstrip())
        return outputList
    
    
def readSetFromFileAtOnce(fileName) -> Set[str]:
    """
    Reads set from a text file.
    
    `fileName` is relative to your cache folder! See :attr:`FEV_KEGG.settings.cachePath`. Reads the file fully at once (more memory intensive). Does **NOT** return the newline characters.
    
    Parameters
    ----------
    fileName : str
        Path and name of the file, in a format your OS understands. Something like 'subfolder/another_folder/myFile.txt' should most likely work.
    Returns
    -------
    Set[str]
        Content of the file as a set. Will be decoded from UTF-8.
    
    Raises
    ------
    ValueError
        Decoding from UTF-8 failed.
    OSError
        File could not be opened.
    """
    with open(os.path.join(cachePath, fileName), 'r', encoding = 'utf_8', errors = 'strict') as file:
        outputList = set()
        for line in file:
            outputList.add(line.rstrip())
        return outputList
    

def readGeneratorFromFileLinewise(fileName) -> Generator[str, None, None]:
    """
    Reads string-generator from a text file.
    
    `fileName` is relative to your cache folder! See :attr:`FEV_KEGG.settings.cachePath`. Returns a generator (can only be read once!) for reading a file line by line (less memory intensive). Does **NOT** return the newline characters.
    
    Parameters
    ----------
    fileName : str
        Path and name of the file, in a format your OS understands. Something like 'subfolder/another_folder/myFile.txt' should most likely work.
        
    Yields
    -------
    str
        Line of content of the file. Will be decoded from UTF-8.
    
    Raises
    ------
    ValueError
        Decoding from UTF-8 failed.
    OSError
        File could not be opened.
    """
    with open(os.path.join(cachePath, fileName), 'r', encoding = 'utf_8', errors = 'strict') as file:
        for line in file:
            yield line.rstrip()


def doesFileExist(fileName) -> bool:
    """
    Does `fileName` already exist AND is a file?
    
    `fileName` is relative to your cache folder! See :attr:`FEV_KEGG.settings.cachePath`.
    
    Parameters
    ----------
    fileName : str
        Path and name of the file, in a format your OS understands. Something like 'subfolder/another_folder/myFile.txt' should most likely work.
    
    Returns
    -------
    bool
        Whether the file already exists.
    """
    return os.path.isfile(os.path.join(cachePath, fileName))

def doesFolderExist(folderName) -> bool:
    """
    Does `fileName` already exist AND is a folder?
    
    `fileName` is relative to your cache folder! See :attr:`FEV_KEGG.settings.cachePath`.
    
    Parameters
    ----------
    fileName : str
        Path and name of the folder, in a format your OS understands. Something like 'subfolder/another_folder' should most likely work.
    
    Returns
    -------
    bool
        Whether the folder already exists.
    """
    return os.path.isdir(os.path.join(cachePath, folderName))


def getFileHandleRead(fileName) -> IOBase:
    """
    A low-level read-only handle to a text file.
    
    `fileName` is relative to your cache folder! See :attr:`FEV_KEGG.settings.cachePath`.
    
    Parameters
    ----------
    fileName : str
        Path and name of the file, in a format your OS understands. Something like 'subfolder/another_folder/myFile.txt' should most likely work.
    
    Returns
    -------
    IOBase
        A handle to the open file. Be sure to close it eventually!
    
    Raises
    ------
    ValueError
        Decoding from UTF-8 failed.
    OSError
        File could not be opened.
    
    Warnings
    --------
    Make sure your code will close this handle eventually! Else, you will waste a lot of resources and eventually prevent your code from opening any other resources!
    """
    return open(os.path.join(cachePath, fileName), 'r', encoding = 'utf_8', errors = 'strict')


def getFileHandleWrite(fileName) -> IOBase:
    """
    A low-level read-write handle to a text file.
    
    `fileName` is relative to your cache folder! See :attr:`FEV_KEGG.settings.cachePath`.
    File will be overwritten completely, if already present.
    
    Parameters
    ----------
    fileName : str
        Path and name of the file, in a format your OS understands. Something like 'subfolder/another_folder/myFile.txt' should most likely work.
    
    Returns
    -------
    IOBase
        A handle to the open file. Be sure to close it eventually!
    
    Raises
    ------
    ValueError
        Decoding from UTF-8 failed.
    OSError
        File could not be opened.
    
    Warnings
    --------
    Make sure your code will close this handle eventually! Else, you will waste a lot of resources and eventually prevent your code from opening any other resources!
    """
    return open(os.path.join(cachePath, fileName), 'w', encoding = 'utf_8', errors = 'strict')







def cache(folder_path, file_name):
    """
    Decorator for caching files on disk.
    
    Checks if result of wrapped function has already been cached into the file specified by :attr:`FEV_KEGG.settings.cachePath`/`folder_path`/`file_name`.
    If yes, read file and return content.
    If no, execute wrapped function, write result to file, and then return result.
    
    Parameters
    ----------
    folder_path : str
        Path of the file, in a format your OS understands. Something like 'subfolder/another_folder/' should most likely work. Remember, this is relative to your :attr:`FEV_KEGG.settings.cachePath`!
    file_name : str
        Name of the file itself.
    
    Returns
    -------
    decorator
        A decorator to automatically cache the decorated function's return value into the file specified by `folder_path`/`file_name`.
    
    Raises
    ------
    Error
        If the decorated function raises any error. Only relevant if the decorated function is executed at all, which it only is if the result has not already been cached.
    OSError
        File could not be opened.
    """
    def decorator(func):
        @functools.wraps(func)
        def caching(*args):

            relativePath = os.path.join(folder_path, file_name)
            fullPath = os.path.join(cachePath, relativePath)
            
            if doesFileExist(relativePath) is True:
                with open(fullPath, 'rb') as file:
                    try:
                        content = pickle.load(file)
                    except Exception:
                        print("\n File causing the exception: " + fullPath + "\n")
                        raise
                return content
            else:
                
                createPath(fullPath) # create folders in path
                result = func(*args)
                with open(fullPath, 'wb') as file:
                    pickle.dump(result, file, protocol = pickle.HIGHEST_PROTOCOL)
                return result

        return caching

    return decorator

class CacheEntry(object):
    
    def __init__(self, isCached: bool, absolutePath: str, result = None):
        """
        Abstract entry of the cache.
        
        While :func:`cache` returns the content of a cache entry directly, :func:`cacheEntry` returns an instance of this class.
        This allows for delayed disk access, because results are always computed within :func:`cache` or :func:`cacheEntry`, but using this class, the actual reading/writing of the cache file is only started upon executing :meth:`getResult`.
        
        Parameters
        ----------
        isCached : bool
            *True* if the underlying cache file already exists and, thus, will be read, not written.
        absolutePath : str
            The absolute path of the underlying cache file.
        result : Object, optional
            Contains the result of the computation to be written into the underlying cache file. Obviously only necessary if `isCached` == *False*.
        
        Note
        ----
        This is especially useful when calculating results within a process pool, to utilise parallelism, but writing results within a thread pool in the main process, keeping the process pool from idling, while waiting for the disk.
        All this is necessary, because neither :class:`ProcessPoolExecutor` nor :class:`ThreadPoolExecutor` support pre-emptive scheduling.
        The caveat of this approach is that the result of the computation, which might be several megabytes, has to be transferred between the background process and the main process via IPC, which is relatively slow.
        Still, on my laptop, this approach proved to be significantly faster overall.
        """
        self.isCached = isCached
        self.absolutePath = absolutePath
        self.result = result
        
    def _readFile(self):
        with open(self.absolutePath, 'rb') as file:
            content = pickle.load(file) 
        return content
    
    def _writeFile(self):
        createPath(self.absolutePath) # create folders in path
        with open(self.absolutePath, 'wb') as file:
            pickle.dump(self.result, file, protocol = pickle.HIGHEST_PROTOCOL)
    
    def getResult(self, noDiskIO=False):
        """
        Result of the cached computation.
        
        Automatically performs read or write operation in this thread, depending on whether the computation was already cached or has only been computed just now.
        
        Parameters
        ----------
        noDiskIO : bool, optional
            If *True*, the result currently available in memory is returned, without any disk access. If the computation was already cached, i.e. an underlying cache file exists, this is most likely *None*! Only use if you know exactly what you are doing!
        
        Returns
        -------
        Object
            Result of the cached computation, no matter if it came from a cache file or from a fresh computation.
        
        Raises
        ------
        OSError
            File could not be opened. Only relevant if `noDiskIO` == False.
        """
        if Parallelism.getShallCancelThreads() is True:
            import concurrent.futures
            raise concurrent.futures.CancelledError()
        
        if noDiskIO is True: # shall not perform disk I/O
            return self.result
        
        elif self.isCached is True: # result is in cache
            return self._readFile()
        
        else: # result is not in cache
            self._writeFile()
            return self.result

def cacheEntry(folder_path, file_name):
    """
    Decorator for caching files on disk, but returning only an intermediate :class:`CacheEntry`.
    
    Similar to cache decorator, but does not read content from file nor write result to file.
    Checks if result of wrapped function has already been cached into the file specified by :attr:`FEV_KEGG.settings.cachePath`/`folder_path`/`file_name`.
    If yes, return appropriate :class:`CacheEntry`. This object provides all information and methods necessary to read/write the desired data from cache.
    If no, execute wrapped function, and then return appropriate :class:`CacheEntry`.
    
    Parameters
    ----------
    folder_path : str
        Path of the file, in a format your OS understands. Something like 'subfolder/another_folder/' should most likely work. Remember, this is relative to your :attr:`FEV_KEGG.settings.cachePath`!
    file_name : str
        Name of the file itself.
    
    Returns
    -------
    decorator
        A decorator to automatically return the delayed disk access in a :class:`CacheEntry` instance.
    
    Raises
    ------
    Error
        If the decorated function raises any error. Only relevant if the decorated function is executed at all, which it only is if the result has not already been cached.
    OSError
        File could not be opened.
    
    Note
    ----
    This is especially handy when splitting computation (in this function) from disk I/O (in :class:`CacheEntry`) in a parallel computation environment.
    """
    def decorator(func):
        @functools.wraps(func)
        def caching(*args):

            relativePath = os.path.join(folder_path, file_name)
            fullPath = os.path.join(cachePath, relativePath)
               
            if doesFileExist(relativePath) is True:
                return CacheEntry(isCached=True, absolutePath=fullPath)
            
            else:
                result = func(*args)
                return CacheEntry(isCached=False, absolutePath=fullPath, result=result)

        return caching

    return decorator

@contextlib.contextmanager
def atomic_write(filename, text=True, keep=False,
                 suffix='.bak', prefix='tmp'):
    """
    Context manager for overwriting a file atomically.
    
    Usage:
    
    ::
    
        with atomic_write("myfile.txt") as f:
            f.write("data")
    
    The context manager opens a temporary file for writing in the same
    directory as `filename`. On cleanly exiting the with-block, the temp
    file is renamed to the given filename. If the original file already
    exists, it will be overwritten and any existing contents replaced.
    (On POSIX systems, the rename is atomic. Other operating systems may
    not support atomic renames, in which case the function name is
    misleading.)
    
    If an uncaught exception occurs inside the with-block, the original
    file is left untouched. By default the temporary file is not
    preserved. To keep the temp file, pass `keep=True`. Any errors in 
    deleting the temp file are ignored.
    By default, the temp file is opened in text mode. To use binary mode,
    pass `text=False` as an argument. On some operating systems, this make
    no difference.
    
    By default, the temp file will have a name starting with "tmp" and
    ending with ".bak". You can vary that by passing strings as the
    `suffix` and `prefix` arguments.
    
    Note
    ----
    Copyright (c) 2017 ActiveState Software Inc. 
    
    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
     
    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
     
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

    """
    # Copyright (c) 2017 ActiveState Software Inc.
    # 
    # Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
    # 
    # The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
    # 
    # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

    path = os.path.dirname(filename)
    fd, tmp = tempfile.mkstemp(suffix=suffix, prefix=prefix, dir=path, text=text)
    replace = os.replace  # Python 3.3 and better.
    
    try:
        if text is True:
            with os.fdopen(fd, 'w', encoding = 'utf_8', errors = 'strict') as f:
                yield f
        else:
            with os.fdopen(fd, 'wb') as f:
                yield f
        # Perform an atomic rename (if possible). This will be atomic on 
        # POSIX systems, and Windows for Python 3.3 or higher.
        replace(tmp, filename)
        tmp = None

    finally:
        if (tmp is not None) and (not keep):
            # Silently delete the temporary file. Ignore any errors.
            try:
                os.unlink(tmp)
            except:
                pass
