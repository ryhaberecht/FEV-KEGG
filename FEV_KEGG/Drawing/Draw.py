import networkx.classes
import networkx.drawing.layout
from FEV_KEGG.Graph import Models
from FEV_KEGG.KEGG import File
from enum import Enum

def toPNG(graph: Models.CommonGraphApi, fileName: 'path/file', layout = 'neato'):
    """
    Draw `graph` and save it as PNG format in a file.
    
    Drawing requires an algorithm to determine the best position für each node and the best path for each edge. This algorithm can be defined via `layout`.
    
    Parameters
    ----------
    graph : Models.CommonGraphApi
        The graph to be drawn.
    fileName : str
        Path and name of the file, the extension '.png' is automatically applied. The path is relative to the current working directory!
    layout : str, optional
        A layout algorithm known to :mod:`pygraphviz`.
    
    Raises
    ------
    ImportError
        If :mod:`pygraphviz` is not installed. This is an optional dependency and **not** installed via pip by deault! PyGraphviz needs Graphviz to function, which is not a python program and has to be installed manually by you!
    NotImplementedError
        If `graph` is not of a NetworkX type.
    """
    
    nxGraph = graph.underlyingRawGraph
    if isinstance(nxGraph, networkx.classes.graph.Graph):
        
        File.createPath(fileName)
        agraph = networkx.nx_agraph.to_agraph(graph.underlyingRawGraph)
        agraph.layout()
        agraph.draw(fileName + '.png', format = 'png', prog = layout, args = '-Tpng')
        
    else:
        raise NotImplementedError

class NetworkxLayout(Enum):
    """
    Enum of layout algorithms known to NetworkX.
    """
    kamada_kawai = networkx.drawing.layout.kamada_kawai_layout
    random = networkx.drawing.layout.random_layout
    shell = networkx.drawing.layout.shell_layout
    spring = networkx.drawing.layout.spring_layout
    spectral = networkx.drawing.layout.spectral_layout
    fruchterman_reingold = networkx.drawing.layout.fruchterman_reingold_layout

def toWindow(graph: Models.CommonGraphApi, layout: NetworkxLayout):
    """
    Draw `graph` and display it in a window.
    
    Drawing requires an algorithm to determine the best position für each node and the best path for each edge. This algorithm can be defined via `layout`.
    
    Parameters
    ----------
    graph : Models.CommonGraphApi
        The graph to be drawn.
    layout : NetworkxLayout
        A layout algorithm known to :mod:`networkx`, as defined in :class:`NetworkxLayout`.
    
    Raises
    ------
    ImportError
        If :mod:`matplotlib` is not installed. This is an optional dependency and **not** installed via pip by deault! Matplotlib needs a working backend to function, which is not a python program and has to be installed manually by you! For a list of backends, see `Matplotlib's website <https://matplotlib.org>`_.
    NotImplementedError
        If `graph` is not of a NetworkX type.
    """
    import matplotlib  # @UnresolvedImport
    nxGraph = graph.underlyingRawGraph
    if isinstance(nxGraph, networkx.classes.graph.Graph):
        positions = layout(nxGraph)
        networkx.drawing.nx_pylab.draw(nxGraph, pos = positions)
        
        edges = nxGraph.edges(keys = True)

        labelDict = dict()
        for n1, n2, key in edges:
            labelDict[(n1, n2)] = key.__str__()
        
        networkx.drawing.nx_pylab.draw_networkx_edge_labels(nxGraph, pos = positions, edge_labels = labelDict)
        matplotlib.pyplot.show()
        
    else:
        raise NotImplementedError