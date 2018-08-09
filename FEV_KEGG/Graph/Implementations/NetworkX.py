import networkx.classes


class NetworkX(object):
    """
    Represents the NetworkX library
    """

class MultiDiGraph(networkx.classes.multidigraph.MultiDiGraph, NetworkX):
    """
    Represents a MultiDiGraph from NetworkX
    """
    
class MultiGraph(networkx.classes.multigraph.MultiGraph, NetworkX):
    """
    Represents a UndirectedMultiGraph from NetworkX
    """