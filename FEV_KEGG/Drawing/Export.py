from FEV_KEGG.Graph import Models
import networkx.classes
import os
from enum import Enum
from FEV_KEGG.Graph.SubstanceGraphs import SubstanceEcGraph
import math
from FEV_KEGG import settings

LABEL_NAME = 'custom_label'
"""
Name of the attribute of the graph to be used for storing the label.
"""
COLOUR_NAME = 'colour'
"""
Name of the attribute of the graph to be used for storing the colour.
"""
DESCRIPTION_NAME = 'custom_description'
"""
Name of the attribute of the graph to be used for storing the description.
"""
ID_NAME = 'custom_id'
"""
Name of the attribute of the graph to be used for storing the id.
"""
REACTION_NAME = 'custom_reaction'
"""
Name of the attribute of the graph to be used for storing the reaction.
"""
MAJORITY_NAME = 'custom_majority'
"""
Name of the attribute of the graph to be used for storing the majority percentage.
"""
URL_NAME = 'URL'
"""
Name of the attribute of the graph to be used for storing the URL.
"""

def addLabelAttribute(nxGraph: networkx.classes.MultiGraph):
    """
    Adds the "custom_label" attribute to each node and edge.
    
    Add an attribute to nodes and edges called "custom_label" (see module variable :attr:`LABEL_NAME`) containing the string of the node's label, or the edge's key's label, repectively.
    A label is defined as either the `.name` field, or if that is *None*, the id.
    This is especially useful if the tool you import this file into does not regularly read the XML's id parameter. For example, Cytoscape does not for edges.
    
    Parameters
    ----------
    nxGraph : networkx.classes.MultiGraph
        A NetworkX graph object.
    """
    if not isinstance(nxGraph, networkx.classes.graph.Graph):
        raise NotImplementedError()
    
    # add label to edges
    edges = nxGraph.edges(keys = True)

    attributeDict = dict()
    for edge in edges:
        try:
            name = edge[2].name
        except AttributeError:
            name = None
        attributeDict[edge] = name if name is not None else edge[2].__str__()
    
    networkx.set_edge_attributes(nxGraph, attributeDict, LABEL_NAME)
    
    # add label to nodes
    nodes = nxGraph.nodes
    
    attributeDict = dict()
    for node in nodes:
        try:
            name = node.name
        except AttributeError:
            name = None
        attributeDict[node] = name if name is not None else node.__str__()
    
    networkx.set_node_attributes(nxGraph, attributeDict, LABEL_NAME)

def addIdAttribute(nxGraph: networkx.classes.MultiGraph):
    """
    Adds the "custom_id" attribute to each node and edge.
    
    Add an attribute to nodes and edges called "custom_id" (see module variable :attr:`ID_NAME`) containing the string of the node's id, or the edge's key's id, repectively.
    
    Parameters
    ----------
    nxGraph : networkx.classes.MultiGraph
        A NetworkX graph object.
    """
    if not isinstance(nxGraph, networkx.classes.graph.Graph):
        raise NotImplementedError()
    
    # add label to edges
    edges = nxGraph.edges(keys = True)

    attributeDict = dict()
    for edge in edges:
        attributeDict[edge] = edge[2].__str__()
    
    networkx.set_edge_attributes(nxGraph, attributeDict, ID_NAME)
    
    # add label to nodes
    nodes = nxGraph.nodes
    
    attributeDict = dict()
    for node in nodes:
        attributeDict[node] = node.__str__()
    
    networkx.set_node_attributes(nxGraph, attributeDict, ID_NAME)

def addDescriptionAttribute(nxGraph: networkx.classes.MultiGraph):
    """
    Adds the "custom_description" attribute to each node and edge.
    
    Add an attribute to nodes and edges called "custom_description" (see module variable :attr:`DESCRIPTION_NAME`) containing the node's or edge's `.description` field, if there is any.
    
    Parameters
    ----------
    nxGraph : networkx.classes.MultiGraph
        A NetworkX graph object.
    """
    if not isinstance(nxGraph, networkx.classes.graph.Graph):
        raise NotImplementedError()
    
    # add description to edges
    edges = nxGraph.edges(keys = True)

    attributeDict = dict()
    for edge in edges:
        try:
            attributeDict[edge] = edge[2].description
        except AttributeError:
            continue
    
    networkx.set_edge_attributes(nxGraph, attributeDict, DESCRIPTION_NAME)
    
    # add description to nodes
    nodes = nxGraph.nodes
    
    attributeDict = dict()
    for node in nodes:
        try:
            attributeDict[node] = node.description
        except AttributeError:
            continue
    
    networkx.set_node_attributes(nxGraph, attributeDict, DESCRIPTION_NAME)
    
def addReactionAttribute(nxGraph: networkx.classes.MultiGraph):
    """
    Adds the "custom_reaction" attribute to each edge.
    
    Add an attribute to edges called "custom_reaction" (see module variable :attr:`REACTION_NAME`) containing the edge's `.reaction` field, if there is any.
    
    Parameters
    ----------
    nxGraph : networkx.classes.MultiGraph
        A NetworkX graph object.
    """
    if not isinstance(nxGraph, networkx.classes.graph.Graph):
        raise NotImplementedError()
    
    # add reaction to edges
    edges = nxGraph.edges(keys = True)

    attributeDict = dict()
    for edge in edges:
        try:
            attributeDict[edge] = edge[2].reaction
        except AttributeError:
            continue
    
    networkx.set_edge_attributes(nxGraph, attributeDict, REACTION_NAME)

def addMajorityAttribute(graph: Models.CommonGraphApi, totalNumberOfOrganisms: int):
    """
    Adds the "custom_majority" attribute to each node and edge.
    
    Add an attribute to nodes and edges called "custom_majority" (see module variable :attr:`MAJORITY_NAME`).
    It contains each node's or edge's majority percentage, if the graph contains counts (`graph.nodeCounts` or `graph.edgeCounts`, see :class:`FEV_KEGG.Graph.Models.CommonGraphApi`).
    For example, say edge A were to occur in 52 of all organisms used to create the `graph`, then ```graph.edgeCounts[A] == 52```.
    Now we have also ```totalNumberOfOrganisms == 76```. Then, we will write ```ceil(52/76*100) = 68``` into the "custom_majority" attribute.
    
    Parameters
    ----------
    graph : Models.CommonGraphApi
        A graph object.
    totalNumberOfOrganisms : int
        Total number of organisms which were involved in creating the `graph`. This is used to calculate the percentage of counts. 100% == `totalNumberOfOrganisms`.
    """
    nxGraph = graph.underlyingRawGraph
    if not isinstance(nxGraph, networkx.classes.graph.Graph):
        raise NotImplementedError("This graph model can not be annotated for majority, yet.")
    
    # add color to edges
    if graph.edgeCounts is not None:
        
        attributeDict = dict()
        for edge in graph.getEdges():
            attributeDict[edge] = math.ceil(graph.edgeCounts.get(edge, 0) / totalNumberOfOrganisms * 100)
        
        networkx.set_edge_attributes(nxGraph, attributeDict, MAJORITY_NAME)
    
    # add color to nodes
    if graph.nodeCounts is not None:
        
        attributeDict = dict()
        for node in graph.getNodes():
            attributeDict[node] = math.ceil(graph.nodeCounts.get(node, 0) / totalNumberOfOrganisms * 100)
        
        networkx.set_node_attributes(nxGraph, attributeDict, MAJORITY_NAME)

def addUrlAttribute(nxGraph: networkx.classes.MultiGraph):
    """
    Adds the "URL" attribute to each edge.
    
    Add an attribute to edges called "URL" (see module variable :attr:`URL_NAME`) containing the edge's or node's URL to KEGG.
    
    Parameters
    ----------
    nxGraph : networkx.classes.MultiGraph
        A NetworkX graph object.
    """
    if not isinstance(nxGraph, networkx.classes.graph.Graph):
        raise NotImplementedError()
    
    # add label to edges
    edges = nxGraph.edges(keys = True)

    attributeDict = dict()
    for edge in edges:
        attributeDict[edge] = edge[2].getUrl()
    
    networkx.set_edge_attributes(nxGraph, attributeDict, URL_NAME)
    
    # add label to nodes
    nodes = nxGraph.nodes
    
    attributeDict = dict()
    for node in nodes:
        attributeDict[node] = node.getUrl()
    
    networkx.set_node_attributes(nxGraph, attributeDict, URL_NAME)


class Colour(Enum):
    BLUE = '#4444FF'
    RED = '#FF5555'
    PINK = '#FF66FF'
    GREEN = '#55FF55'
    YELLOW = '#FFFF55'
    TURQUOISE = '#55FFFF'

def addColourAttribute(graph: Models.CommonGraphApi, colour : Colour, nodes = False, edges = False):
    """
    Adds the "colour" attribute to selected nodes/edges.
    
    If both, `nodes` and `edges` are *False*, nothing is coloured.
    Adds an attribute to nodes and edges called "colour" (see module variable :attr:`COLOUR_NAME`) containing the string of a hex value of a colour in RGB, see :class:`Colour`.
    
    Parameters
    ----------
    graph : Models.CommonGraphApi
        A graph object.
    colour : Colour
        Colour to use.
    nodes : Iterable[NodeView] or bool, optional
        Nodes to be coloured, i.e. annotated with the "colour" attribute containing `colour`. If *True*, all nodes are coloured.
    edges : Iterable[EdgeView] or bool, optional
        Edges to be coloured, i.e. annotated with the "colour" attribute containing `colour`. If *True*, all edges are coloured.
    
    Warnings
    ----
    Any operation on the resulting graph, e.g. :func:`FEV_KEGG.Graph.Models.CommonGraphApi.intersection`, removes the colouring. It actually removes all attributes.
    """
    nxGraph = graph.underlyingRawGraph
    if not isinstance(nxGraph, networkx.classes.graph.Graph):
        raise NotImplementedError("This graph model can not be coloured, yet.")
    
    # add color to edges
    if edges is not False:
        # colour something
        if edges is True:
            # colour everything
            edgesToColour = nxGraph.edges
        else:
            # colour what is in `edges`
            edgesToColour = edges
        
        attributeDict = dict()
        for edge in edgesToColour:
            attributeDict[edge] = colour.value
        
        networkx.set_edge_attributes(nxGraph, attributeDict, COLOUR_NAME)
    
    # add color to nodes
    if nodes is not False:
        # colour something
        if nodes is True:
            # colour everything
            nodesToColour = nxGraph.nodes
        else:
            # colour what is in `nodes`
            nodesToColour = nodes
        
        attributeDict = dict()
        for node in nodesToColour:
            attributeDict[node] = colour.value
        
        networkx.set_node_attributes(nxGraph, attributeDict, COLOUR_NAME)


def toGraphML(graph: Models.CommonGraphApi, file, inCacheFolder = False):
    """
    Export `graph` to `file` in GraphML format.
    
    Each of the graph element's attributes is translated into a column.
    
    Parameters
    ----------
    graph : Models.CommonGraphApi
        The graph to be exported.
    file : str
        Path and name of the exported file. See `inCacheFolder`.
    inCacheFolder : bool, optional
        If *True*, interpret `file` relative to the cache folder. See :attr:`FEV_KEGG.settings.cachePath`.
        If *False*, interpret `file` relative to the current working directory.
        
    Raises
    ------
    NotImplementedError
        If `graph` is not of a NetworkX type.
    """
    nxGraph = graph.underlyingRawGraph
    if isinstance(nxGraph, networkx.classes.graph.Graph):
        
        if inCacheFolder is True:
            file = os.path.join(settings.cachePath, file)
        
        dirName = os.path.dirname(file)
        if not os.path.isdir(dirName) and dirName != '':
            os.makedirs(os.path.dirname(file))
        
        networkx.write_graphml(nxGraph, file + '.graphml', prettyprint=False)
        
    else:
        raise NotImplementedError()


def toGML(graph: Models.CommonGraphApi, file, inCacheFolder = False):
    """
    Export `graph` to `file` in GML format.
    
    Each of the graph element's attributes is translated into a column.
    
    Parameters
    ----------
    graph : Models.CommonGraphApi
        The graph to be exported.
    file : str
        Path and name of the exported file. See `inCacheFolder`.
    inCacheFolder : bool, optional
        If *True*, interpret `file` relative to the cache folder. See :attr:`FEV_KEGG.settings.cachePath`.
        If *False*, interpret `file` relative to the current working directory.
     
    Raises
    ------
    NotImplementedError
        If `graph` is not of a NetworkX type.
    """
    nxGraph = graph.underlyingRawGraph
    if isinstance(nxGraph, networkx.classes.graph.Graph):
        
        if inCacheFolder is True:
            file = os.path.join(settings.cachePath, file)
        
        dirName = os.path.dirname(file)
        if not os.path.isdir(dirName) and dirName != '':
            os.makedirs(os.path.dirname(file))
        
        networkx.write_gml(nxGraph, file + '.gml', lambda x: x.__str__())
        
    else:
        raise NotImplementedError()
    
def forCytoscape(graph: Models.CommonGraphApi, file, inCacheFolder = False, addDescriptions = True, totalNumberOfOrganisms: int = None):
    """
    Export `graph` to `file` in GraphML format, including some tweaks for Cytoscape.
    
    Parameters
    ----------
    graph : Models.CommonGraphApi
        The graph to be exported.
    file : str
        Path and name of the exported file. See `inCacheFolder`.
    inCacheFolder : bool, optional
        If *True*, interpret `file` relative to the cache folder. See :attr:`FEV_KEGG.settings.cachePath`.
        If *False*, interpret `file` relative to the current working directory.
    addDescriptions : bool, optional
        If *True*, downloads additional descriptions for substance
    totalNumberOfOrganisms : int, optional
        Total number of organisms which were involved in creating the `graph`. This is used to calculate the percentage of counts. 100% == `totalNumberOfOrganisms`.
        If *None*, no majority attribute is added. Only relevant if `graph` has counts, which it usually does not.
        See :func:`addMajorityAttribute` for more info.
            
    Raises
    ------
    NotImplementedError
        If `graph` is not of a NetworkX type.
    
    Warnings
    --------
    Adding extra descriptions, when `addDescriptions` == *True* (default!) is a lengthy process and can take some minutes.
    """
    # fetch extra descriptions, so they can be saved in attributes
    if addDescriptions is True:
        graph.addSubstanceDescriptions()
        if isinstance(graph, SubstanceEcGraph):
            graph.addEcDescriptions()
    
    # add certain attributes to the graph, so they can be stored in the resulting XML
    addLabelAttribute(graph.underlyingRawGraph)
    addIdAttribute(graph.underlyingRawGraph)
    addDescriptionAttribute(graph.underlyingRawGraph)
    addUrlAttribute(graph.underlyingRawGraph)
    if addDescriptions is True and isinstance(graph, SubstanceEcGraph):
        addReactionAttribute(graph.underlyingRawGraph)
    
    if totalNumberOfOrganisms is not None:
        addMajorityAttribute(graph, totalNumberOfOrganisms)
    
    toGraphML(graph, file, inCacheFolder=inCacheFolder)
    