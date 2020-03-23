from Node import Node

class AbstractNetwork:
    """Abstract network definition, can not be instantiated"""
    
    def __init__(self, amount_nodes, amount_links):
        """
        Creates empty nodelist and call createNetwork of the extending class
        """
        self.nodes = {}
        self.__createNetwork__(amount_nodes, amount_links)

    def __createNetwork__(self, amount_nodes, amount_links):
        """
        Method overwritten by subclasses, nothing to do here
        """
        raise NotImplementedError

    def appendNode(self, node):
        """
        Appends node to network
        """
        self.nodes[node.id] = node

    def maxDegree(self):
        """
        Returns the maximum degree in this network
        """
        return max([x.degree() for x in self.nodes.values()])

    def size(self):
        """
        Returns network size
        """
        return len(self.nodes)

    def __str__(self):
        return "This network has %6d nodes." % (len(self.nodes))

    def getNode(self, identifier):
        """
        Returns node according to key
        """
        if identifier not in self.nodes:
            self.nodes[identifier] = Node(identifier)
            
        return self.nodes[identifier]
    