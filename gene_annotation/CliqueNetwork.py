from AbstractNetwork import AbstractNetwork
from Node import Node
from random import import choice
import csv

class CliqueNetwork(AbstractNetwork):
    def __init__(self, filename):
        """
        Create a network from a file
        """
        self.nodes = {}
        self.cliques = {}
        
        with open(filename, 'r') as input:
            reader = csv.reader(input, delimiter='\t')
            for row in reader:
                node1 = self.getNode(row[0])
                node2 = self.getNode(row[1])
                if not node1.hasLinkTo(node2):
                    node1.addLinkTo(node2)
                    node2.addLinkTo(node1)
                    self.addCliqueOfSize(Clique([node1, node2]), 2)
    
    
    def findCliquesUpToN(self, n):
        for i in range(3, n+1):
            self.extendCliques()
            
        
        for size in self.cliques:
            if size > 2:
                for clique in self.cliques[size]:
                    if not clique.contained_in_larger_clique:
                        print(clique)
                
    
    
    """
    Try to extend the cliques with the highest n by 1
    """
    def extendCliques(self):
        previous_n = max(list(self.cliques))
        biggest_cliques = self.cliques[previous_n]
        for clique in biggest_cliques:
            extended_cliques = clique.extend()
            self.addCliquesOfSize(extended_cliques, previous_n + 1)
            
    
    def addCliquesOfSize(self, cliques, size):
        for clique in cliques:
            
            if not size in self.cliques:
                self.addCliqueOfSize(clique, size)
            elif not clique in self.cliques[size]:
                self.addCliqueOfSize(clique, size)
    
    
    def addCliqueOfSize(self, clique, size):
        if size in self.cliques:
            self.cliques[size].append(clique)
        else:
            self.cliques[size] = [clique]
    
    def evolveNetwork(self, t):
        for t in range(0, t):
            """
            choosing edge to delete
            """
            to_delete = choice(self.cliques[2])
            self.addRandomEdge()
            
            self.cliques[2].remove(to_delete)
            to_delete.nodes[0].removeNode(to_delete.nodes[1])
            to_delete.nodes[1].removeNode(to_delete.nodes[0])
    
    def addRandomEdge(self):
        node1 = choice(self.nodes)
        node2 = choice(self.nodes)
        if node1 == node2 or node1.hasLinkTo(node2):
            self.addRandomEdge()
        else:
            node1.addLinkTo(node2)
            node2.addLinkTo(node1)


class Clique():
    def __init__(self, nodes):
        self.nodes = nodes
        self.contained_in_larger_clique = False
    
    """
    Try to extend clique by 1 and returns all possible extensions
    """
    def extend(self):
        for i in range(0, len(self.nodes)):
            """
            Determine all nodes which have an edge to all nodes in the current clique
            """
            if i == 0:
                intersection = self.nodes[i].nodelist
            else:
                intersection= [node for node in self.nodes[i].nodelist if node in intersection]
        
        new_cliques = []
        if intersection:
            self.contained_in_larger_clique = True
            for node in intersection:
                tmp=self.nodes.copy()
                tmp.append(node)
                
                new_cliques.append(Clique(tmp))
            
        return(new_cliques)
    
    def __eq__(self, other):
        if isinstance(self, other.__class__):
            return frozenset(self.nodes) == frozenset(other.nodes)
        else:
            return False
        
    def __hash__(self):
        return hash(frozenset(self.nodes))
    
    
    def __str__(self):
        return(", ".join(str(x) for x in self.nodes))
        