from AbstractNetwork import AbstractNetwork
from Node import Node
from collections import defaultdict
from collections import Counter
import csv
import math

class InteractionNetwork(AbstractNetwork):
    def __init__(self, filename):
        """
        Create a network from a file
        """
        self.nodes = {}
        with open(filename, 'r') as input:
            reader = csv.reader(input, delimiter='\t')
            self.number_of_pairs = 0
            for row in reader:
                self.number_of_pairs += 1
                node1 = self.getNode(row[0])
                node2 = self.getNode(row[1])
                if not node1.hasLinkTo(node2):
                    node1.addLinkTo(node2)
                    node2.addLinkTo(node1)
PPI = InteractionNetwork('human_network.tsv')
 
class UniProtMapping(AbstractNetwork):
    def __init__(self, filename):   
        self.UniProtID = []
        self.AlternativeNames = []
        with open(filename, 'r') as input:
            reader = csv.reader(input, delimiter='\t')
            for row in reader:
                self.UniProtID.append(row[0])
                self.AlternativeNames.append(row[4])
        UniProtMap = dict(zip(self.UniProtID, self.AlternativeNames))
        #print('UniProtMap:',len(UniProtMap))
        
        PPI = InteractionNetwork('pig_network.tsv')
        self.name_list = defaultdict(list)
        for node in PPI.nodes:
            for k, v in UniProtMap.items():
                if node in v or node in k:
                    self.name_list[node].append(k)
        print ('name list',len(self.name_list))
        
UniProt = UniProtMapping('human_uniprot.tsv')


class GOID(AbstractNetwork):    
    def __init__(self, filename):       
        #PART A: GET GOANNOTATION BY PROTEIN ID
        self.GOAnnotation = defaultdict(list)
        with open(filename, 'r') as input:
            reader = csv.reader(input, delimiter='\t')
            for row in reader:
                if row[0] == "UniProtKB" and row[-9] == "P":
                    for node in UniProt.name_list:
                        for UniProtID in UniProt.name_list[node]:
                            if UniProtID == row[1]:
                                self.GOAnnotation[node].append(row[4])                        
                else:
                    pass
        print (self.GOAnnotation)
        print ('Annotation: ',len(self.GOAnnotation))

    def protein_by_annotation(self):
        self.protein_by_annotation = defaultdict(list)
        for ID in self.GOAnnotation.keys():
            print (ID)
            for value in self.GOAnnotation[ID]:
                self.protein_by_annotation[value].append(ID)        
        print ('Protein by annotation: ',self.protein_by_annotation)
    
        
        self.protein_per_annotation = {}
        count = 0
        for item in self.protein_by_annotation:
            self.protein_per_annotation[item] = len(self.protein_by_annotation[item])
            if self.protein_per_annotation[item] == 1:
                count +=1
        print ('Number of unique annotations in the network: ',count)
        print ('Highest no of proteins per annotation: ',max(self.protein_per_annotation.values()))
        print ('Smallest no of proteins per annotation: ',min(self.protein_per_annotation.values()))
        print ('Average no of proteins per annotation: ',sum(self.protein_per_annotation.values())/len(self.protein_per_annotation))
        

    def Overview(self):  
        #PART B: GENERATING OVERVIEW
        print ('Number of protein in the network: ', len(PPI.nodes))
        print ('Number of protein without annotation: ',len(PPI.nodes)-len(self.GOAnnotation))
        print ('Percentage of protein without annotation: ',(len(PPI.nodes)-len(self.GOAnnotation))/len(PPI.nodes)*100,'%')
        
        self.annotation_per_protein = {}
        for node in self.GOAnnotation:
            self.annotation_per_protein[node] = len(self.GOAnnotation[node])
        print ('Annotation per protein:',self.annotation_per_protein)
        print ('Highest no of annotations per protein: ',max(self.annotation_per_protein.values()))
        print ('Smallest no of annotations per protein: ',min(self.annotation_per_protein.values()))
        print ('Average no of annotations per protein: ',sum(self.annotation_per_protein.values())/len(self.annotation_per_protein))
        
        
    def Most_common_annotation(self,n):   
        #PART C: GET N MOST COMMON ANNOTATIONS:
        protein_per_annotation = self.protein_per_annotation
        print ('Most common 5 annotations:',Counter(protein_per_annotation).most_common(n))
        print ('Least common 5 annotations:',Counter(protein_per_annotation).most_common()[-n:])
                
    def AnnotationEnrichment(self,n):    
        protein_per_annotation = self.protein_per_annotation()
        protein_by_annotation = self.protein_by_annotation()
        PPI = InteractionNetwork('chicken_network1.tsv')
        #PART D: 
        #KA
        KA = defaultdict(dict)
        f = math.factorial
        for annot in protein_by_annotation:
            if protein_per_annotation[annot] >= 2:
                KA[annot] = f(protein_per_annotation[annot])//(f(2)*(f(protein_per_annotation[annot]-2)))
            else:
                KA[annot] = 0
        print ('KA: ', KA)
        #N
        N = f(len(PPI.nodes))//(f(2)*f(len(PPI.nodes)-2))
        print ('N: ', N)
        #n
        n = PPI.number_of_pairs
        print ('n: ',n)
        
        #ka
        ka = defaultdict(dict)
        count = 0
        for annot in protein_by_annotation:
            count = 0
            if KA[annot]>=1:
                for i in range(len(protein_by_annotation[annot])-1):
                    for j in range(i+1,len(protein_by_annotation[annot])):
                        if (protein_by_annotation[annot][i] in PPI.nodes) and (protein_by_annotation[annot][j] in PPI.nodes) and (PPI.nodes[protein_by_annotation[annot][i]].hasLinkTo(PPI.nodes[protein_by_annotation[annot][j]])):
                            count += 1
                        ka[annot] = count
            else: 
                ka[annot] = 0
        print ('ka: ',ka)
        
        #Calculate pA and report finding:
        pA = {}
        for annot in protein_by_annotation:
            K = KA[annot]
            i = ka[annot]
            prob = f(K)/(f(i)*f(K-i))*f(N-K)//(f(n-i)*f(N-K-n+i))/(f(N)/(f(N-n)*f(n)))
            pA[annot] = prob
        print (pA)
        
        count_prob1 = 0
        count_prob2 = 0 
        count_prob3 = 0
        for annot in pA:
            if pA[annot] < 0.05:
                count_prob1 += 1
            if 0.5 < pA[annot]:
                count_prob2 += 1
            if pA[annot] > 0.95:
                count_prob3 += 1
        print (count_prob1,count_prob2,count_prob3)
        print (count_prob1/len(protein_by_annotation)*100,'% ',count_prob2/len(protein_by_annotation)*100,'% ',count_prob3/len(protein_by_annotation)*100,'% ')    
     
        
        pA = sorted(pA.items(), key=lambda x: x[1])
        print('Lowest 5: ',pA[:n])
        print('Highest 5: ',pA[-n:])

a = GOID('human_GO.gaf')
print(a.Overview())
print(a.protein_by_annotation())
print(a.Most_common_annotation(5))