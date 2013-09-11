#    Copyright (C) 2013 by
#    David Combe (david.combe@gmail.com), Universite Jean Monnet, Saint-Etienne
#    All rights reserved.
#    BSD license. 

"""
Save the best partitions seen so far.

"""

from pprint import pprint

class State :
    
    mapping=None
    
    dicOfProperties = None
    
    
    def __init__(self, mapping, dicOfProperties):
        self.dicOfProperties = dicOfProperties
        self.mapping = mapping
    
    
    
    
class StateSaver :
    """
    Save different partitions. Output the n best.
    """
    
    
    numberOfPartitionsToSave=None
    
    groundTruth=None
    
    tab = None
    
    
    
    
    
    def __init__(self,numberOfPartitionsToSave=5, groundTruth=None) :
        self.numberOfPartitionsToSave=numberOfPartitionsToSave
        
        self.groundTruth=groundTruth
        self.tab = []
        
        
        
    def add(self, mapping, dicOfProperties):
        """
        Add a partition.
        """
        
        
        self.tab.insert(0, State(mapping, dicOfProperties))
        
        
        try:
            self.tab.pop(self.numberOfPartitionsToSave)
        except:
            pass
        if len(self.tab)>11:
            print "if len(self.tab)>11"
            exit()
    
    def saveOutput(self):
        """
        Exports the mappings of the n best partitions.
        """
        
        print "The best partition is at the beginning.\n"
        for i in self.tab :
            pprint(i.dicOfProperties)
            pprint(i.mapping)
        print "StatusSaver: END"
        
        
        
        
        
        
        











