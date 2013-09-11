#    Copyright (C) 2013 by
#    David Combe (david.combe@gmail.com), Universite Jean Monnet, Saint-Etienne
#    All rights reserved.
#    BSD license. 



"""
Convenience functions in order to change the representation of the partitions (dictionary - array).
"""



from scipy import zeros
from pprint import pprint

def createInvertedAuthorIndex(globAuthorIndex):
    """
    aut to indice toward authors array (i2aut)
    """
    
    
    resList = {}#zeros([len(globAuthorIndex)], str)
    
    for i in globAuthorIndex:
        #print str(i) + str(globAuthorIndex[i])
        resList[int(globAuthorIndex[i])] = i
        #print str(globAuthorIndex[i])+" "+str(i)+"\n"
    #pprint(resList)
    return resList

def autArrayToi2aut(globAuthorIndex):
    """
    aut to indice toward authors array (i2aut)
    """
    print "empty function !"
    exit()
    
    
    
def dicToArray(dic, globAuthorIndex, translator={"IR":3, "Agents":1, "DB":4, "ML":6, "HCI":5,"AI":2}):
    """
    Transform ground truth from dic to array using translator
    """
    print "dic glob translator"
    pprint(dic)
    pprint(globAuthorIndex)
    pprint(translator)
    res=zeros([len(globAuthorIndex)],int)
    for i in globAuthorIndex:
        res[globAuthorIndex[i]]=translator[dic[i]]
    return res

def dicToArraySimple(dic, globAuthorIndex):
    """
    Transform ground truth from dic to array
    """
    res=zeros([len(globAuthorIndex)],int)
    for i in globAuthorIndex:
        res[globAuthorIndex[i]]=dic[i]
    return res

def arrayToDic_b(partitionArray, invertedIndex):
    """
    Change partition notation array -> dic
    
    """
    dic={}
    for i in invertedIndex:
        dic[invertedIndex[i]]=partitionArray[i]
    #pprint(dic)
        
    return dic
    
    


def arrayToDic(partitionArray):
    """
    Change partition notation array -> dic
    
    """
    dic={}
    global inv
    
    #print "inv"
    #pprint(inv)
    
    for i in inv:
    
        #print str(i) 
        #print "inv "+inv[i]
        #print " -- "+str( partitionArray[i])
    
        dic[inv[i]]=partitionArray[i]
    #pprint(dic)
        
    return dic
    
# etirement (et renversement de la mesure de distance du texte en similarite)
def splineText(inputValueAsDistance, distanceOutput=False):
    if distanceOutput:
        return 1-splineText(inputValueAsDistance)

    #if(1-math.pow(inputValueAsDistance,30)==0):
    #     return 0.000001
    #else:
    #     return 1-math.pow(inputValueAsDistance,30)
    return 1-inputValueAsDistance
    
def arrayToDicIndex(authorIndex):
    authorIndexIndex={}
    i=0
    for aut in authorIndex:
        authorIndexIndex[aut]=i
        i=i+1
    return authorIndexIndex
    
