#    Copyright (C) 2013 by
#    David Combe (david.combe@gmail.com), Universite Jean Monnet, Saint-Etienne
#    All rights reserved.
#    BSD license. 


from communityThreshold import *
from math import *
from pprint import pprint


def testDistanceEucEntreCentroids():
    '''
        test de la distance euclidienne
    '''
    assert(distanceEucEntreCentroids([0,0,0], [0,0,0])==0)
    assert(distanceEucEntreCentroids([0,0,0], [0,0,1])==1)
    assert(distanceEucEntreCentroids([0,1,0], [0,0,0])==1)
    assert(distanceEucEntreCentroids([0,0,0], [1,1,1])==sqrt(3))
    assert(distanceEucEntreCentroids([0,0,0], [2,2,2])==2*sqrt(3))
    assert(distanceEucEntreCentroids([-1,-1,-1], [1,1,1])==2*sqrt(3))
    
def test_moyenne_points(): 
    '''
        test du calcul du centre de gravite
    '''
    globTfIdfTab=[
                    [-1,-1,1,1],
                    [1,1,-1,-1],
                    [1,-1,-1,1],
                    [-1,1,1,-1]
    ]
    setGlobTfIdfTab(globTfIdfTab)
    assert((moyenne_points([0,1,2,3])==[0.,0.,0.,0.]).all())
    assert((moyenne_points([0,1])==[0.,0.,0.,0.]).all())
    assert((moyenne_points([1,2])==[1.,0.,-1.,0.]).all())

def testVarianceInter():
    '''
        A IMPLEMENTER
    '''
    print "Variance non testee"
    pass

def testVarianceTotale():
    '''
        A IMPLEMENTER
    '''
    print "Variance non testee"
    pass

    
    
    
    
    
def graphExampleCorrele(): 
    global globTfIdfTab
    globTfIdfTab=[
                            [0],
                            [0],
                            [100],
                            [100]
                    ]
    
    global globAuthorIndex
    globAuthorIndex={"alpha": 0,
                     "beta": 1,
                     "gamma": 2,
                     "delta": 3
                     }
    G=nx.Graph()
    G.add_node("alpha") 
    G.add_node("beta")
    G.add_node("gamma")
    G.add_node("delta")
    
    G.add_edge("alpha","beta")
    G.add_edge("delta","gamma")
    return G
    
    
    
    
def graphExampleDecorrele(): 
    global globTfIdfTab
    globTfIdfTab=[
                            [0],
                            [0],
                            [100],
                            [100]
                    ]
    
    global globAuthorIndex
    globAuthorIndex={"alpha": 0,
                     "beta": 1,
                     "gamma": 2,
                     "delta": 3
                     }
    G=nx.Graph()
    
    G.add_node("alpha") 
    G.add_node("beta")
    G.add_node("gamma")
    G.add_node("delta")
    
    G.add_edge("alpha","gamma")
    G.add_edge("delta","beta")
    return G
    
    
    
    
    
    
    
def graphExample(): 
    global globTfIdfTab
    globTfIdfTab=[
                            [2,1],
                            [3,1],
                            [22,21],
                            [20,22]
                    ]
    
    global globAuthorIndex
    globAuthorIndex={"alpha": 0,
                     "beta": 1,
                     "gamma": 2,
                     "delta": 3
                     }
    G=nx.Graph()
    
    G.add_node("alpha") 
    G.add_node("beta")
    G.add_node("gamma")
    G.add_node("delta")
    
    G.add_edge("alpha","beta")
    G.add_edge("beta","gamma")
    G.add_edge("delta","gamma")
    G.add_edge("alpha","gamma")
    return G
    
    

def testMainAvecExemple():
    '''
    Teste la partition renvoyee avec le graphe exemple a 4 sommets (alpha beta / gamma delta)
    '''
    partition=testWithExample()
    pprint(partition)  
    assert(partition['alpha']==partition['beta'])
    assert(partition['gamma']==partition['delta'])
    assert(partition['alpha']!=partition['delta'])
    
def testWithBiggerExample():
    graph = graphBiggerExample()
    pprint(graph.nodes())
    partition = best_partition(graph, 8, None, None, None)
    print >> sys.stderr, str(modularity(partition, graph))
    for elem, part in partition.iteritems() :
        print str(elem) + " " + str(part)
        
    return partition
        
        
def testWithExample():
        #setDebug(0)
        
        
            
        graph=graphExampleCorrele()
        partition = best_partition(graph, 8, None, None, None)
        for elem, part in partition.iteritems() :
            print str(elem) + " " + str(part)
            
        graph=graphExampleDecorrele()
        partition = best_partition(graph, 8, None, None, None)
        for elem, part in partition.iteritems() :
            print str(elem) + " " + str(part)
            
            
            
        print "Petit exemple"
        graph = graphExample()
        partition = best_partition(graph, 8, None, None, None)
        print >> sys.stderr, str(modularity(partition, graph))
        for elem, part in partition.iteritems() :
            print str(elem) + " " + str(part)
            
        return partition
    
def testWithKarate():
    #ouverture de karate.bin
    #affectation de poids orthogonaux
    
    #assert(partition['alpha']==partition['beta'])
    #assert(partition['gamma']==partition['delta'])
    #assert(partition['alpha']!=partition['delta'])
    '''
        A IMPLEMENTER
    '''
    print "Karate non testee"
    # tester la meilleure modularite que l'on trouve avec le reseau karate
    # load_bin
    
    
def launchTests():
    print "Beginning tests"
    testDistanceEucEntreCentroids()
    test_moyenne_points()
    testVarianceInter()
    testVarianceTotale()
    testMainAvecExemple()
    testWithKarate()
    # tester avec karate
    print "Tests succeeded"
    
def __main() :
    """Main function in tests file"""
    launchTests()


if __name__ == "__main__" :
    __main()
