#    Copyright (C) 2013 by
#    David Combe (david.combe@gmail.com), Universite Jean Monnet, Saint-Etienne
#    All rights reserved.
#    BSD license. 


from pprint import pprint 
import networkx as nx # pour les exemples
import sys
import random

from scipy.spatial.distance import cosine, euclidean
from scipy import zeros
import numpy as np
from scipy.spatial.distance import cosine, euclidean



# Les fonctions a produire sont
# - variInter( node, dest)
# - (variTotale)
# - hypothese de deplacement: moveNode(self, node, comSource, comTarget)
# - deplacement effectif: reallyMoveNode(self, node, comSource, comTarget)



#global globTfIdfTab
#globTfIdfTab=[]


global node2com 
node2com = {}


def graphExample(): 
    global globTfIdfTab
    
    globTfIdfTab=[
                            [2,1],
                            [3,1],
                            [22,21],
                            [20,22]#,
                            #[21,22] 
                    ]
    
    global globAuthorIndex
    globAuthorIndex={"alpha": 0,
                     "beta": 1,
                     "delta": 2,
                     "gamma": 3
                     }  


    G=nx.Graph()
    
    G.add_node("alpha") 
    G.add_node("beta")
    G.add_node("delta")
    G.add_node("gamma")
    #G.add_node(4)

    G.add_edge("alpha","beta")
    G.add_edge("beta","delta")
    G.add_edge("delta","gamma")
    G.add_edge("alpha","delta")
    return G
        
        
def graphBiggerExample(): 

    global globTfIdfTab
    
    globTfIdfTab=[
                            [2,1],
                            [3,1],
                            [2,2],
                            [1,1],
                            [22,21],
                            [20,22],
                            [1,2]    
                    ]
        
    global globAuthorIndex
    globAuthorIndex={0: 0,
                     1: 1,
                     2: 2,
                     3: 3,
                     4: 4,
                     5: 5,
                     6: 6
                     }  
    
    G=nx.Graph()
	
    G.add_node(0) 
    G.add_node(1)
    G.add_node(2)
    G.add_node(3)
    G.add_node(4)
    G.add_node(5)
    G.add_node(6)
    
    G.add_edge(0,1)
    G.add_edge(1,2)
    G.add_edge(1,3)
    G.add_edge(2,3)
    G.add_edge(3,4)
    G.add_edge(4,5)
    G.add_edge(4,6)
    G.add_edge(6,2)
    G.add_edge(5,0)
    
    return G
        
        
        

def getCom(node):
    global node2com 
    return node2com[node]

def setCom(node,com):
    global node2com 
    node2com[node]=com
    




class APartition :
    node2com={}
    
    def __init__(self, partition) :
        self.node2com = partition
        
        
        
    def dicInverse(self):
        dicinv={}
        for (k,v) in self.node2com.items(): 
                a=dicinv.get(v, set())
                a.add(k)
                dicinv[v]=a
                #print "k:"+str(k)+"  v:"+str(v)
        #print "dic_inv"
        #pprint(dicinv)
        return dicinv
        
        
    def __hash__(self):
        dicinv=self.dicInverse()
        res=0
        for k in dicinv:
            res=res+(','.join(dicinv[k])).__hash__()
        return res
        
        
    def __eq__(self, b):
        """
        Deux partitions sont identiques si...
        """
        commonDic={}
        for i in self.node2com:
            if (self.node2com[i] not in commonDic):
                commonDic[self.node2com[i]]=b.node2com[i]
            else:
                if (commonDic[self.node2com[i]]!=b.node2com[i]):
                    return False
        if len(set(self.node2com.values())) != len(set(b.node2com.values())):
            return False
        return True
    
    def __str__(self):
        return str(self.node2com)
    
    
    
    
    
class VarianceStatus :
    """
    State of variance
    """
    
    globalTfIdfTab=[]
    centreDeGraviteTotal=[]
    currentVarInter=-1
    
    # pour chaque communaute, le sous-compte de variance inter
    #varInterSubCount = {}
    varInterSubGravityCenter = {}
    
    currentPartition=None
    
    globAuthorIndex=None
    
    classesNumber = 0
    
    def __init__(self, graph, globalTfIdfTab, partition, globAuthorIndex ) :
        """
        graph, globalTfIdfTab, partition, globAuthorIndex
        """
        self.globAuthorIndex = globAuthorIndex
        self.varInterSubCount = dict([])
        self.globalTfIdfTab=globalTfIdfTab
        self.centreDeGraviteTotal=self.moyenne_points(range(len(globalTfIdfTab)))
        
        print "init: cdg total"
        pprint(self.centreDeGraviteTotal)
        
        sum={}
        count={}
        
        #try: #la partition est un Status de Aynaud
        for k,v in partition.node2com.items():
            count[v]=count.get(v, 0)+1
            a=sum.get(v,  zeros([len(self.globalTfIdfTab[0])])  )
            b=self.globalTfIdfTab[self.globAuthorIndex[k]]
            sum[v]= np.sum([a,b]  , axis=0   )
            self.currentPartition = partition
            
        
                
                
        #calcul des sous-centres de gravite
        for k,v in sum.items():
            #print str(k) + "<->" + str(v)
            #print "count: "+str(count[k])
            self.varInterSubGravityCenter[k]=v/count[k]
        
        #calcul couteux total de la variance inter
        self.currentVarInter=self.varianceInter()
        
        self.classesNumber = len(globalTfIdfTab)
        
        print "VarIncr structure initialised"
        #for i in self.currentPartition.dicInverse():
        #    self.varInterSubCount[i]=self.__partialVarianceInterOfClass( i)
        
        
        
    def __partialVarianceInterOfClass(self, cla):
        """
        Carre de la distance
        """
        #print "cdg partiel"
        #pprint(self.varInterSubGravityCenter[cla])
        #print "cdg total"
        #pprint(self.centreDeGraviteTotal)
        #pprint(self.varInterSubGravityCenter)
        #pprint(cla)
        #pprint(self.varInterSubGravityCenter[cla])
        #pprint(self.centreDeGraviteTotal)
        
        
        dis=euclidean(self.varInterSubGravityCenter[cla],self.centreDeGraviteTotal)
        #print "distance "+ str(dis)
        #print "coef: "+str(len(self.currentPartition.dicInverse()[cla]))+"/"+str(len(self.globalTfIdfTab))+"  -----   "+str((len(self.currentPartition.dicInverse()[cla])/float(len(self.globalTfIdfTab))))
        #print (len(self.currentPartition.dicInverse()[cla])/float(len(self.globalTfIdfTab)))*(dis)**2
        return dis**2
        #dic=currentPartition.dicInverse()
        #listIndices=dic[cla]
    
    
    def varianceInter(self):
        """
        Calcule la variance inter de la partition courante
        """
        sum=0
        for i in self.currentPartition.dicInverse().keys():
            sum=sum+ ( (len(self.currentPartition.dicInverse()[i]))   *   self.__partialVarianceInterOfClass(i))
        return sum / float(len(self.globalTfIdfTab))
        
    
    def setCentreDeGraviteTotal(cdg):
        self.centreDeGraviteTotal=cdg
    
    
    def calculatePartOfVarInterWithGivenGC(self, gc ):
        res=euclidean(gc, self.centreDeGraviteTotal)
        return res**2
        
        
    def reallyMoveNode(self, node, comSource, comTarget   ):
        """
        Met a jour la partition, les centres de gravite, la var inter courante.
        """
        assert(self.currentPartition.node2com[node] == comSource)
        
        (varInter, cdgSource, cdgCible, classesNumber) = self.moveNode(node, comSource, comTarget)
        
        self.currentVarInter=varInter
        
        if (cdgSource!=None): 
            #classesNumber = self.classesNumber
            pass
        else: #si une communaute a disparu
            #classesNumber = self.classesNumber -1
            self.classesNumber = self.classesNumber - 1
        
        # pour chaque communaute, le sous-compte de variance inter
        
        self.varInterSubGravityCenter[comSource] = cdgSource
        if (cdgSource == None):
            del self.varInterSubGravityCenter[comSource]
        self.varInterSubGravityCenter[comTarget] = cdgCible
        
        self.currentPartition.node2com[node]=comTarget
        
        assert(self.currentPartition.node2com[node] == comTarget)
        return self.currentVarInter

    def moveNode(self, node, comSource, comTarget   ): #, varInterSubCountComSource,   varInterSubCountComTarget    ):
        """
        On imagine que l'on deplace un noeud et on renvoie la variance qu'il y aurait dans cette configuration
        On renvoie la varInter qu'il y aurait dans cette configuration, le cdg source et le cdg cible
        """
        #on enleve le carre de la distance avec les anciens centres de gravite qu'on soustrait
        if (comSource == comTarget ):
            return self.currentVarInter,self.varInterSubGravityCenter[comSource] ,self.varInterSubGravityCenter[comTarget] ,self.classesNumber
        assert(self.currentPartition.node2com[node]==comSource)
        
        
        # on refait un calcul de style centre de gravite:
        # on prend la variance inter courante (une moyenne)
        res=self.currentVarInter
        #print "Depart "+str(res)

        #print "moveNodeVariance moins"
        # on enleve la variance inter d'avant de la communaute source * truc/n
        res=res-  self.__partialVarianceInterOfClass( comSource)* (len(self.currentPartition.dicInverse()[comSource]))   /len(self.globalTfIdfTab)   
        #print res
        
        # on enleve la variance inter d'avant de la communaute cible * truc2/n
        res=res-  self.__partialVarianceInterOfClass( comTarget)* (len(self.currentPartition.dicInverse()[comTarget]))   /len(self.globalTfIdfTab)   
        #print res
        #print "3 __partialVarianceInterOfClass( comTarget)"
        #print self.__partialVarianceInterOfClass( comTarget)* (len(self.currentPartition.dicInverse()[comTarget]))   /len(self.globalTfIdfTab)   
        
        
        #print "devrait etre zero: "+str(res)
        #exit()
        
        
        # calcul des 2 nouveaux gc
        gc1dest, gc2source = self.moveGravityCenters( self.varInterSubGravityCenter[comTarget]      ,      self.globalTfIdfTab[self.globAuthorIndex[node]]    , (len(self.currentPartition.dicInverse()[comTarget]))  , self.varInterSubGravityCenter[comSource], (len(self.currentPartition.dicInverse()[comSource])))
        if (gc2source!=None): 
            classesNumber = self.classesNumber
            
        else: #si une communaute a disparu
            classesNumber = self.classesNumber -1
        
        #print "moveNodeVariance plus"
        # on ajoute la variance inter de la communaute source * truc-1/n
        if (gc2source!=None): # si gc2source n'a pas disparu
            res=res+     self.calculatePartOfVarInterWithGivenGC(gc2source)  * (len(self.currentPartition.dicInverse()[comSource])-1) / float(len(self.globalTfIdfTab))
        #print res
        # on ajoute la variance inter de la communaute cible * truc2+1/n
        res=res+     self.calculatePartOfVarInterWithGivenGC(gc1dest)  *   (len(self.currentPartition.dicInverse()[comTarget])+1) / float(len(self.globalTfIdfTab))
        #print res
        
        
        
        
        
        
        #on prend les nouveau centre et on calcule les nouveaux carres qu'on ajoute            
        
        assert(self.currentPartition.node2com[node]==comSource)
        
        return res, gc2source, gc1dest, classesNumber
        
        
    def variance(self):
        global node2com
        #pour chaque communaute
        print node2com.values()
        
        
        return 
    
    def __str__(self) :
        out="Partition:\n"
        out=out + str(self.currentPartition)
        out=out + "\nSous-compte de la variance inter:\n"
        #for (k,v) in self.varInterSubCount.items():
        #    out=out+str(k)+","+str(v)
        for i in self.currentPartition.dicInverse():
            out=out+str(i)+" - "+str(self.__partialVarianceInterOfClass(i))+"\n"
        out = out + "\nVariance inter: "+str(self.varianceInter())+"\n"
        
        #print "\nVariance totale: "+str(self.)
        
        return out
        
       


    def moveGravityCenters(self,initialCom1GravityCenter, newPoint, initialAmountOfNodesInCom1, initialCom2GravityCenter, initialAmountOfNodesInCom2):
        """
        On donne le nouveau point, le GC1 de destination et le poids de la com 1 et on retourne le nouveau centre de gravite
        Idem avec la communaute source
        """
        newGC1=zeros([len(initialCom1GravityCenter)])
        newGC2=zeros([len(initialCom1GravityCenter)])
        
        for i in range(len(initialCom1GravityCenter)):
            newGC1[i]=(initialAmountOfNodesInCom1*initialCom1GravityCenter[i]+newPoint[i])/(initialAmountOfNodesInCom1+1)
            assert(newGC1[i]!=float("inf"))
            assert(newGC1[i]!=float("-inf"))
            #si la com2 n'existe plus:
            if((initialAmountOfNodesInCom2-1)!=0):
                newGC2[i]=((initialAmountOfNodesInCom2)*initialCom2GravityCenter[i]-newPoint[i])/(initialAmountOfNodesInCom2-1)
                assert(newGC2[i]!=float("inf"))
                assert(newGC2[i]!=float("-inf"))
            else:
                newGC2=None
        return (newGC1, newGC2)


    def moyenne_points(self,indices):
        '''
        Calcule le centre de gravite d'un ensemble de points,
        a savoir la moyenne des differents attributs
        Les vecteurs sont recuperes dans la variable globale globTfIdfTab
        '''
        nb_att=len(self.globalTfIdfTab[indices[0]])
        tab=zeros([nb_att])
        
        #pprint(tab)
        #pprint(indices)
        
        for i in indices:
            for j in range(nb_att):
                tab[j]=tab[j]+self.globalTfIdfTab[i][j]
        for j in range(nb_att):
            tab[j]=tab[j]/len(indices)
        return tab
    
    def induced_varStatus(self):
        """
        On remodele le varstatus selon les communautes.
        
        """
        pprint(self)
        print "Not implemented"
        exit()
        #partition discrete
        currentPartition = APartition()
        globalTfIdfTab=[]
        centreDeGraviteTotal=[]
        currentVarInter=-1
        
        # pour chaque communaute, le sous-compte de variance inter
        #varInterSubCount = {}
        varInterSubGravityCenter = {}
        
        currentPartition=None
        
        globAuthorIndex=None
        
        #classesNumber = 0 Ne change pas
        
        
        
    
    
    

#def testMoveGravityCenter():
#    global globTfIdfTab
#    G=graphExample()
#    print G.nodes()
#    
#    global node2com
#    node2com["alpha"]= 0
#    node2com["beta"] = 0
#    node2com["delta"]= 1
#    node2com["gamma"]= 1
#    
#    vs=VarianceStatus(None, globTfIdfTab, partition)
#    
#    gc1=self.moyenne_points([2,3])
#    gc2=self.moyenne_points([0,1])
#    pprint(gc1)
#    pprint(gc2)
#    res1,res2= (self.moveGravityCenters(gc1, globTfIdfTab[1], 2, gc2, 2))
#    pprint(res2)
#    assert (res2==[2., 1.]).all()
#    print "Test de bougeage de cdg passe avec succes"



def testMoveNode():
    part=APartition({"alpha":0,"beta":0,"delta":0,"gamma":1})
    global globTfIdfTab
    G=graphExample()
    vs = VarianceStatus(G, globTfIdfTab, part)
    
    print(vs.__str__())
    
    #print "move node"
    var, rien1, rien2 = vs.moveNode( "gamma", 1, 0)#, varInterSubCountComSource,   varInterSubCountComTarget    )
    #print "var inter avec noeud deplace: "+str(var)
    # la var inter de la partition unitaire est zero
    assert(var<0.01)
    assert(var>-0.01)
    print "Test de bougeage de noeud passe avec succes"


def __main():
    part=APartition({"alpha":0,"beta":0,"delta":1,"gamma":1})
    #part.init({"a":0,"b":0,"c":0,"d":1})
    
    part2=APartition({"beta":1,"delta":0,"gamma":0,"alpha":1})
    #part2.init({"a":0,"b":0,"c":0,"d":2})
    
    assert(part.__hash__()== part2.__hash__())


    print part.__eq__(part)

    #testMoveGravityCenter()
    
    testMoveNode()
    
    
    
    
    
    global globTfIdfTab
    G=graphBiggerExample()
    part=APartition({0: 0,
                     1: 0,
                     2: 5,
                     3: 5,
                     4: 6,
                     5: 6,
                     6: 6
                     }  )
    vs = VarianceStatus(G, globTfIdfTab, part)
    
    
    
    
    print vs.reallyMoveNode(0, 0, 6)
    print vs.reallyMoveNode(1, 0, 6)
    print vs.reallyMoveNode(2, 5, 6)
    print vs.reallyMoveNode(3, 5, 6)
    
    print "-----------fin-----------"
    print str(vs)
    
    #print vs.reallyMoveNode(3, 5, 6)
    
    
    
    part=APartition({0: 5,
                     1: 5,
                     2: 5,
                     3: 5,
                     4: 6,
                     5: 6,
                     6: 5
                     }  )
    vs = VarianceStatus(G, globTfIdfTab, part)
    
    print str(vs)
    
    
    
    
    
    
    
    
    exit()
    
    global globTfIdfTab
    G=graphExample()
    print G.nodes()
        
    print globTfIdfTab
    


if __name__ == "__main__" :
    __main()




