#    Copyright (C) 2013 by
#    David Combe (david.combe@gmail.com), Universite Jean Monnet, Saint-Etienne
#    Thomas Aynaud (thomas.aynaud@lip6.fr), 
#    All rights reserved.
#    BSD license. 


#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
This module implements community detection.
"""

from __future__ import division

import rpy2.robjects as robjects

from collections import deque
from sklearn import metrics
from varIncrementPhase2 import MetaVarIncrement, APartition, metaNodesHistory

from stateSaver import StateSaver
from partitionRepresentation import dicToArraySimple

import time as t

from pprint import pprint 
import networkx as nx
import sys
import random
import math

from scipy.spatial.distance import cosine, euclidean
from scipy import zeros

from scipy.stats import f

from gexf import write_gexf


try:
    import psyco
    psyco.full()
except ImportError:
    pass

__author__ = """Thomas Aynaud (thomas.aynaud@lip6.fr), David Combe (david.combe@gmail.com), Universite Jean Monnet, Saint-Etienne """



#nb de tests de deplacements de noeuds maximum
DC_MAX_APPAIRING=10000
DEBUG=0
COEFF_VARIANCE_INTER=1

# 1 pour utiliser que la formule
# 2 pour utiliser la super methode probabiliste
#UTILISER_CALINSKI=19 #100 #1 #101 #105 #19 #20 # 19 pour la formule du tableau de CL corrigee
UTILISER_CALINSKI=6
USE_MODULARITY_ONLY=0
USE_CALINSKI_ONLY=0
ECRIRE_LOG_COMPLET=1
CSV_SEPARATOR=","
UNE_EVAL_SUR=100000

EXPORTER_LEVELS_GEXF=1


if(UTILISER_CALINSKI==7):
    global coefProbaCal
    coefProbaCal = 1
if(UTILISER_CALINSKI==17):
    global coefProbaCal
    coefProbaCal = 0.35


PASS_MAX = 10
MIN = 0.0000001

logCount=0 

# dire si l'on veut se servir du processus de sauvegarde du top des partitions
STATESAVER=1
state_saver = None


currentGlobAuthorIndex=None






#global globDistanceMatrix
#global globAuthorIndex




def setDebug(debugLevel):
    global DEBUG
    DEBUG=debugLevel

def partition_at_level(dendogram, level) :
    """Return the partition of the nodes at the given level

    Level 0 is the first partition, and the best is len(dendogram) - 1

    :param dendogram: a list of partitions, ie dictionnaries where keys of the i+1 are the values of the i.
    :type dendogram: list of dictionary
    :param level: an integer which belongs to [0..len(dendogram)-1]
    :type level: integer
    :rtype: dictionary
    :return: a dictionary where keys are the nodes and the values are the set it belongs to
    """
    
    partition = dendogram[0].copy()
    for index in range(1, level + 1) :
        for node, community in partition.iteritems() :
            partition[node] = dendogram[index][community]
    return partition


def modularity(partition, graph) :
    """Compute the modularity of a partition of a graph

    :param partition: the partition of the nodes, i.e a dictionary where keys are their nodes and values the communities
    :type partition: dictionary
    :param graph: the networkx graph which is decomposed
    :type graph: networkx graph
    :rtype: float
    :return: The modularity
    """

    inc = dict([])
    deg = dict([])
    links = graph.size(weight = True)
    if links < 1 :
        return -2
    for node in graph :
        com = partition[node]
        deg[com] = deg.get(com, 0.) + graph.degree(node, weight = True)
        for neighbor, datas in graph[node].iteritems() :
            weight = datas.get("weight", 1)
            if partition[neighbor] == com :
                if neighbor == node :
                    inc[com] = inc.get(com, 0.) + float(weight)
                else :
                    inc[com] = inc.get(com, 0.) + float(weight) / 2.

    res = 0.
    for com in set(partition.values()) :
        print com, links, inc.get(com, 0.), deg.get(com, 0.)
        res += (inc.get(com, 0.) / links) - (deg.get(com, 0.) / (links))**2
    return res


# def distance(a, b):
   # global globDistanceMatrix
   # global globAuthorIndex
   # return globDistanceMatrix[globAuthorIndex[a]][globAuthorIndex[b]]


#tous les index des noeuds d'une communaute
def nodesOfCommunity(com, stat):
    global globAuthorIndex
    
    liste=[]
    #print "Com:"+str(com)
    #pprint(stat.node2com.items())
    for v,c in stat.node2com.items():
        if c == com:
            #print str(v) +" est dans "+str(c)
            liste.append(globAuthorIndex[v])
    return liste






#  max de la variance intra ou la minimisation de la variance interclasse.



def com2nodesIndex(communaute, stat):
        "renvoie les index des noeuds presents dans une communaute selon le Status actuel"
        return nodesOfCommunity(communaute, stat)


def scoreQuandNoeudAvecSaCommunaute(communauteOrigine, nomNoeud, communauteVoisine, stat):
        '''
        communauteOrigine et communauteVoisine doivent etre differentes
        '''
        #communauteOrigine=-1 
        
        global globAuthorIndex
        indexNoeud=globAuthorIndex[nomNoeud]
        
        #pprint(com2nodesIndex(communauteOrigine, stat))
        #print "identifiant com voisine: "+str(communauteVoisine)
        #print "Voisine:"
        #pprint(com2nodesIndex(communauteVoisine, stat))
        com1=[]
        com1.extend(com2nodesIndex(communauteOrigine, stat))
        com2=[]
        com2.extend(com2nodesIndex(communauteVoisine, stat))
        com1.append(indexNoeud)
        # la variance interclasse doit etre la plus grande possible : 0: memes attributs
        # max de la variance: infini ?
        #return (varianceInterClasse2Classes(com1, com2))
        return varianceInterClasse2ClassesSurTout(com1, com2, stat)


    
def scoreQuandNoeudAvecCommunauteVoisine(communauteOrigine, nomNoeud, communauteVoisine, stat):
        '''
        Calcul de la variance inter-classes uniquementsur deux classes (le reste s'annule quand on fait la difference entre la 
        variance de deux configurations)
        communauteOrigine et communauteVoisine doivent etre differentes
        '''
        #print "Debut Calcul Score scoreQuandNoeudAvecCommunauteVoisine\nOrig"
        #print "identifiant com orig: "+str(communauteOrigine)
        
        #Grosse modif qui se voit
        communauteOrigine=-1
        
        
        global globAuthorIndex
        indexNoeud=globAuthorIndex[nomNoeud]
        
        #pprint(com2nodesIndex(communauteOrigine, stat))
        #print "identifiant com voisine: "+str(communauteVoisine)
        #print "Voisine:"
        #pprint(com2nodesIndex(communauteVoisine, stat))
        com1=[]
        com1.extend(com2nodesIndex(communauteOrigine, stat))
        com2=[]
        com2.extend(com2nodesIndex(communauteVoisine, stat))
        com2.append(indexNoeud)
        # la variance interclasse doit etre la plus grande possible : 0: memes attributs
        # max de la variance: infini ?
        #return (varianceInterClasse2Classes(com1, com2))
        return varianceInterClasse2ClassesSurTout(com1, com2, stat)


def tousLesIndices(stat):
    """
    Returns the indices of all the elements of the graph
    """
    return stat.node2com.keys()
    #s=stat.node2com.values()
    global globAuthorIndex
    global globTfIdfTab
    
    #pprint(globAuthorIndex)
    #pprint(stat.node2com.values())
    #glob  node->index
    return [globAuthorIndex[x] for x in stat.node2com]
    #return stat.node2com.values()
    #def varianceGroupe():
    #def distanceListePointsCentre(indexsCommunaute, centre):


def varianceInterClasse2ClassesSurTout(com1, com2, stat):
        '''
        Différence de la moyenne des variances, les autres points / communautes sont omis
        '''
        #l=[]
        #for i in tousLesIndices(stat):
        #    l.append()
        global centreDeGraviteTot
        
        centreDeGraviteCom1 = moyenne_points(com1)
        centreDeGraviteCom2 = moyenne_points(com2)
        
        #distanceCosEntreCentroids(centreDeGraviteTot, centreDeGraviteCom1)
        
        nbTotalDePoints = len(stat.node2com.keys())
        
        distancePourCom1=distanceEucEntreCentroids(centreDeGraviteTot, centreDeGraviteCom1)
        distancePourCom2=distanceEucEntreCentroids(centreDeGraviteTot, centreDeGraviteCom2)
        
        var1 = (( len(com1)/ nbTotalDePoints)* distancePourCom1* distancePourCom1 ) 
        var2 = (( len(com2)/ nbTotalDePoints)* distancePourCom2* distancePourCom2 )
        
        return var1 + var2
        
        
        
        

def varianceInterClasse2Classes(com1, com2):
        '''
        Plus utilise
        '''
        #print "Debut variance intra 2 classes"
        #print "com1: "
        #pprint(com1)
        #print "com2: "
        #pprint(com1)
        # http://www.math-info.univ-paris5.fr/smel/cours/ts/node14.html
        # la variance des moyennes (variance inter-classes)
        cardCom1=len(com1)
        cardCom2=len(com2)
    
        #distance ?
        
        moyenneCom1=moyenne_points(com1)
        moyenneCom2=moyenne_points(com2)
        #print "long com1: "+str(len(com1))
        #for i in com1:
        #        print "moy1 "
        #        pprint(globTfIdfTab[i])
        #print "long com2: "+str(len(com2)) 
        #for i in com2: 
        #        print "moy2 "
        #        pprint(globTfIdfTab[i],depth=5000)       
        out=(cardCom2/(cardCom1+cardCom2)) * distanceEucEntreCentroids(moyenneCom1,moyenneCom2) * distanceEucEntreCentroids(moyenneCom2,moyenneCom1)
        #print "Variance: "+str(out)+" distanceCosEntreCentroids: "+str(distanceCosEntreCentroids(moyenneCom1,moyenneCom2)) +" premiere partie: "+str((cardCom2/(cardCom1+cardCom2)))
        return out



def recalculateTfIdfValues(statusList, currentStatus):
        """
        DEPRECATED
        # calcul des nouveaux attributs apres un level du dendogramme
        # pour chaque communaute
        # on cree une ligne qui est la somme de tous les elts
        """
        global status_list
        global globTfIdfTab
        global globAuthorIndex
        global currentGlobAuthorIndex
        nouveauxAttributs=[]
        compteNoeuds=[]
        
        if(DEBUG > 4):
            pprint(status_list)
            pprint(status_list[-1])
        dernierDico = status_list[-1]#.node2com
        
        if(DEBUG > 4):
            pprint(currentStatus.node2com)
        
        print str(len(globTfIdfTab[0]))+" attributs"
        
        for k in currentStatus.node2com.keys():
            nouveauxAttributs.append(zeros([len(globTfIdfTab[0])]))
            compteNoeuds.append(0)
            
        if(DEBUG > 4):
            pprint( nouveauxAttributs )
            pprint(compteNoeuds)
        
        #somme des attributs des differents noeud d'une nouvelle communaute
        for i in dernierDico: #boucle sur les noeuds
            for j in range(len(globTfIdfTab[0])): #boucle sur les attributs
                #!aCorriger
                nouveauxAttributs[dernierDico[i]][j]+=globTfIdfTab[globAuthorIndex[i]][j]
            compteNoeuds[dernierDico[i]] = compteNoeuds[dernierDico[i]]+1 
        
        # moyenne
        for i in range(len(nouveauxAttributs[0])):#boucle sur les attributs
            for j in range(len(nouveauxAttributs)):#boucle sur les noeuds
                 nouveauxAttributs[j][i]=nouveauxAttributs[j][i]/compteNoeuds[j]
            
        
        
        #globAuthorIndex, globTfIdfTab
        
        nbAut=len(currentGlobAuthorIndex)
        
        #reinit
        #globAuthorIndex=[]
        currentGlobAuthorIndex={}
        for i in range(nbAut):
            currentGlobAuthorIndex[i]=i
            
        
        
        
        
        
        if(DEBUG > 4):
            pprint(compteNoeuds)
            pprint(nouveauxAttributs)
        globTfIdfTab = nouveauxAttributs
        
        return nouveauxAttributs




def varianceInterClasseToutesClasses(partition):
        """
        Calcule la variance interclasse pour une configuration de status donnee.
        """
        print "Calcul de la variance interclasse en mode lent..."
        #print "dest point courant: "+str(destinationPointCourant)
        #print "status.node2com.keys()"
        #pprint(status.node2com.items())
        nbTotalDePoints=len(status.node2com.keys())
        communautes=set([])
        for k,v in status.node2com.items():
                #print(str(k)+" <->  "+str(v))
                #if v==(-1):
                #        communautes.add(destinationPointCourant)
                #pointCourant=k
                        
                #else:
                        communautes.add(v)
        
        varCommunaute = 0
        global centreDeGraviteTot
        
        #pour chaque communaute
        for i in communautes:
                #print("Communaute "+str(i))
                #pprint(com2nodesIndex(i, status))
                #assert(len(com2nodesIndex(i, status))>0)
                indexDesPointsDeLaCom = com2nodesIndex(i, status)
                #if(i == destinationPointCourant):
                #        indicePointCourant = globAuthorIndex[pointCourant]
                #        indexDesPointsDeLaCom.append(indicePointCourant)
                centreDeGraviteCom = moyenne_points(indexDesPointsDeLaCom)
                distancePourCom=distanceEucEntreCentroids(centreDeGraviteTot, centreDeGraviteCom)
                com=com2nodesIndex(i, status)
                #if destinationPointCourant==i:
                #        com.append(globAuthorIndex[pointCourant])
                #print "Etape calcul var: "+str(( len(com)/ nbTotalDePoints)* distancePourCom* distancePourCom)
                varCommunaute = varCommunaute + (( len(com)/ nbTotalDePoints)* distancePourCom* distancePourCom)
        print "OK: "+str(varCommunaute)
        return varCommunaute




def varianceInterClasseToutesClasses(status):#, destinationPointCourant):
        """
        Calcule la variance interclasse pour une configuration de status donnee.
        """
        print "Calcul de la variance interclasse en mode lent..."
        #print "dest point courant: "+str(destinationPointCourant)
        #print "status.node2com.keys()"
        #pprint(status.node2com.items())
        nbTotalDePoints=len(status.node2com.keys())
        communautes=set([])
        for k,v in status.node2com.items():
                #print(str(k)+" <->  "+str(v))
                #if v==(-1):
                #        communautes.add(destinationPointCourant)
                #pointCourant=k
                        
                #else:
                        communautes.add(v)
        
        varCommunaute = 0
        global centreDeGraviteTot
        
        #pour chaque communaute
        for i in communautes:
                #print("Communaute "+str(i))
                #pprint(com2nodesIndex(i, status))
                #assert(len(com2nodesIndex(i, status))>0)
                indexDesPointsDeLaCom = com2nodesIndex(i, status)
                #if(i == destinationPointCourant):
                #        indicePointCourant = globAuthorIndex[pointCourant]
                #        indexDesPointsDeLaCom.append(indicePointCourant)
                centreDeGraviteCom = moyenne_points(indexDesPointsDeLaCom)
                distancePourCom=distanceEucEntreCentroids(centreDeGraviteTot, centreDeGraviteCom)
                com=com2nodesIndex(i, status)
                #if destinationPointCourant==i:
                #        com.append(globAuthorIndex[pointCourant])
                #print "Etape calcul var: "+str(( len(com)/ nbTotalDePoints)* distancePourCom* distancePourCom)
                varCommunaute = varCommunaute + (( len(com)/ nbTotalDePoints)* distancePourCom* distancePourCom)
        print "OK: "+str(varCommunaute)
        return varCommunaute
        
        
        

def setGlobTfIdfTab(myGlobTfIdfTab):
    '''
        Utilise par le fichier de tests unitaires
    '''
    global globTfIdfTab 
    globTfIdfTab=myGlobTfIdfTab


def moyenne_points(indices):
    '''
    Calcule le centre de gravite d'un ensemble de points,
    a savoir la moyenne des differents attributs
    Les vecteurs sont recuperes dans la variable globale globTfIdfTab
    '''
    global globTfIdfTab 
    nb_att=len(globTfIdfTab[indices[0]])
    tab=zeros([nb_att])
    
    for i in indices:
        for j in range(nb_att):
            tab[j]=tab[j]+globTfIdfTab[i][j]
    for j in range(nb_att):
        tab[j]=tab[j]/len(indices)
    return tab



def varianceTotale(status):
        """
        Variance de tous les points V_T
        """
        # Si la variance a deja ete calculee on renvoie sa valeur
        global laVarianceTotale
        
        if laVarianceTotale!= None:
            return laVarianceTotale
        # try:
            # laVarianceTotale
            # return laVarianceTotale
        # except NameError:
                
        global centreDeGraviteTot
        moyenne = centreDeGraviteTot
        somme = 0
        for i in globTfIdfTab:
                distancePourCom=distanceEucEntreCentroids(moyenne, i)
                somme = somme + distancePourCom * distancePourCom
        
        varianceTot = somme / len(globTfIdfTab)
        
        #pprint(varianceTot)
        #exit()
        
        laVarianceTotale = varianceTot
        
        
        if(2==0): # ???
        
                #moyenne des carrés moins le carré des moyennes
                varCommunaute = varCommunaute + (( len(com)/ nbTotalDePoints)* distancePourCom* distancePourCom)
                somme=distanceEucEntreCentroids(centreDeGraviteTot, centreDeGraviteCom)
                
                #on fait un tableau de tableau des noeuds selon la communaute d'appartenance
                global globAuthorIndex
                partition=[]
                nbTotalDePoints=0
                imax=-1
                #for v,c in stat.node2com.items():
                #pprint(status.node2com) 
                for (k,v) in status.node2com.items():
                      if v>imax:
                                imax=v
                      nbTotalDePoints=nbTotalDePoints+1
                for i in range(imax+1):
                        partition.append([])
                
                for k,v in status.node2com.iteritems():
                        partition[v].append(globAuthorIndex[k])
                
                #Calcul de la variance
                variance=0
                
                
                global centreDeGraviteTot
                centreDeGraviteCom1 = moyenne_points(partition[i])
                
                for i in range(imax+1):
                        distancePourCom1=distanceEucEntreCentroids(centreDeGraviteTot, centreDeGraviteCom1) 
                        variance=variance+(( len(partition[i])/ nbTotalDePoints)* distancePourCom1* distancePourCom1 ) 
                
                #pprint(variance)
                #distanceCosEntreCentroids(centreDeGraviteTot, centreDeGraviteCom1)
                
                #nbTotalDePoints = len(stat.node2com.keys())
                laVarianceTotale=variance
                return variance
        
        return varianceTot


    
def distanceEucEntreCentroids(centro1,centro2):
        '''
        Calcule avec scipy la distance du cosinus entre deux vecteurs.
        '''
        return euclidean(centro1,centro2)    
    
    
    
    
#def varianceInterClasse(com1, com2):    
        # http://www.math-info.univ-paris5.fr/smel/cours/ts/node14.html
        # la variance des moyennes (variance inter-classes)
    
        
    
#    somme=0
#    for com in communautesVoisinesDuPoint:
#            somme= somme + ()()()
#    V_
    
#    pass    
    
    
    
def ecrireLog(currentLevel, com, node, associatedWith, varianceInter, mod, critGlob, varTot, estMeilleur,composanteTexte,valeurCalinski, ari, homogeneity, completeness, vMeasure, nbDeClasses):
        global log
        global logCount
        global CSV_SEPARATOR
        log.write( str(logCount) + CSV_SEPARATOR + str(currentLevel) +CSV_SEPARATOR+ str(com) +CSV_SEPARATOR+ str(node) +CSV_SEPARATOR+ str(associatedWith) +CSV_SEPARATOR+ str(varianceInter) +CSV_SEPARATOR+ str(mod) +CSV_SEPARATOR+ str(critGlob) + CSV_SEPARATOR +  str(varTot) + CSV_SEPARATOR +str(estMeilleur)+   CSV_SEPARATOR+str(composanteTexte)+CSV_SEPARATOR+str(varianceInter/varTot)+ CSV_SEPARATOR+str(valeurCalinski)+ CSV_SEPARATOR+str(ari)+ CSV_SEPARATOR+str(homogeneity)+ CSV_SEPARATOR+str(completeness)+ CSV_SEPARATOR+str(vMeasure)+ CSV_SEPARATOR+str(nbDeClasses)+    "\n" )
        logCount=logCount+1
        #log.flush()
        
        

def getListeIndexsPointsCommunaute(com):
        return nodesOfCommunity(com, stat)

def getListeIndexsPointsCommunautePlusPoint(com, pointAdditionnel):
        
        return nodesOfCommunity(com, stat)

#def distanceCommunautesTexte(communauteOriginaleDeA, a, b, stat):
#    print "Distance variance entre "
#    nodesOfCommunity(communauteOriginaleDeA, stat)
#    print "et\n"
#    nodesOfCommunity(a, stat)
#    global globDistanceMatrix
#    global globAuthorIndex
#    pprint(a)
#    print "Moy"
#    
#    pprint(a)
#    indexDeA=globAuthorIndex[a]
#    listePoints=[]
#    listePoints.append(indexDeA)
#    
#    
#    
#    
#    pprint(varianceIntraClasse())
#    
#    pprint(moyenne_points(listeNomsToListeIndexNoeuds(nodesOfCommunity(communauteOriginaleDeA, stat))))
#    pprint(b)
#    #pprint(globDistanceMatrix)
#    #pprint(globAuthorIndex[a])
#    #pprint(globAuthorIndex[b])
#    return globDistanceMatrix[globAuthorIndex[a]][globAuthorIndex[b]]

def listeNomsToListeIndexNoeuds(names):
        global globAuthorIndex
        indexs=[]
        for i in names:
                indexs.append(globAuthorIndex[i])
        return indexs


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
    #G.add_edge(3,4)
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
    
        
def graphExample2(): 
    global globTfIdfTab
    
    globTfIdfTab=[
                            [2,1],
                            [3,1],
                            [22,21],
                            [20,22]#,
                            #[21,22] 
                    ]
        
    global globAuthorIndex
    globAuthorIndex={0: 0,
                     1: 1,
                     2: 2,
                     3: 3#,#4: 4
                     }  


    G=nx.Graph()
    
    G.add_node(0) 
    G.add_node(1)
    G.add_node(2)
    G.add_node(3)
    #G.add_node(4)

    G.add_edge(0,1)
    G.add_edge(1,2)
    G.add_edge(2,3)
    #G.add_edge(3,4)
    G.add_edge(0,2)
    return G
    
def calinskiGeneratorProbability(classesNumber, varInterClasses, varTot, n):
        """
        Donne une probabilite pour une partition donnee
        Entrees:
        - v1 et v2, les degres de liberte, cad nb de classes et nb d'elements remixes
        - var inter
        - var totale
        - x, la valeur observee, cad le resultat de la formule de calinski
        http://labh-curien.univ-st-etienne.fr/wiki-SocialMining/index.php/09/02/2012_R%C3%A9union_MG_CL
        """
        looseAssertlt(varInterClasses,varTot)
        
        assert(classesNumber<=n)
        
        #print "classesNumber, varInterClasses, varTot, n"
        #print classesNumber, varInterClasses, varTot, n
        
        #Si on est dans les extremes: -3.14
        if (n == classesNumber):
            print "Premier coup: extremum --> -3.14"
            return -31401
        if (classesNumber == 1):
            print "extremum: too far --> -3.14"
            return -31402
        
        #print "classesNumber, varInterClasses, varTot, n"
        #print classesNumber, varInterClasses, varTot, n
        
        #petit ajustement pour eviter div par zero
        calinski = ((   ((varInterClasses-0.000001/(classesNumber-1))/((varTot-varInterClasses+0.000001)/(n-classesNumber)))     ))
        #print calinski
        
        # logcdf version scipy
        res = f.logcdf(   calinski    ,   classesNumber-1, n-classesNumber  , scale=10000000)
        
        # logcdf version R
        # res=r("pf(calinski    ,   classesNumber-1, n-classesNumber  ,log.p=TRUE)")
        #rcall = "pf("+str(calinski)   +" , "+  str(classesNumber-1)+", "+ str(n-classesNumber) +"  ,log.p=TRUE)"
        #res=robjects.r(rcall)[0] #"pf(calinski    ,   classesNumber-1, n-classesNumber  ,log.p=TRUE)")

        
        
        #print "Resultat probabilite Calinski: "+str(res)
        
        if (res==0):#Problem
            print "calinski a zero"
            print "("+str(varInterClasses)+" / "+str(classesNumber)+"-1)/(("+str(varTot)+"-"+str(varInterClasses)+")/("+str(n)+"-"+str(classesNumber)+"))"
            exit()
        if (res != res): #infini
            print "calinski a NaN"
            return -31401
            exit()
            
            
        #print (1+math.log(res,10)/10)
        #raw_input("ert")
        return res #1+math.log(res,10)/10
        
        
        #return res#1-res
        
        
        
        
#def calculDuScoreAdditionVarianceModularitePondere(status, candidat, communauteCandidat, varInterClasses, mod=0):
def calculDuScoreAdditionVarianceModularitePondere(classesNumber, varInterClasses, mod, status, modAttr):
        '''
        Choix du mix des scores de modularite et de variance
        Renvoie automatiquement -1 si l'on a affaire a la partition discrete ou a la partition jointe 
        dans le cas ou l'on utilise Calinski.
        '''
        
        
        return ((mod)+(modAttr)),modAttr
        #return ((modAttr+1)),modAttr
        
        
        #print "------->"
        #print varInterClasses
        assert(classesNumber>0)
        
        assert(mod>-1.001)
        
        assert(mod<1.001)
        
        #print str(classesNumber)+"Nb de classes"
        global laVarianceTotale
        global varStatus
        #print "in calculDuScoreAdditionVarianceModularitePondere"
        #calcul de la variance
        #varInterClassesSlow=varianceInterClasseToutesClasses(status)#, communauteCandidat)
        
        
        #print "Verif ecart var: slow "+str(varInterClassesSlow)+" <> inc "+str(varInterClasses)

        
        #assert( abs(varInterClasses-varInterClassesSlow)<0.001   )
        
        #return varInterClasses
        
        
        #calcul de la modularite
        #mod=__modularity(status)
        
        #calcul de la variance totale
        varTot=laVarianceTotale
        
        #pprint(varInterClasses)
        #pprint(varTot)
        looseAssertlt(varInterClasses,varTot)
        
        
        # la variance est à diviser par la variance totale (intra-classe + inter-classe)
        #print str(mod)+" <-> "+str(varInterClasses)+" <-> "+str(varTot)
        
        #nb d'elements
        global globTfIdfTab
        n=len(globTfIdfTab)
        
        #nb de classes
        c=classesNumber
        #c=numberOfClasses(status)-1
        if c<2:
            c==2
        
        #return mod+COEFF_VARIANCE_INTER*varInterClasses/varTot
        #varInterClasses
        #varTot
        
        #print "varInterClasses: "+str(varInterClasses)
        #print "mod: "+str(mod)
        #print "vic"
        #print varInterClasses
        
        composanteTexte=-3.14
        
        if (USE_MODULARITY_ONLY==1):
            try:
                composanteTexte=calinskiGeneratorProbability(classesNumber, varInterClasses, varTot, n)
            except:
                pass
            return mod,composanteTexte
        if (USE_CALINSKI_ONLY==1):
            composanteTexte=calinskiGeneratorProbability(classesNumber, varInterClasses, varTot, n)
            return composanteTexte,composanteTexte
        
        # =============================================================================== #
        
        
        if (UTILISER_CALINSKI==1):
            if (n == c):
                return -1,composanteTexte
            if (c == (1)):
                return -1,composanteTexte
            res =  mod+ ((   COEFF_VARIANCE_INTER*((varInterClasses/(c-1))/((0.0000001+varTot-varInterClasses)/(n-c)))     ))
        
        elif (UTILISER_CALINSKI==100):
            if (n == c):
                return -1,composanteTexte
            if (c == (1)):
                return -1,composanteTexte
            res =   ((   COEFF_VARIANCE_INTER*((varInterClasses/(c-1))/((varTot-varInterClasses)/(n-c)))     ))
        
            #print str(mod)+"+ ((   COEFF_VARIANCE_INTER*(("+str(varInterClasses)+"/("+str(c)+"-1))/(("+str(varTot)+"-"+str(varInterClasses)+")/("+str(n)+"-"+str(c)+")))     )))"
            
        elif(UTILISER_CALINSKI==17):
                global coefProbaCal
                # Coef sur mod et proba calinski
                if (n == c):
                    return (-1,composanteTexte)
                #print "nbClasses, vic, vt, n"
                #print classesNumber, varInterClasses, varTot, n
                scoreProba=calinskiGeneratorProbability(classesNumber, varInterClasses, varTot, n)
                assert(scoreProba == scoreProba)
                #print "scoreProba"
                #print scoreProba
                
                # manage negative values
                #if (mod<0):
                #    if(scoreProba<0):
                #        scoreProba=-scoreProba
                
                #composanteTexte=scoreProba
                composanteTexte=coefProbaCal*(scoreProba/400 +1)
                res=coefProbaCal*(scoreProba/400 +1) + (1-coefProbaCal)*mod
                #res= scoreProba * mod
            
        elif(UTILISER_CALINSKI==7):
                global coefProbaCal
                # Coef sur mod et proba calinski
                if (n == c):
                    return (-1,composanteTexte)
                #print "nbClasses, vic, vt, n"
                #print classesNumber, varInterClasses, varTot, n
                scoreProba=calinskiGeneratorProbability(classesNumber, varInterClasses, varTot, n)
                assert(scoreProba == scoreProba)
                #print "scoreProba"
                #print scoreProba
                
                # manage negative values
                #if (mod<0):
                #    if(scoreProba<0):
                #        scoreProba=-scoreProba
                
                composanteTexte=scoreProba
                
                res=coefProbaCal*scoreProba + (1-coefProbaCal)*mod
                #res= scoreProba * mod
            
        elif(UTILISER_CALINSKI==6):
                if (n == c):
                    return -1,composanteTexte
                #print "nbClasses, vic, vt, n"
                #print classesNumber, varInterClasses, varTot, n
                scoreProba=calinskiGeneratorProbability(classesNumber, varInterClasses, varTot, n)
                assert(scoreProba == scoreProba)
                #print "scoreProba"
                #print scoreProba
                
                # manage negative values
                #if (mod<0):
                #    if(scoreProba<0):
                #        scoreProba=-scoreProba
                
                res=scoreProba 
                #res= scoreProba * mod
        elif(UTILISER_CALINSKI==2):
                if (n == c):
                    return -1,composanteTexte
                #print "nbClasses, vic, vt, n"
                #print classesNumber, varInterClasses, varTot, n
                scoreProba=calinskiGeneratorProbability(classesNumber, varInterClasses, varTot, n)
                assert(scoreProba == scoreProba)
                #print "scoreProba"
                #print scoreProba
                
                # manage negative values
                if (mod<0):
                    if(scoreProba<0):
                        scoreProba=-scoreProba
                    
                res=scoreProba + mod
                #res= scoreProba * mod
        elif(UTILISER_CALINSKI==3):
                    res = mod + varInterClasses/varTot
                    
        elif(UTILISER_CALINSKI==16):
            moyVar=varStatus.meanOfVariancesOfGroups()
            #print varTot,mod,classesNumber,varInterClasses,moyVar
            res= ((varTot*mod)/(classesNumber*varInterClasses*(moyVar+.00001) ))
        elif(UTILISER_CALINSKI==18):
            #moyVar=varStatus.meanOfVariancesOfGroups()
            #print varTot,mod,classesNumber,varInterClasses,moyVar
            res= ((varTot*mod)/(classesNumber*varInterClasses ))
        elif(UTILISER_CALINSKI==19):
            #moyVar=varStatus.meanOfVariancesOfGroups()
            #print varTot,mod,classesNumber,varInterClasses,moyVar
            
            res= ((varInterClasses*mod)/(classesNumber*varTot ))
        elif(UTILISER_CALINSKI==105):
            #moyVar=varStatus.meanOfVariancesOfGroups()
            #print varTot,mod,classesNumber,varInterClasses,moyVar
            #res= ((varInterClasses*mod)/(varTot )) + n - math.log((1/classesNumber+classesNumber)) #* log(classesNumber)
            res= ((varInterClasses*mod)/(varTot )) + n - abs(  math.log(n,10) - classesNumber )    # ((1/classesNumber+classesNumber)) #* log(classesNumber)      
        elif(UTILISER_CALINSKI==104):
            #moyVar=varStatus.meanOfVariancesOfGroups()
            #print varTot,mod,classesNumber,varInterClasses,moyVar
            
            res= ((varInterClasses*mod)/(classesNumber*varTot )) + n - abs(4 - classesNumber)
        elif(UTILISER_CALINSKI==101):
            #moyVar=varStatus.meanOfVariancesOfGroups()
            #print varTot,mod,classesNumber,varInterClasses,moyVar
            
            res= ((varInterClasses*mod)/(classesNumber*varTot )) + n - abs(3 - classesNumber)
        elif(UTILISER_CALINSKI==27):
            #pour mieux pénaliser un grand nombre de classes que le critère de base (du tableau)
            #moyVar=varStatus.meanOfVariancesOfGroups()
            #print varTot,mod,classesNumber,varInterClasses,moyVar
            
            res= ((varInterClasses*mod)/(varTot )) - 1.0/classesNumber
        elif(UTILISER_CALINSKI==20):    
            assert(varTot != 0)
            res= ((varInterClasses*mod)/(varTot ))
        elif(UTILISER_CALINSKI==21):    
            assert(varTot != 0)
            res= (varInterClasses*(1+mod)**2)
        elif(UTILISER_CALINSKI==5):
                        #res =  mod+((        COEFF_VARIANCE_INTER*varInterClasses/varTot        )*2-1)
                        res =  (mod+   COEFF_VARIANCE_INTER*((  varInterClasses/varTot   )*2-2)   ) / (1+COEFF_VARIANCE_INTER)
        else:
            print "UTILISER_CALINSKI value unknown"
            exit()
        #print res
        #assert(res>-10000)
        #assert(res<10000)
        assert(not math.isnan(res))
        return (res, composanteTexte)
        
        
        
        
        
        
def best_partition(graph, threshold,  authorIndex, tfIdfTab, partition = None, groundTruthArray=None) :
    '''Compute the partition of the graph nodes which maximises the modularity
    (or try..) using the Louvain heuristices

    :param graph: the networkx graph which is decomposed
    :type graph: networkx graph
    :param partition: the algorithm will start using this partition of the nodes. It's a dictionary where keys are their nodes and values the communities
    :type partition: dictionary, optional
    :rtype: dictionary
    :return: The partition, with communities numbered from 0 to number of communities
    
    
    DC
    Takes moreover
    - index of authors
    - tfIdfTab (attributs)
    - distanceMatrix (plus utilise)
    - threshold (plus utilise)
    '''
    #pprint(authorIndex)
    #pprint(groundTruthArray)
    #for i in authorIndex:
    #    print i, ",", groundTruthArray[authorIndex[i]][0:1]
    #pprint(authorIndex)
    #exit()
    
    # log file initialisation
    global log
    global CSV_SEPARATOR
    global state_saver
    global STATESAVER
    
    global centreDeGraviteTot
    #try:
    centreDeGraviteTot=None
    #except:
    #    pass
        
    global laVarianceTotale
    
    laVarianceTotale=None
    #del laVarianceTotale
    
    nameOfLog = "logMethod"+str(UTILISER_CALINSKI)+".csv"
    if USE_MODULARITY_ONLY==1:
        nameOfLog="drivenByMod_"+nameOfLog
    if USE_CALINSKI_ONLY==1:
        nameOfLog="drivenByCalinsinski_"+nameOfLog
    if ECRIRE_LOG_COMPLET==0:
        nameOfLog="onlyBetterPartitions_"+nameOfLog
    log = open(nameOfLog, "w")
    log.write("#i"+CSV_SEPARATOR+"pass"+CSV_SEPARATOR+"com"+CSV_SEPARATOR+"node"+CSV_SEPARATOR+"associatedWith"+CSV_SEPARATOR+"varianceInter"+CSV_SEPARATOR+"mod"+CSV_SEPARATOR+"critGlob"+CSV_SEPARATOR+"totalVariance"+CSV_SEPARATOR+"isBetter"+CSV_SEPARATOR+"composanteTexte"+CSV_SEPARATOR+"varInterSurVarTot"+CSV_SEPARATOR+ "calinski"+CSV_SEPARATOR+"ari"+CSV_SEPARATOR+"homogeneity"+CSV_SEPARATOR+"completeness"+CSV_SEPARATOR+"vMeasure"+CSV_SEPARATOR+"numberOfClasses\n")
    
    
    
    
    
    # on ne veut pas le petit exemple
    if(tfIdfTab!=None):
            print "Index des auteurs"
            global globAuthorIndex
            globAuthorIndex=authorIndex
            #pprint(globAuthorIndex)
            
            print "tfIdfTab"
            global globTfIdfTab
            globTfIdfTab = tfIdfTab
    
    else:
            #Example simplifie
            graph=graphExample()
    
    global globTfIdfTab
    assert(len(globTfIdfTab) == len(graph.nodes()))
    print "Lancement Louvain etendu\n"
    dendo = generate_dendogram(graph, threshold, part_init=partition, groundTruthArray=groundTruthArray)
    
    print "Fin Louvain etendu"
    
    
    global EXPORTER_LEVELS_GEXF
    
    graph2=graph.copy()
    try:
        for node in graph2.nodes():
            del graph2.node[node]['att']
    except:
        pass
        
    i=0
    partitions=[]
    if EXPORTER_LEVELS_GEXF:
        while(i<len(dendo)):
            partitions.append( partition_at_level(dendo, i))
            #print "mod",str(modularity(partition_at_level(dendo, i), graph2))

            print "Append info to graph2"
        
            for node in graph2.nodes():
                #print type(node), type(i), type(partitions[i][node])
                graph2.node[node]['louvain'+str(i)]=str(partitions[i][node])
                
            #pprint(graph2.node[node])


            i=i+1
        #pprint(partitions)
        print "Exporting to "+"exportLevelsTotem.gexf"+"\n"
        write_gexf(graph2,"exportLevelsTotem.gexf")
    
    
    
    if(STATESAVER==1):
        
        state_saver.saveOutput()
    
    return partition_at_level(dendo, len(dendo) - 1 )








def load_binary(data) :
    """Load binary graph as used by the cpp implementation of this algorithm

    :param data: the file containing the data
    :type data: string or file
    :rtype: networkx.Graph
    :return: The graph

    """

    import types
    if type(data) == types.StringType :
        data = open(data, "rb")
    import array
    reader = array.array("I")
    reader.fromfile(data, 1)
    num_nodes = reader.pop()
    reader = array.array("I")
    reader.fromfile(data, num_nodes)
    cum_deg = reader.tolist()
    num_links = reader.pop()
    reader = array.array("I")
    reader.fromfile(data, num_links)
    links = reader.tolist()
    graph = nx.Graph()
    graph.add_nodes_from(range(num_nodes))
    prec_deg = 0
    for index in range(num_nodes) :
        last_deg = cum_deg[index]
        neighbors = links[prec_deg:last_deg]
        graph.add_edges_from([(index, int(neigh)) for neigh in neighbors])
        prec_deg = last_deg
    return graph



def renumber(dictionary) :
    """Renumber the values of the dictionary from 0 to n

    :param dictionary: the partition
    :type dictionary: dictionary
    :rtype: dictionary
    :return: The modified partition and the number of classes

    """

    count = 0
    ret = dictionary.copy()
    new_values = dict([])
    for key in dictionary.keys() :
        value = dictionary[key]
        new_value = new_values.get(value, -1)
        if new_value == -1 :
            new_values[value] = count
            new_value = count
            count = count + 1
        ret[key] = new_value
    return ret, count






def generate_dendogram(graph, threshold, part_init = None, groundTruthArray=None) :
    """Find communities in the graph and return the associated dendogram

    :param graph: the networkx graph which will be decomposed
    :type graph: networkx graph
    :param part_init: the algorithm will start using this partition of the nodes. It's a dictionary where keys are their nodes and values the communities
    :type part_init: dictionary, optional
    :rtype: list of dictionaries
    :return: a list of partitions, ie dictionnaries where keys of the i+1 are the values of the i. and where keys of the first are the nodes of graph
    
    """
    global globAuthorIndex
    global currentGlobAuthorIndex 
    global globTfIdfTab
    global varStatus
    
    global state_saver, STATESAVER
    if STATESAVER==1:
        state_saver = StateSaver(numberOfPartitionsToSave=10, groundTruth=groundTruthArray)
    
    
    currentGlobAuthorIndex=globAuthorIndex
    
    current_graph = graph.copy()
    status = Status()
    
    
    from relabel import relabel_nodes
    
    #pprint(current_graph.nodes())
    
    
    
    current_graph=relabel_nodes(current_graph, globAuthorIndex)
    
    
    
    #H2=current_graph.copy()
    #print "Exporting to "+"exportScriptGetoor.gexf"+"\n"
    #write_gexf(H2,"exportScriptGetoor.gexf")
    
    #pprint(H2.nodes())
    
    
    #si on a fourni une partition par laquelle commencer:
    if part_init == None :
        partInit={}
        for i in range(len(globAuthorIndex)):
            partInit[i]=i
        status.init(current_graph, partInit)
        
        if(DEBUG > 1):
            print "Initialisation de la structure pour la variance..."
        varStatus=MetaVarIncrement(current_graph, globTfIdfTab, status.node2com, globAuthorIndex, metaNodesHistory(status.node2com))
    else: #sinon
        partInit=part_init
        status.init(current_graph, partInit)
        
        if(DEBUG > 1):
            print "Initialisation de la structure pour la variance..."
        varStatus=MetaVarIncrement(current_graph, globTfIdfTab, part_init, globAuthorIndex, metaNodesHistory(status.node2com))
    
    mod = __modularity(status)
    
    
    
    
    
    global status_list
    status_list = list()
    
    status_list.append(globAuthorIndex)
    
    
    
    
    # variable globale pour savoir dans quel level on est
    global currentLevel
    currentLevel = len(status_list)
    
    
    
    #assert(len(status.node2com) == len(varStatus.comGC))
    #pprint(status.node2com)
    #pprint(varStatus.comGC)
    
    
    
    
    #print str(varStatus)
    #print len(graph.nodes())
    #pprint(status.node2com)
    
    
    
    new_curCritGlob = one_level(current_graph, status, threshold, status_list, groundTruthArray=groundTruthArray)
    #new_mod = __modularity(status)
    
    #if False: #pour ne faire qu'une passe
    
    # si il y a des classes qui ont disparu, on decale les etiqutettes de classe
    partition,_ = renumber(status.node2com)
    #print "renumberedPartition modularity"
    #pprint(partition)
    status_list.append(partition)
    
    #! enlever ca
    #return status_list[:]
    
    # preparer la situation pour un tour de phase 2
    # on change le currentGlobAuthorIndex
    #currentGlobAuthorIndex
    
    # on change le graphe -> inducedGraph
    
    
    curCritGlob=new_curCritGlob 
    
    #mod = new_mod
    current_graph = induced_graph(partition, current_graph)
    
    # induced varstatus
    
    
    #print(varStatus.currentPartition)
    varStatus.inducedVarStatus(partition)
    #pprint(partition)
    #print(len(partition))
    #print(len(varStatus.comGC))
    #print(len(status.node2com))
    #assert(len(status.node2com) == len(varStatus.comGC))
    
    
    currentGlobAuthorIndex=[]
    for i in range(len(current_graph.nodes())):
        currentGlobAuthorIndex.append(i)
    
    
    status.init(current_graph)
    #print "-----"
    #pprint(status.node2com)
    #print "varStatus.currentPartition.node2com"
    #print(varStatus.currentPartition.node2com)
    #print "hierarchy"
    #print(varStatus.hierarchy)
    #print "status_list"
    #print(status_list)
    #raw_input("induced")
    #recalcul des vecteurs Tf-Idf
    #recalculateTfIdfValues(status_list, status)
    
    while True :
        currentLevel = len(status_list)
        
        new_curCritGlob = one_level(current_graph, status, threshold,  status_list, groundTruthArray=groundTruthArray)
        new_mod = __modularity(status)
        #print new_curCritGlob, new_mod
        
        if (   math.isnan(new_curCritGlob)):
            print "nan"
            exit()
        if new_curCritGlob - curCritGlob <MIN:
            break
        #if new_mod - mod < MIN :
        #    break
        partition,_ = renumber(status.node2com)
        status_list.append(partition)
        curCritGlob = new_curCritGlob
        #mod = new_mod
        current_graph = induced_graph(partition, current_graph)

        varStatus.inducedVarStatus(partition)
        #assert(len(status.node2com) == len(varStatus.comGC))
        
        #recalculateTfIdfValues(status_list, status)
        status.init(current_graph)
        
        
        
    return status_list[:]


def induced_graph(partition, graph) :
    """Produce the graph where nodes are the communities

    there is a link of weight w between communities if the sum of the weights
    of the links between their elements is w

    :param partition: a dictionary where keys are graph nodes and  values the part the node belongs to
    :type partition: dictionary
    :param graph: the initial graph
    :type graph: networkx graph
    :rtype: networkx.Graph
    :return: a networkx graph where nodes are the parts
    """
    #pprint(partition)
    #exit()
    
    
    new_graph = nx.Graph()
    new_graph.add_nodes_from(partition.values())
    for node1, node2, datas in graph.edges_iter(data = True) :
        weight = datas.get("weight", 1)
        com1 = partition[node1]
        com2 = partition[node2]
        weight_prec = new_graph.get_edge_data(com1, com2, {"weight":0}).get("weight", 1)
        new_graph.add_edge(com1, com2, weight=weight_prec + weight)
        #print (weight_prec + weight)
    
    
    return new_graph







def one_level(graph, status, threshold, status_list, groundTruthArray=None) :
    """Compute one level of communities

    :param graph: the graph we are working on
    :type graph: dictionary
    :param status: a named tuple with node2com, total_weight, internals, degrees set
    :type status: Status
    :return: nothing, the status is modified during the function
    """
    
    
    #pprint(status_list)
    #raw_input("induced")
    
    global logCount #used for log file
    global currentLevel
    global globTfIdfTab
    global currentGlobAuthorIndex
    global varStatus
    
    global UNE_EVAL_SUR
    
    global STATESAVER
    global state_saver
    
    ari=-2
    homogeneity=-2
    completeness=-2
    vMeasure=-2
    #pprint(globTfIdfTab)
    
    appairing_count=0
    
    global centreDeGraviteTot
    #try:
    #        centreDeGraviteTot
    #except NameError:
    if centreDeGraviteTot == None:      
            centreDeGraviteTot = moyenne_points(tousLesIndices(status))
    
    
    #! Ne calculer ca qu'une fois
    global laVarianceTotale
    laVarianceTotale=varianceTotale(status)
    if(DEBUG > 0):
        print "Var totale: "+str(laVarianceTotale)
    
    
    modif = True
    nb_pass_done = 0
    print str(status)
    cur_mod = __modularity(status)
    print cur_mod
    
    
    new_mod = cur_mod
    
    
    # creation d'une partition discrete
    partitionDisjointe = {}
    j=0
    for i in graph.nodes():
        partitionDisjointe[i] = j
        j=j+1
    
    ## initialisation de la structure pour la variance, avec la partition disjointe
    
    print "---------------------------------------------------"
    print "|        Etat initial                             |"
    print "| mod :                   "+str(cur_mod)
    print "---------------------------------------------------"
    #partition=APartition(partitionDisjointe)
    #varStatus=VarianceStatus(graph, globTfIdfTab, partition, globAuthorIndex)
    
    #print "varStatus:"
    #print str(varStatus)
    #assert(abs(varStatus.currentVarInter-laVarianceTotale)<0.1)
    
    cur_critereGlobal=-2
    cur_varInter = laVarianceTotale
    
    
    
    
    #exit()
    #print(varStatus)
    
    cur_varInter2=varStatus.currentVarInter
    
    
    
    #pprint(status.node2com)
    #print(str(varStatus))
    #assert(len(status.node2com) == len(varStatus.comGC))
    
    #exit()
    
    while modif and nb_pass_done != PASS_MAX and (appairing_count<DC_MAX_APPAIRING):
        new_crit = cur_critereGlobal
        
        #! enlever les refs a la mod la ou elles ne sont pas necessaires
        cur_mod = new_mod
        modif = False
        nb_pass_done = nb_pass_done + 1
        
        #if(DEBUG > 4):
        #    pprint(status_list)
        
        
        
        
        
        #varInterClassesSlow=varianceInterClasseToutesClasses(status)
        
        
        #print "Verif slow "+str(varInterClassesSlow)+" <> inc: "+str(cur_varInter)
        #assert(abs(varInterClassesSlow - cur_varInter )<.0001)
        
        looseAssertlt(cur_varInter,laVarianceTotale)
        
        
        
        #print numberOfClasses(status), cur_varInter, cur_mod, status
        
        #cur_critereGlobal = calculDuScoreAdditionVarianceModularitePondere(status, None, None, cur_varInter, mod=cur_mod)
        modAttr=-2
        cur_critereGlobal, composanteTexte = calculDuScoreAdditionVarianceModularitePondere(numberOfClasses(status), cur_varInter, cur_mod, status, modAttr)
        
        #print(varStatus)
        #print(numberOfClasses(status))
        
        
        
        #assert(cur_critereGlobal!=-1)
        #print "Critere global: "+str(cur_critereGlobal)
        #assert(cur_critereGlobal <100000)
        #assert(cur_critereGlobal > -10000)
        #cur_varInter=varianceInterClasseToutesClasses(status, None)
        
        
        
        
        #print "Variance inter-classes de la partition: "+str(calculDuScoreAdditionVarianceModularitePondere(status, candidat, communauteCandidat))
        if(DEBUG > 4):
            print "Variance inter-classes de la partition: "+str(cur_varInter)
            print "Variance totale : "+str(laVarianceTotale)
            print "Modularite: "+str(cur_mod)
        #cur_critereGlobal = calculDuScoreAdditionVarianceModularitePondere(status, node, dnc, varInterIncrementale, mod=mod_incr)
        #cur_mod+cur_varInter*COEFF_VARIANCE_INTER/laVarianceTotale
        if(DEBUG > 4):
            print "Critere global: "+str(cur_critereGlobal)
        #for i in range(5):
        #        print "pour info "+str( __neighcom(i, graph, status, threshold))
        
        start2 = t.time()
        
        iNodes=0
        #print "Passe: pour tous les noeuds"
        #pprint(graph.nodes())
        #pprint(status.node2com)
        
        for node in graph.nodes() :
                if ((iNodes % 50) == 49 ):
                    print iNodes
                iNodes=iNodes+1
                
                if(DEBUG >1):
                    print "\nEtape noeud="+str(node)
                if (appairing_count<DC_MAX_APPAIRING): 
                # 
                    com_node = status.node2com[node]
                    
                    #degres internes a la communaute
                    degc = status.gdegrees.get(node, 0.)
                    
                    #somme des poids
                    totw = status.total_weight * 2.
                    #communautes voisines
                    neigh_communities = __neighcom(node, graph, status, threshold)
                    
                    
                    
                    #__remove(node, com_node,
                    #        neigh_communities.get(com_node, 0.), status)
                    best_com_mod = [com_node]
                    best_increase_mod = -10
                    
                    #critereGlobalSansFusion = cur_critereGlobal # valeur du critere global si on ne change rien
                    meilleurCritereGlobal=cur_critereGlobal # valeur du critere global si on ne change rien
                    meilleureCommunauteCritereGlobal = com_node
                    
                    #cur_mod = __modularity(status)
                    #print "Avant d'enumerer: "+str(cur_mod)
                    if(DEBUG > 4):
                        pprint(status.node2com)
                    #print "Voisins Candidats: "
                    #for i,j in neigh_communities.iteritems():
                    #    print(str(i))
                    
                    #currentBestValues
                    cvarInterIncrementale=None
                    cx1=None
                    cx2=None
                    cx3=None
                    cx4=None
                    cclassesNumber=None
                    
                    
                    #pprint(neigh_communities.iteritems())
                    
                    
                    # Pour chaque communaute voisine com de node 
                    for com, dnc in neigh_communities.iteritems() :
                        #start = t.time() 
                        if(DEBUG > 2): 
                            print "com_node:"+str(com_node)
                            print str(com)+"<>"+str(dnc)
                        
                        if (com_node == com) : 
                            continue
                        # on met le noeud dans cette communaute et on teste
                        
                        #deg = neigh_communities.get(com, 0.)
                        #__insert(node, com, deg, status)
                        
                        
                        # calcul de la nouvelle modularite avec ce noeud
                        #deltaModularity(m, degWithC, totalDegC, degI)
                        #pprint(status.node2com)
                        #quitte
                        #incr2=deltaModularity(status.total_weight, deg, status.degrees.get(com_node, 0.), status.gdegrees.get(node, 0.))
                        
                        #print "Situation: "
                        #pprint(status.node2com)
                         #print "On bouge "+str(node)+" qui vient de "+str(com_node)+" et qui va pe dans "+str(com)
                        #print "Pour le moment la mod calculee est de "+str(cur_mod)
                        
                        
                        # on enleve de la communaute actuelle com_node
                        incr2=deltaModularity(status.total_weight, neigh_communities.get(com_node, 0.), status.degrees.get(com_node, 0.)-status.gdegrees.get(node, 0.), status.gdegrees.get(node, 0.))
                        #degre_with_nodecom, total degree nodecom
                        #print "incr2: "+str(incr2)
                        # on rejoint la communaute com
                        #pprint( status.total_weight) 
                        #pprint( dnc)
                        #pprint( status.degrees.get(com, 0.))
                        if (com_node == com) : 
                            degreTotDeC = status.degrees.get(com, 0.) - status.gdegrees.get(node, 0.)
                        else:
                            degreTotDeC = status.degrees.get(com, 0.)
                        #pprint( status.gdegrees.get(node, 0.))
                        #incr1=deltaModularity(status.total_weight, dnc,  status.degrees.get(com, 0.), status.gdegrees.get(node, 0.))
                        incr1=deltaModularity(status.total_weight, dnc,  degreTotDeC, status.gdegrees.get(node, 0.))
                        #print "incr1: "+str(incr1)
                        
                        #print (incr1-incr2)
                        #print incr2
                        
                        
                        # la variance dans cette configuration
                        
                        #if varStatus.currentPartition.node2com[node]!=com_node:
                        #    return 
                        
                        looseAssertlt(cur_varInter,laVarianceTotale)
                        #print varStatus.currentPartition.node2com[node],com_node
                        
                        #pprint(status.node2com)
                        #print "status_list"
                        #pprint(status_list)
                        
                        #pprint( varStatus.currentPartition.node2com)
                        assert(varStatus.currentPartition.node2com[node]==com_node)
                        
                        #com_node_varStatusVersion = #varStatus.currentPartition.node2com[node]
                        #com_varStatusVersion = #varStatus.currentPartition.node2com[node]
                        
                        
                        #print com_node_varStatusVersion
                        #print(str(varStatus))
                        #print com
                        #print(str(status))
                        
                        
                        #print "com_node",com_node
                        #pprint(varStatus.nbOfNodesInComs)
                        #assert(varStatus.nbInitialNodes(com_node_varStatusVersion)>0)
                        assert(varStatus.nbInitialNodes(com_node)>0)
                        
                        #print "comSource", com_node_varStatusVersion
                        #print "comTarget", com_varStatusVersion
                        
                        #assert(len(status.node2com) == len(varStatus.comGC))
                        assert(com in varStatus.comGC)
                        
                        (modAttr, varInterIncrementale, x1, x2,  classesNumber) =  varStatus.moveGroup(node, com_node, com, damax=laVarianceTotale )
                        
                        #varStatus.moveGroup(node, com_node_varStatusVersion, com_varStatusVersion , damax=laVarianceTotale )
                        
                        #pprint(varInterIncrementale)
                        #pprint(laVarianceTotale)
                        #exit()
                        looseAssertlt(varInterIncrementale,laVarianceTotale)
                        
                        
                        assert(classesNumber>0)
                        
                        #print node
                        #print com_node
                        #print com
                        #print "on bouge le noeud "+str(node)+" de la communaute "+str(com_node)+" a la communaute "+str(com)+"\nOn obtient une varInterIncr de "+str(varInterIncrementale)
                        #assert(varInterIncrementale)
                        
                        
                        #print "cur_mod vaut avant ajout de l'inc: "+str(cur_mod)
                        mod_incr=cur_mod+incr1-incr2
                        #print "-->",cur_mod,incr1,incr2
                        #pprint(status.node2com)
                        if(DEBUG > 0): 
                            print cur_mod,incr1,incr2
                            print "cur_mod: "+str(cur_mod)
                            print "mod_incr: "+str(mod_incr)
                        
                        #print "Avec la cible la mod calculee est de "+str(mod_incr)
                        
                        #F
                        #__remove(node, com_node,
                        #         neigh_communities.get(com_node, 0.), status)
                        #mod_slow = __modularity(status)
                        #deg = neigh_communities.get(com, 0.)
                        #__insert(node, com_node, deg, status)
                        
                        #print "Avant chgmt:"
                        #pprint(status.node2com)
                        #print str(mod_slow)+"<>"+str(mod_incr)
                        #print "En slow elle sera de "+str(mod_slow)
                        
                        #raw_input("stop")
                        
                        # comparaison avec la modularite calculee a chaque fois
                        #neigh_communities = __neighcom(node, graph, status, threshold)
                        #__remove(node, com_node,
                        #         neigh_communities.get(com_node, 0.), status)
                        
                        #deg = neigh_communities.get(com, 0.)
                        #__insert(node, com, deg, status)
                        
                        #mod_slow = __modularity(status)
                        #print "Comparons: "+str(com)
                        #pprint(status.node2com)
                        #print str(mod_slow)+"<>"+str(mod_incr)
                        
                        #pprint(varInterIncrementale)
                        #pprint(laVarianceTotale)
                        #assert(math.floor(10000000*varInterIncrementale)/10000000<=laVarianceTotale)
                        looseAssertlt(varInterIncrementale,laVarianceTotale)
                        
                        n=len(globTfIdfTab)
                        #print "n",n
                        valeurCalinski = 314105#calinskiGeneratorProbability(classesNumber, varInterIncrementale, laVarianceTotale, n)
                        valeurCritereGlobal, composanteTexte = calculDuScoreAdditionVarianceModularitePondere( classesNumber    , varInterIncrementale, mod_incr, status, modAttr)
                        #print classesNumber    , varInterIncrementale, mod_incr, valeurCritereGlobal#, status
                        #pour  trouver infos veriteTerrain
                        #regularPartition=metaPartition2RegularPartition(status, status_list)
                        #calculAri=metrics.adjusted_rand_score(groundTruthArray, regularPartition)
                        
                        #valeurCritereGlobal=calculAri
                        
                        
                        #assert(abs(mod_incr-mod_slow)<0.001)
                        
                        #neigh_communities = __neighcom(node, graph, status, threshold)
                        #__remove(node, com,
                        #         neigh_communities.get(com, 0.), status)
                        #deg = neigh_communities.get(com_node, 0.)
                        #__insert(node, com_node, deg, status)
                        
                        
                        #print "mod: "+str(mod_incr)+" varInter: "+str(varInterIncrementale)

                        
                        
                        #assert(valeurCritereGlobal < 100000)
                        #assert(valeurCritereGlobal > -10000)
                        
                        
                        #print "Verif ecart var: slow "+str(valeurCritereGlobal)+" <> inc "+str(varInterIncrementale)
                        #assert(  abs(valeurCritereGlobal - varInterIncrementale) < 0.0001  )
                        
                        
                    
                        # la variance dans cette configuration
                        #varInter = valeurCritereGlobal - modNewSituation
                        
                        
                        #varianceInterClasseToutesClasses(status, com)
                        
                        
                        
                        #neigh_communities = __neighcom(node, graph, status, threshold)
                        #__remove(node, com,
                        #         neigh_communities.get(com, 0.), status)
                        #deg = neigh_communities.get(com_node, 0.)
                        #__insert(node, com_node, deg, status)
                        
                        
                        
                        #Amelioration formule louvain
                        cop = [] #copie de la partition courante, que l'on va transformer en tableau pour faire l'evaluation
                        
                        #if (valeurCritereGlobal>meilleurCritereGlobal) or (ECRIRE_LOG_COMPLET==1):
                            
                        
                        
                        
                        if(DEBUG > 2):
                                
                                #pprint(status.node2com)
                                print str(node)+" -> "+str(com)
                                
                                print "#"+str(logCount)
                                print "      ------------------------------------------- "
                                print "      Pour communaute     | "+str(com)
                                # Modularite
                                #poids internes à la communaute
                                #totc = status.degrees.get(com, 0.)
                                #amelioration de modularite (exactement l'amelioration?)
                                #incr =  dnc  - totc * degc / totw
                                #print "Incr: "+str(incr)
                                
                                #vraie inc
                                #incr = -
                        
                                #print "Mod block: "+str(mod_block)
                                
                                #print "incr: "+str(incr)
                                # la modularite dans cette configuration
                                #modNewSituation = cur_mod + incr
                                print "      Modularite incr:    | "+str(mod_incr)
                                #print "      Var inter:          | "+str(varInterClasses)
                                print "      Var interIncrem:    | "+str(varInterIncrementale)
                                
                                # la somme des deux
                                
                                print "      valeurCritereGlobal | "+str(valeurCritereGlobal)
                        
                        
                        
                        #logging
                        if (valeurCritereGlobal>meilleurCritereGlobal) or (ECRIRE_LOG_COMPLET==1):
                        # si le score est meilleur ou que l'on a decide de tout logguer
                            
                            if(DEBUG > 2):
                                print "log "+str(numberOfClasses(status))
                            #print "Calcul de la V-mesure"
                            if groundTruthArray!= None: #si on a une verite terrain
                                #pas top, prend du temps, a ne pas faire aussi souvent
                                regularPartition=metaPartition2RegularPartition(status, status_list)
                                
                                # pour ne pas evaluer systematiquement
                                if ((numberOfClasses(status) % UNE_EVAL_SUR) == 0) or (ari == -2) or (numberOfClasses(status)<50):
                                        # on evalue une fois sur UNE_EVAL_SUR ou on initialise les scores d'evaluation
                                        print "ok mais",(numberOfClasses(status) % UNE_EVAL_SUR)
                                        #if((numberOfClasses(status) % UNE_EVAL_SUR) == 0) or (ari == -2) or (numberOfClasses(status)<50):
                                        cop = status.node2com.copy()
                                        cop[node]=com
                                        # on est dans un cas ou on veut garder la trace et evaluer
                                        print "beginning evaluation"
                                        #on recupere la partition sous-jacente
                                        regularPartition=metaPartition2RegularPartition(status, status_list)
                                        
                                        
                                        assert(len(groundTruthArray)==len(regularPartition))
                                        #pprint(groundTruthArray)
                                        #pprint(regularPartition)
                                        homogeneity, completeness, vMeasure = metrics.homogeneity_completeness_v_measure(groundTruthArray, regularPartition) #.1,.1,.1#
                                        print "ari"
                                        ari=metrics.adjusted_rand_score(groundTruthArray, regularPartition)
                                        print "end evaluation"
                                    
                                regularPartition=metaPartition2RegularPartition(status, status_list)
                                state_saver.add(regularPartition, {"hybridScore":valeurCritereGlobal, "numberOfClasses":numberOfClasses(status)})
                                
                                
                            else : # on a pas de verite terrain
                                ari=-1
                                homogeneity=-1
                                completeness=-1
                                vMeasure=-1
                                
                                
                            
                            # on ecrit le log, eventuellement avec de vieilles valeurs
                            
                            ecrireLog(currentLevel, com_node, node, com, varInterIncrementale, mod_incr, valeurCritereGlobal, laVarianceTotale, valeurCritereGlobal>meilleurCritereGlobal, composanteTexte,valeurCalinski, ari, homogeneity, completeness, vMeasure, numberOfClasses(status))
                        
                            #print "OK"
                            #exit()
                        
                        
                        
                        
                        if valeurCritereGlobal>meilleurCritereGlobal:
                                meilleurCritereGlobal=valeurCritereGlobal
                                meilleureCommunauteCritereGlobal=com
                                cvarInterIncrementale=varInterIncrementale
                                
                                #Pour ne pas refaire le calcul
                                cvarInterIncrementale=varInterIncrementale
                                cx1=x1
                                cx2=x2
                                
                                cclassesNumber=classesNumber
                        
                        #__remove(node, com,
                        #         neigh_communities.get(com, 0.), status)
                        #print "Un bougeage: "+str(t.time() - start)
                        
                    # Fin du test avec tous les voisins du noeud    
                    
                    best_com=meilleureCommunauteCritereGlobal
                    
                    #print "On l'affecte avec "+str(meilleureCommunauteCritereGlobal)
                    neigh_communities = __neighcom(node, graph, status, threshold)
                    modAttr,cur_varInter = varStatus.reallyMoveGroup( node, com_node, best_com, damax=laVarianceTotale, cnewVarianceIntra=cvarInterIncrementale,cx1=cx1, cx2=cx2,  cclassesNumber=cclassesNumber)
                    
                    #with open("reallyMovedNodes.txt", "a") as myfile:
                    #    myfile.write(str(node)+":"+ str(com_node)+" -> "+str( best_com)+"\n")
                    
                    
                    
                    looseAssertlt(cur_varInter,laVarianceTotale)
                    
                    __remove(node, com_node,
                                 neigh_communities.get(com_node, 0.), status)
                    
                    if(DEBUG > 4):
                        print "--Affectation avec la communaute "+str(best_com)
                    
                    
                    
                    deg = neigh_communities.get(best_com, 0.)
                    __insert(node, best_com, deg, status)
                    
                    #assert(varStatus.classesNumber == numberOfClasses(status))
                    #print varStatus.classesNumber , numberOfClasses(status)
                    #print com_node_varStatusVersion
                    #print(str(varStatus))
                    #print com
                    #print(str(status))
                    
                    
                    
                    
                    if best_com != com_node :
                        modif = True
                    
                    
                    #"\n   Variance inter-classe de la partition: "+str(varianceInterClasseToutesClasses(status, None))+"\n   Mod: "+str(cur_mod)
                    #valeurCritereGlobal=calculDuScoreAdditionVarianceModularitePondere(status, None, None)
                    cur_critereGlobal=meilleurCritereGlobal
                    if(DEBUG > 4):
                        print "   Critere global: "+str(meilleurCritereGlobal)
                    #assert(meilleurCritereGlobal <100000)
                    #assert(meilleurCritereGlobal > -10000)
                    new_mod = __modularity(status)
                    assert(new_mod < 1.001)
                    cur_mod = new_mod
                    if(DEBUG > 0):
                        print "Fin: "+str(cur_mod)
                    if(DEBUG > 4):
                        print "   Modularite: "+str(new_mod)
                        print "   Rappel - Variance totale des points: "+str(laVarianceTotale)
                    #if(DEBUG > 2):
                        #print(varStatus)
                    
        print "Fin passe: pour tous les noeuds: temps: "  +str(t.time() - start2)                  
        # si la modularite s'est trop peu amelioree, on arrete la "passe"
        #new_mod = __modularity(status)
        #if new_mod - cur_mod < MIN :
        
        
        
        if cur_critereGlobal - new_crit < MIN :
            if(DEBUG > 0):
                print "------------------------\n        Fin du level\n------------------------"
                #\nVariance inter-classe de la partition: "+str(varianceInterClasseToutesClasses(status, None))+"\nMod: "+str(cur_mod)
                #valeurCritereGlobal=calculDuScoreAdditionVarianceModularitePondere(status, None, None)
                #exit() #corriger ligne au-dessus
                print "Critere global: "+str(valeurCritereGlobal)
                #assert(valeurCritereGlobal < 100000)
                #assert(valeurCritereGlobal > -10000)
                print "Rappel - Variance totale des points: "+str(laVarianceTotale)
                
                #pprint( recalculDesAttributs(status) )
            #pprint(status.node2com)
            #print(str(varStatus))
            #exit()
            #raw_input("xxx")
            return cur_critereGlobal
                
        if(DEBUG > 0):
            print "/------------------------------------------\n-----------------------------------PASS DONE\n\\-------------------------------------------"
    print "<<<<<<<<<<<Cnul>>>>>>>>>>>"
    #raw_input("Cnul")
    #assert(not (math.isnan(cur_critereGlobal)))
    return cur_critereGlobal

def numberOfClasses(status):
    """
    Retourne le nombre de classe en ce moment dans la structure pour la modularite
    """
    return len(set(status.node2com.values()))


def deltaModularity(m, degWithC, totalDegC, degI):
    #on mult par 2 le degWithC donc on enleve le 2 au 2m
    #print "m: "+str(m)
    #print "deg avec C: "+str(degWithC)
    #print "deg total de C: "+str(totalDegC)
    #print "deg de i: "+str(degI)
    #degWithC=degWithC*2
    #print autre formule
    return (   2*degWithC/(2*m)  -  ((totalDegC+degI)/(2*m))**2 )   +  (totalDegC/(2*m))**2  +   (degI/(2*m))**2
    
     
    




class Status :
    """
    To handle several data in one struct.

    Could be replaced by named tuple, but don't want to depend on python 2.6
    """
    node2com = {}
    total_weight = 0
    internals = {}
    degrees = {}
    gdegrees = {}
    def __init__(self) :
        self.node2com = dict([])
        self.total_weight = 0
        self.degrees = dict([])
        self.gdegrees = dict([])
        self.internals = dict([])
        self.loops = dict([])
    def __str__(self) :
        return ("node2com : " + str(self.node2com) + " degrees : "
            + str(self.degrees) + " internals : " + str(self.internals)
            + " total_weight : " + str(self.total_weight))

    def copy(self) :
        """Perform a deep copy of status"""
        new_status = Status()
        new_status.node2com = self.node2com.copy()
        new_status.internals = self.internals.copy()
        new_status.degrees = self.degrees.copy()
        new_status.gdegrees = self.gdegrees.copy()
        new_status.total_weight = self.total_weight

    def init(self, graph, part = None) :
        #pprint(part)
        """Initialize the status of a graph with every node in one community"""
        count = 0
        self.node2com = dict([])
        self.total_weight = 0
        self.degrees = dict([])
        self.gdegrees = dict([])
        self.internals = dict([])
        self.total_weight = graph.size(weight = "weight")
        if part == None :
            for node in graph.nodes() :
                self.node2com[node] = count
                deg = float(graph.degree(node, weight = "weight"))
                self.degrees[count] = deg
                self.gdegrees[node] = deg
                #self.total_weight += deg
                #boucles / on double le poids du coup sur les boucles ?
                self.loops[node] = float(graph.get_edge_data(node, node, {"weight":0}).get("weight", 1))
                self.internals[count] = self.loops[node]
                count = count + 1
        else :
            for node in graph.nodes() :
                com = part[node]
                self.node2com[node] = com
                deg = float(graph.degree(node, weight = "weight"))
                
                self.degrees[com] = self.degrees.get(com, 0) + deg
                self.gdegrees[node] = deg
                #self.total_weight += deg
                inc = 0.
                for neighbor, datas in graph[node].iteritems() :
                    #pprint(graph[node])
                    #print "node",node, neighbor
                    #pprint(datas)
                    weight = datas.get("weight", 1)
                    #si la classe du voisin est la meme que le noeud examine
                    if part[neighbor] == com :
                        
                        if neighbor == node :
                            inc += float(weight)
                        else :
                            # on ajoute la moitie du poids de l'arete
                            inc += float(weight) / 2.
                self.internals[com] = self.internals.get(com, 0) + inc



def __neighcom(node, graph, status, threshold) :
    """
    Compute the communities in the neighborood of node in the graph given
    with the decomposition node2com
    """
    weights = {}
    for neighbor, datas in graph[node].iteritems() :
        if (neighbor != node) :
            #
            #if distance(neighbor, node)> threshold:
            #
                
                weight = datas.get("weight", 1)
                neighborcom = status.node2com[neighbor]
                weights[neighborcom] = weights.get(neighborcom, 0) + weight
    return weights

def __remove(node, com, weight, status) :
    """ Remove node from community com and modify status"""
    status.degrees[com] = ( status.degrees.get(com, 0.)
                                    - status.gdegrees.get(node, 0.) )
    status.internals[com] = float( status.internals.get(com, 0.) -
                weight - status.loops.get(node, 0.) )
    status.node2com[node] = -1

def __insert(node, com, weight, status) :
    """ Insert node into community and modify status"""
    status.node2com[node] = com
    status.degrees[com] = ( status.degrees.get(com, 0.) +
                                status.gdegrees.get(node, 0.) )
    status.internals[com] = float( status.internals.get(com, 0.) +
                        weight + status.loops.get(node, 0.) )

def __modularity(status) :
    """
    Compute the modularity of the partition of the graph
    """
    links = float(status.total_weight) # = m
    if(DEBUG > 4):
        print "Links(somme des poids): " + str(links)
    result = 0.
    for community in set(status.node2com.values()) :
        #le in_degree est la moitie du nombre d'aretes (ou du poids qui 
        #leur est affecte) qui ont les deux extremites dans la classe examinee
        in_degree = status.internals.get(community, 0.)
        degree = status.degrees.get(community, 0.)
        if links > 0 :
            result = result + in_degree / links - ((degree / (2*links))**2)
            #print  str(in_degree) +"/"+ str(links )+"-(("+str(degree) +"/ (2*"+str(links)+"))**2)"
            #result = result + in_degree / m - ((degree / (2*m))**2)
            # in_degree : le double du nb d'aretes dans la communaute
            #print "in:"+str(in_degree)+" adj:"+str(degree)+" modIntern: "+str(in_degree / links - ((degree / (2.*links))**2))
            if(DEBUG > 4):
                print "in:"+str(in_degree)+" adj:"+str(degree)+" modIntern: "+str(in_degree / links - ((degree / (2.*links))**2))
                
    if(DEBUG > 4):
        print "Result: "+str(result)
    #exit()
    #pprint(status.node2com)
    assert(result < 1)
    return result
    
def metaPartition2RegularPartition(status,dendogram):
    #pprint(dendogram)
    if len(dendogram)==0:
        part=dicToArraySimple(status.node2com, currentGlobAuthorIndex)
        return part
        
    else:
        level = len(dendogram) - 1
        partition = dendogram[0].copy() #on prend la liste originale des noeuds
        for index in range(1, level + 1) :
            for node, community in partition.iteritems() :
                partition[node] = dendogram[index][community]
        #pprint(partition)
        for node, community in partition.iteritems() :
            partition[node] = status.node2com[community]
            #print str(node)+" -> "+str(status.node2com[community])
        #assert(partition['3']==3)
        #print "pprint"
        #pprint(partition)
        #pprint(globAuthorIndex)
        #print "Dendogram length"
        #print len(dendogram)
        #pprint(status.node2com)
        #pprint(dendogram)
        part=dicToArraySimple(partition, globAuthorIndex)
        return part
        #return partition
        
        
        
        #under_p = partition_at_level(dendogram, len(dendogram) - 1 )
        #status.node2com
        #return status.node2com

def looseAssertlt(v1,v2):    
    #pprint(v1)
    #pprint(v2)
    #assert(v1<=v2)
    #assert(math.floor(100000000*v1)/100000000<=v2)
    pass
    

def __main() :
    """Main function"""
    from testsCommunityThreshold import launchTests, testWithExample, testWithBiggerExample
    testWithExample()
    exit()
    try :
        #filename = sys.argv[1]
        testWithExample()
        #testWithBiggerExample()
    except (IndexError, IOError):
        print "Usage : ./community"
        print "Performs communauty detection for a small (5 nodes) graph"


if __name__ == "__main__" :
    __main()


