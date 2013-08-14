#    Copyright (C) 2013 by
#    David Combe (david.combe@gmail.com), Universite Jean Monnet, Saint-Etienne
#    All rights reserved.
#    BSD license. 


from __future__ import division

from pprint import pprint 
import networkx as nx # pour les exemples
import sys
import random 
import math

#import modAttrEsp as ma


from scipy.spatial.distance import cosine, euclidean
from scipy import zeros, ones
from copy import deepcopy
import numpy as np


from varIncrement import *

from  scipy.stats.mstats import gmean
from  scipy.stats.mstats import hmean
from numpy import mean
#Que manque t il ?
# APartition : une partition

#global logTruc
#logTruc = open("logTruc.csv", 'w')

#def logguerUnTruc(truc):
#    global logTruc
#    
#    logTruc.write(str(truc)+"\n")


class APartition :
    """
    Classe qui contient un mapping entre les groupes et les metagroupes
    le dictionnaire de la partition dans les 2 sens
    le compte des noeuds dans les metagroupes
    """
    node2com={}
    com2node={}
    

    def __init__(self,partition=None, nbComs=None) :
        
        if partition==None:
            #partition={}
            #count=0
            #for i in self.node2com:
            #    partition[i] = count
            #    count=count+1
            assert(nbComs != None)
            partition={}
            count=0
            for i in range(0,nbComs):
                partition[i] = count
                count=count+1
        self.node2com = partition
        self.com2node=self.initDicInverse()
        
        
    def initDicInverse(self):
        """
        cree le dictionnaire inverse
        """
        dicinv={}
        for (k,v) in self.node2com.items(): 
                a=dicinv.get(v, set())
                a.add(k)
                dicinv[v]=a
                #print "k:"+str(k)+"  v:"+str(v)
        #print "dic_inv"
        #pprint(dicinv)
        return dicinv
        
    def dicInverse(self):
        return self.com2node
    
    def update(self,node, com):
        self.node2com[node]=com
        self.com2node[com]=node
    
    def __hash__(self):
        dicinv=self.dicInverse()
        res=0
        for k in dicinv:
            res=res+(','.join(dicinv[k])).__hash__()
        return res
        
        
    def numberOfNodesInCom(self,com):
        #pprint(self.com2node)
        return len(self.com2node[com])
        
        
    def __eq__(self, b):
        """
        Deux partitions sont identiques si...
        """
        if b==None:
            return False
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


class metaNodesHistory:
    '''
    hierarchie des groupements
    
    '''
    partitions=[] # une liste de APartition
    numberOfNodes={}
    
    def __init__(self,initialNode2com): # a changer pour imposer une partition de depart
        self.partitions.append(APartition(initialNode2com))
        #numberOfNodes=zeros([len(initialNode2com)])
        for i in initialNode2com:
            self.numberOfNodes[i]=1
        
    
    def getLatestPartition(self):
        assert(len(self.partitions)>0)
        return self.partitions[len(self.partitions)-1]
    
    def nodesInACom(self, com):
        '''
        Donne les noeuds presents dans nimporte quelle partition de n'importe quel niveau.
        '''
        i=len(self.partitions)-1
        while i>1:
            val=self.partitions[i].node2com[com]
            i=i-1
        return val 
    
    
    def infererLaPartitionDeTousLesSommets(self, currentPartition=None):
        """
        Combine les partitions de tous les niveaux pour obtenir la partition "brute"
        """
        
        #print(self)
        
        
        tabAffectation = np.zeros([len(self.partitions[0].node2com)])
        
        for i in range(len(self.partitions[0].node2com)):
            tabAffectation[i] = self.partitions[0].node2com[i]
        
        
        
        n=1
        #n=len(self.partitions)-1
        while n<len(self.partitions): #>1:
            for i in range(len(self.partitions[0].node2com)):
            
                #pprint(self.partitions)
                tabAffectation[i] = self.partitions[n].node2com[tabAffectation[i]]
                
            #val=self.partitions[n].node2com[com]
            n=n+1
        #pprint(tabAffectation)
            
        
        
        if currentPartition != None :
            
            for i in range(len(self.partitions[0].node2com)):
                tabAffectation[i] = currentPartition[tabAffectation[i]]
            
        #pprint(tabAffectation)
            
        
        
        return tabAffectation
        
    
    def nbNodesInAGroup(self, group):
        #le nombre de noeuds d'un certain group dans la derniere partition enregistree
        
        #pprint(self.partitions)
        #return self.partitions[len(self.partitions)-1].numberOfNodesInCom(com)
        return self.numberOfNodes[group]



    def ajouterUnePartition(self, p):
        """
        dans la hierarchie
        """
        #print str(self)
        #print "corriger numberOfNodes"
        
        #pprint(self.numberOfNodes)
        #pprint(p)
        
        ap = APartition(p)
        inverse = ap.dicInverse()
        #pprint(inverse)
        nbNodes={}
        for c in inverse:
            for n in inverse[c]:
                #print c,n
                nbNodes[c]=nbNodes.get(c,0) + self.nbNodesInAGroup(n)
        #calcul du dict de nombre de noeuds par communaute de p
        
        #for i in p:
        #    print "i",i
        #    comI=p[i]
        #    print "comI",comI
        #    print nbNodes.get(comI,0) 
            #print self.numberOfNodes[len(self.numberOfNodes)-1][comI]
        #    nbNodes[comI] = nbNodes.get(comI,0) + self.partitions[-1][comI]
        #    pprint(self.numberOfNodes)
            #! le str ci-dessous est pas logique
            #nbNodes[comI] = nbNodes.get(comI,0) + self.numberOfNodes[str(comI)]
            
        #pprint(nbNodes)
        self.numberOfNodes = nbNodes
        #self.numberOfNodes.append(nbNodes)
        self.partitions.append(ap)
        

    def __str__(self):
        out= "metaHistory\nlistes des partitions dans l'historique\n"
        j=0
        for i in self.partitions:
            out = out + "Partition "+str(j)+"\n"
            out=out+str(i)+"\n"
            
        out=out+"numberOfNodes\n"+str(self.numberOfNodes)+"\n"
        
        
        
        return out+"\nFin de l'historique\n"



class MetaVarIncrement :
    '''
    Meta: les noeuds deviennent des groupes,
    les communautes restent des communautes.
    '''
    #underlyingGroups={} # de 1 noeud vers son group
    firstGlobalTfIdfTab=None
    hierarchy=None
    
    # pour un group en particulier
    groupVariance={}
    #groupMoyenne={}
    groupGC={}
    
    comGC={}
    comVar={}
    
    #moyenne des variances des communautes
    VarComsMean=None
    
    
    #la partition courante
    currentPartition=None
    classesNumber=0
    
    currentVarInter=0L
    currentModAttr=0L
    
    #nbElements=[]
    
    #centre de gravite et variance de chaque groupe
    #groupVariance=[]
    
    
    
    #centre de gravite et variance de chaque communaute
    # plutot cdg
    varInterSubGravityCenter=[]
    
    centreDeGraviteTotal=[]
    
    
    #currentNode2Com={}
    
    euclideanDistanceOfComsToGlobGC=[]
    
    
    nbInitialNodesInGroups={}
    
    # preciser ce truc
    nbOfNodesInComs=[]
    
    
    def __init__(self, underlyingGroups, globalTfIdfTab, partition, globAuthorIndex, hierarchy ) :
        '''
        On definit les groupes.
        '''
        
        #self.underlyingGroups=underlyingGroups
        
        self.hierarchy=hierarchy
        #On enregistre la variance et les GC de tous les groupes
        
        self.firstGlobalTfIdfTab=globalTfIdfTab
        
        #on cree la partition des groupes en communautes discrete
        self.currentPartition=APartition(partition=partition) 
        # a changer si on veut imposer une partition de depart
        #pprint(self.currentPartition.node2com)
        
        #pprint(globalTfIdfTab)
        #pprint(len(globAuthorIndex))
        # calcul des cdg
        
        self.VarComsMean=0
        
        self.varInterSubGravityCenter=zeros([len(globAuthorIndex), len(globalTfIdfTab[0])])
        #pprint(self.varInterSubGravityCenter)
        

        for i in range(len(globAuthorIndex)):
            #pprint(globalTfIdfTab[i])
            #pprint(self.varInterSubGravityCenter[i])
            self.varInterSubGravityCenter[i]= globalTfIdfTab[i].copy()
        
        
        # calcul des GC
        #pprint(partition)
        for i in partition:
            #print i
            #pprint(globAuthorIndex[i])
            #pprint(globalTfIdfTab[globAuthorIndex[i]])
            #self.groupGC[i] = globalTfIdfTab[globAuthorIndex[i]]
            self.groupGC[i] = globalTfIdfTab[i]
        #pprint(self.groupGC)
        
        #pprint(partition)
        #pprint(globalTfIdfTab)

        
        
        # calcul des variances
        # a verifier
        # structure pour isoler les valeurs presentes dans chaque groupe :
        
        
        #lesValeurs[attribut]{classe}[sommets]
        
        
        lesValeurs = []
        for j in range(len(globalTfIdfTab[0])):
            #pprint(lesValeurs)
            lesValeurs.append({})
            for k, v in partition.items():
                #pprint(lesValeurs)
                #print j,k,v
                #print lesValeurs[j].get(v,[])
                tmp= lesValeurs[j].get(v,[])
                tmp.append(globalTfIdfTab[k][j])
                lesValeurs[j][v] = tmp
        
        #pprint(lesValeurs)
        #initialisation des variances et des moyennes
        for i in set(partition.values()): # pour chaque classe
            #self.comVar[i] = zeros([len(globalTfIdfTab[0])])
            self.comGC[i] = zeros([len(globalTfIdfTab[0])])
            for j in range(len(globalTfIdfTab[0])): #pour chaque attribut
                #print lesValeurs[j][i]
                #self.comVar[i][j] = np.var(lesValeurs[j][i])
                self.comGC[i][j] = np.mean(lesValeurs[j][i])
            #pprint(self.comVar)
            #pprint(self.comGC)
            
        
        
        # for k,v in partition.items():
            # #print k
            # #print v
            # #print "com "+str(v) +" -> "+str(globalTfIdfTab[globAuthorIndex[k]])
            
            # #self.comGC[v]=globalTfIdfTab[globAuthorIndex[k]]
            # self.comGC[v]=globalTfIdfTab[k]
        # #pprint(self.comGC)
        
        
        
        
        
        
        
        
        for i in self.currentPartition.node2com:
            #print i
            #print ">"
            #print globalTfIdfTab[globAuthorIndex[i]]
            #print globalTfIdfTab[i]
            #self.groupGC[i] = globalTfIdfTab[globAuthorIndex[i]]
            self.groupGC[i] = globalTfIdfTab[i]
            #self.comMoy[i] = globalTfIdfTab[globAuthorIndex[i]]
        
        #pprint(self.groupGC)
        
        
        self.classesNumber = len(set(partition.values()))  #globalTfIdfTab)
        
        
        
        #initialiser node2com a la partition discrete
        
        #self.currentNode2Com=APartition() 
        
        try:
            nbCommunities=len(hierarchy.getLatestPartition().node2com.keys)
        except:
            nbCommunities=len(globAuthorIndex)
        nbCommunities=len(set(partition.values()))
        self.nbOfNodesInComs=zeros([nbCommunities],int) # a changer si on veut imposer une partition de depart
        
        #compte combien de sommets dans chaque classe
        for k,v in lesValeurs[0].items():
            
            #lesValeurs[attribut]{classe}[sommets]
            
        
            self.nbOfNodesInComs[k] = len(v)
        #pprint(self.nbOfNodesInComs)    
        #exit()
        
        self.centreDeGraviteTotal=self.moyenne_points(range(len(globalTfIdfTab)))
        
        #self.nbInitialNodesInGroups=ones(len(globAuthorIndex), dtype=int)
        self.nbInitialNodesInGroups={}
        for i in partition:
            self.nbInitialNodesInGroups[i]=1
        
        #print "nbCommunities"
        #print nbCommunities
        #for i in range(0,nbCommunities):
        #    self.nbOfNodesInComs[i] = 1
            
        nbAttr = len(globalTfIdfTab[0])
            
        for i in self.currentPartition.node2com:
        #    self.groupMoyenne = globalTfIdfTab[globAuthorIndex[i]]
            self.groupVariance[i]=zeros([nbAttr])
            
        #pprint(self.groupMoyenne)
        
        #!
        #print "euclideanDistanceOfComsToGlobGC a initialiser"
        #euclideanDistanceOfComsToGlobGC=
        self.euclideanDistanceOfComsToGlobGC=zeros([len(globAuthorIndex)])
        #for i in self.currentPartition.node2com:
        #    self.euclideanDistanceOfComsToGlobGC[i]=0
        #distGCglobToNewSourceCom = euclidean(self.centreDeGraviteTotal,newGCcomSource)
        
        #initialisation de la var inter
        # la variance de la distance des cdg au cdg global
        
        # a changer si on peut forcer la partition
        out=0
        #print "cdgTotal"
        #pprint(self.centreDeGraviteTotal)
        #pprint(globAuthorIndex)
        


        # for i in globAuthorIndex.values():
            # #print "------"
            # #print (euclidean(globalTfIdfTab[i],self.centreDeGraviteTotal))
            # #print (euclidean(globalTfIdfTab[i],self.centreDeGraviteTotal))**2
            # #print (1.0/ len(globAuthorIndex) )*(euclidean(globalTfIdfTab[i],self.centreDeGraviteTotal))**2
            # out = out+(1.0/ len(globAuthorIndex) )*(euclidean(globalTfIdfTab[i],self.centreDeGraviteTotal))**2
            # self.euclideanDistanceOfComsToGlobGC[i] = (euclidean(globalTfIdfTab[i],self.centreDeGraviteTotal))
        
        
        # #initialisation des variances et des moyennes
        # for i in partition.values(): # pour chaque classe
            # #self.comVar[i] = zeros([len(globalTfIdfTab[0])])
            # #self.comGC[i] = zeros([len(globalTfIdfTab[0])])
            # for j in range(len(globalTfIdfTab[0])): #pour chaque attribut
                # #print lesValeurs[j][i]
                # self.comVar[i][j] = np.var(lesValeurs[j][i])
                # self.comGC[i][j] = np.mean(lesValeurs[j][i])
            # pprint(self.comVar)
            # pprint(self.comGC)
        
            
        
        
        # pour toutes les communauts
        for i in set(partition.values()): # pour chaque classe
            #on calcule la distance du centre au centre de gravite total
            out = out+(len(lesValeurs[0][i])/ len(globAuthorIndex) )*(euclidean(self.comGC[i],self.centreDeGraviteTotal))**2
            #print "_",out
          
        
        self.currentVarInter = out
        self.verifCoherence()
        
        #print self
        #exit()
        
    
    def verifCoherence(self):
        '''
        Verifie que le decompte de nbOfNodesInComs est correcte vis a vis du nombre total de noeuds,
        que la variance inter est positive.
        '''
        out = 0
        
        #print "self.nbOfNodesInComs"
        #pprint(self.nbOfNodesInComs)
        
        for v in self.nbOfNodesInComs:
            out=out + v
        #print "len(self.firstGlobalTfIdfTab)",len(self.firstGlobalTfIdfTab)
        assert(out == len(self.firstGlobalTfIdfTab))
        
        #que la variance inter est positive
        #assert(self.currentVarInter>0)
        
        """def validatePartition(self):
        #end of phase 2
        hierarchy.ajouterUnePartition(currentNode2Com)"""
    
    def nbInitialNodes(self,com):
        '''
        Le nb d'atomes dans chaque communaute.
        '''
        #print str(com)+"<-com"
        #pprint(self.nbOfNodesInComs)
        return self.nbOfNodesInComs[com]
    
    #ne fonctionne que sur les premiers points
    def moyenne_points(self,indices): 
        '''
        Calcule le centre de gravite d'un ensemble de points,
        a savoir la moyenne des differents attributs
        Les vecteurs sont recuperes dans la variable globale globTfIdfTab
        '''
        nb_att=len(self.firstGlobalTfIdfTab[indices[0]])
        tab=zeros([nb_att])
        
        #pprint(tab)
        #pprint(indices)
        
        for i in indices:
            for j in range(nb_att):
                tab[j]=tab[j]+self.firstGlobalTfIdfTab[i][j]
        for j in range(nb_att):
            tab[j]=tab[j]/len(indices)
        return tab
        
     #def getMeanOfVariancesOfGroups():
        
        
     #   return
        
        
    def meanOfVariancesOfGroups(self,typeMean="usual"):
        if(typeMean=="usual"):
            #print "usual"
            return mean(self.groupVariance.values())
        elif(typeMean=="harmonic"):  #That is: n / (1/x1 + 1/x2 + ... + 1/xn)
            #print "harmonic"
            return hmean(self.groupVariance.values())
        elif(typeMean=="geometric"):  #That is: n-th root of (x1 * x2 * ... * xn)
            #print "geometric"
            return gmean(self.groupVariance.values())
        #groupVariance={}
        print "Error: meanOfVariancesOfGroups"
        exit()
        return None
        
    def inferNewMeanOfVariancesOfGroups():
        self.VarComsMean
        pass
        
        
    def moveGroup(self, group, comSource, comTarget, damax=-1):
        '''
        toutes les valeurs sont celles "initialement" presentes
        #On renvoie la varInter qu'il y aurait dans cette configuration, le cdg source et le cdg cible, et le nombre de classes
        #target est com2
        '''
        #self.verifCoherence()
        #print "previous vic"
        #pprint( self.currentVarInter)
        #pprint(damax)
        looseAssertlt(self.currentVarInter,damax)
        
        #print "Debut moveGroup"
        
        if (comSource == comTarget ):
            #(newVarianceIntra, newGCcomSource, newVarcomSource, newGCcomCible, newVarcomCible , classesNumber)
            
            #print "comSource=comTarget -> deviation"
            #return self.currentVarInter,    self.comGC[comSource]     ,self.comVar[comSource] ,  self.comGC[comTarget]  ,self.comVar[comTarget] ,self.classesNumber
            return self.currentModAttr, self.currentVarInter,    self.comGC[comSource]    ,  self.comGC[comTarget],  self.classesNumber
            
            
        #print "###"    , group
        #pprint(self.currentPartition.node2com)
        #pprint(comSource)
        
        assert(self.currentPartition.node2com[group]==comSource)
        
        
        
        
        nbNodesGroup = self.hierarchy.nbNodesInAGroup(group)
        
        
        #initialement
        #print "comSource"
        #print comSource
        #pprint(self.groupGC)
        
        #print self.comGC
        #pprint(self.comGC)
        
        #print "gcComSource"
        gcComSource=self.comGC[comSource]
        #print gcComSource
        
        
        #print "Debug"
        #pprint(self.comGC)
        #pprint( self.groupGC) 
        #print comTarget
        
        #print "gcComTarget"
        try:
            gcComTarget=self.comGC[comTarget]
        except:
            pprint(self.comGC)
            print comTarget
            
            print(self, group, comSource, comTarget, damax)
            print(str(self))
            exit()
        #print gcComTarget
        #print "varSource"
        #varSource= self.comVar[comSource]
        #print varSource
        #print "varCible"
        #print "comTarget",comTarget
        #pprint(self.comVar)
        #print self, group, comSource, comTarget
        #varCible= self.comVar[comTarget]
        #print varCible
        
        
        nbNodesSource= self.nbInitialNodes(comSource) 
        nbNodesTarget= self.nbInitialNodes(comTarget)
        assert(nbNodesSource>0)
        assert(nbNodesTarget>0)
        
        
        assert(nbNodesTarget==nbNodesTarget)
        
        #print "&&&&",nbNodesSource
        
        
        #print self.centreDeGraviteTotal,gcComSource
        #anciennes distances au CDG global de comSource et comTarget
        distGCglobToAncSourceCom = euclidean(self.centreDeGraviteTotal,gcComSource)
        
        #print distGCglobToAncSourceCom
        #print distGCglobToAncSourceCom**2/4
        
        
        
        #exit()
        #print self.centreDeGraviteTotal,gcComTarget
        distGCglobToAncTargetCom = euclidean(self.centreDeGraviteTotal,gcComTarget)
        
        # processing
        #print ">>>>>"
        
        
        
        #print "-----****-----moveGravityCentersVersionGroups"
        #print group, nbNodesGroup, gcComSource, nbNodesSource, gcComTarget, nbNodesTarget
        (newGC1, newGC2) = self.moveGravityCentersVersionGroups(group, nbNodesGroup, gcComSource, nbNodesSource, gcComTarget, nbNodesTarget)
        #print newGC1, newGC2
        #print "-----****-----"
        
        #pprint(self.currentPartition)
        #pprint(self.currentPartition.node2com)
        
        
        #pprint(group)
        
        #pprint(self.groupMoyenne)
        
        
        
        # calcul intermediaires
        # se rapporte au group deplace
        moyB= self.groupGC[group]
        
        #pprint(self.groupVariance)
        #pprint(group)
        
        #print "group",group
        #print "self.groupVariance"
        #pprint(self.groupVariance)
        #varB = self.groupVariance[group]
        
        #pprint(self.nbInitialNodesInGroups)
        
        
        nB = self.nbInitialNodesInGroups[group]
        
        moyGroup = self.groupGC[group]
        #varGroup = self.groupVariance[group]
        
        
        moyComTarget = self.comGC[comTarget]
        #varComTarget = self.comVar[comTarget]
        
        
        # ce que l'on veut calculer
        
        newGCcomCible = newGC1 # l'ancien GC -1/(dif entre ancien GC et attribut du group deplace ) pondere apr le nombre de noeuds dedans
        # la var de la com source - la var du group bouge
        #assert(varB!=None)
        #newVarcomSource = self.__varianceDifferenceTwoSets(gcComSource, varSource, moyB, varB, nbNodesSource, nB)
        newGCcomSource = newGC2
        
        #print moyGroup, varGroup, moyComTarget, varComTarget, nbNodesGroup, nbNodesTarget
        #newVarcomCible = self.varianceUnionTwoSets(moyGroup, varGroup, moyComTarget, varComTarget, nbNodesGroup, nbNodesTarget)


        
        
        
        
        
        #pprint(varComTarget)
        #raw_input("pprint(varComTarget)")
        #logguerUnTruc(varComTarget)
        
        #pprint(newVarcomCible)
        #if(newVarcomCible[0] != 0):
        #    raw_input("newVarcomCible != 0")
        
        # calcul des 2 nouveaux gc
        #gc1dest, gc2source = self.moveGravityCentersVersionGroups( self.varInterSubGravityCenter[comTarget]      ,      group    , (len(self.currentPartition.dicInverse()[comTarget]))  , self.varInterSubGravityCenter[comSource], (len(self.currentPartition.dicInverse()[comSource])))
        #gc1dest, gc2source = self.moveGravityCentersVersionGroups(  group    , self.comGC[comTarget]      ,      (len(self.currentPartition.dicInverse()[comTarget]))  , self.varInterSubGravityCenter[comSource], (len(self.currentPartition.dicInverse()[comSource])))
        # calcul des 2 nouveaux variance
        
        
        #print "ancGCcomSource->newGCcomSource"
        #print gcComSource, newGCcomSource
        
        #print "Doit disparaitre !"
        
        #print "ancVarcomSource->newVarcomSource"
        #print varSource, newVarcomSource
        
        
        #print "ancGCcomCible->newGCcomCible"
        #print gcComSource, newGCcomCible
        
        #print "varSource->newVarcomCible"
        #print varSource, newVarcomCible
        
        
        n=len(self.firstGlobalTfIdfTab)
        
        #logguerUnTruc(newVarcomCible)
        
        #distGCglobToAncSourceCom = euclidean(self.centreDeGraviteTotal,gcComSource)
        #distGCglobToAncTargetCom
        
        # on calcule la distance au GC des nouvelles com
        
        distGCglobToNewTargetCom = euclidean(self.centreDeGraviteTotal,newGCcomCible)

        #print self.centreDeGraviteTotal,newGCcomCible
        #print distGCglobToNewSourceCom, distGCglobToNewTargetCom
        
        
        if (nbNodesSource-nbNodesGroup==0):
            newComposanteComSource=0
            #print "disparait"
        else:
            #print "disparaitPas"
            distGCglobToNewSourceCom = euclidean(self.centreDeGraviteTotal,newGCcomSource)
            #print distGCglobToNewSourceCom
            #print self.centreDeGraviteTotal,newGCcomSource
            newComposanteComSource=((nbNodesSource-nbNodesGroup)/n)*distGCglobToNewSourceCom**2
            #print nbNodesSource,nbNodesGroup,n,distGCglobToNewSourceCom
            #print distGCglobToNewSourceCom
            assert(newComposanteComSource>=0)

        
        
        #print "on calcule la distance au GC des nouvelles com"
        #print distGCglobToNewTargetCom
        
        #print "previous vic"
        #pprint( self.currentVarInter)
        #print "newVarianceIntra"
        # la ponderation des distances au carre du nombre de noeuds concernes par un cdg
        
        #print nbNodesTarget,n
        
        assert(self.currentVarInter>=0)
        assert((nbNodesSource/n)*distGCglobToAncSourceCom**2>=0)
        assert((nbNodesTarget/n) * distGCglobToAncTargetCom**2>=0)
        assert(newComposanteComSource>=0)
        
        assert(((nbNodesTarget+nbNodesGroup)/n)*distGCglobToNewTargetCom**2 >=0)
        
        
        #newVarianceIntra= self.currentVarInter - (nbNodesSource/n)*distGCglobToAncSourceCom**2 -(nbNodesTarget/n) * distGCglobToAncTargetCom**2 + newComposanteComSource + ((nbNodesTarget+nbNodesGroup)/n)*distGCglobToNewTargetCom**2
        newVarianceIntra= self.currentVarInter - (float(nbNodesSource)/float(n))*float(distGCglobToAncSourceCom)**2.0 -(float(nbNodesTarget)/float(n)) * float(distGCglobToAncTargetCom)**2.0 + float(newComposanteComSource) + ((float(nbNodesTarget)+float(nbNodesGroup))/float(n))*(float(distGCglobToNewTargetCom)**2.0)
        
        #logguerUnTruc("+"+str(self.currentVarInter)+"-"+str( (float(nbNodesSource)/float(n))*float(distGCglobToAncSourceCom)**2)+"-"+str( (float(nbNodesTarget)/float(n)) * float(distGCglobToAncTargetCom)**2 )+"+"+str( float(newComposanteComSource) )+"+"+str(((float(nbNodesTarget)+float(nbNodesGroup))/float(n))*float(distGCglobToNewTargetCom)**2))
        #logguerUnTruc(newVarianceIntra)
        
        #print "la variance etait de ",self.currentVarInter
        #print "la com source valait", (float(nbNodesSource)/float(n))*float(distGCglobToAncSourceCom)**2
        #print "la com cible valait", (float(nbNodesTarget)/float(n)) * float(distGCglobToAncTargetCom)**2
        #print "la com source vaudra", str( float(newComposanteComSource) )
        #print "la com cible vaudra", str(((float(nbNodesTarget)+float(nbNodesGroup))/float(n))*float(distGCglobToNewTargetCom)**2)
        
        
        #print "+",self.currentVarInter,"-", (float(nbNodesSource)/float(n))*float(distGCglobToAncSourceCom)**2,"-", (float(nbNodesTarget)/float(n)) * float(distGCglobToAncTargetCom)**2 ,"+", float(newComposanteComSource) ,"+", ((float(nbNodesTarget)+float(nbNodesGroup))/float(n))*float(distGCglobToNewTargetCom)**2
        
        #pprint(("+",self.currentVarInter,"-", (nbNodesSource/n)*distGCglobToAncSourceCom**2,"-", (nbNodesTarget/n) * distGCglobToAncTargetCom**2 ,"+", newComposanteComSource ,"+", ((nbNodesTarget+nbNodesGroup)/n)*distGCglobToNewTargetCom**2))
        
        #print nbNodesTarget,nbNodesGroup,n,distGCglobToNewTargetCom
        #pprint(newVarianceIntra)
        #pprint(damax)
        
        
        #looseAssertlt(newVarianceIntra,damax)
        #!
        #print "pas bien !"
        if False: #newVarianceIntra>damax:
            newVarianceIntra=damax
            
            print "------------------------\nresultats de moveGroup"
            #print (newVarianceIntra, newGCcomSource, newVarcomSource, newGCcomCible, newVarcomCible , classesNumber)
            # for i in self.comVar:
                # print i, self.comVar[i], self.comGC[i], self.nbOfNodesInComs[i]
            
            # pprint(self.comVar[comSource])
            # pprint(self.comVar[comTarget])
            print comSource, comTarget
            pprint(newVarianceIntra)
            print "((",float(nbNodesTarget),"+",float(nbNodesGroup),",)/",float(n),",)*(",float(distGCglobToNewTargetCom),"**2.0)"
            print distGCglobToNewTargetCom,self.centreDeGraviteTotal,newGCcomCible
            exit()
            
            
        # faux 
        #newVarianceIntra = self.currentVarInter - (nbNodesSource/n)*varSource - (nbNodesTarget/n) * varSource + newComposanteComSource + ((nbNodesTarget+nbNodesGroup)/n)*newVarcomCible
        #newVarianceIntra=3.14159
        
        
        # renvoyer nouveau nb de noeuds...
        
        
        #classesNumber depend du nombre de groupes dans la communaute source
        
        #if(newVarcomCible[0] != 0):
        #    assert(newVarianceIntra<480000)
        
        
        
        assert(nbNodesSource>=nB)
        if  nbNodesSource == nB: # une classe en moins
            classesNumber = self.classesNumber-1
        else :
            classesNumber = self.classesNumber
        
        #!
        #assert(newVarianceIntra<=191.375)
        #print "<><>"
        #pprint(newGCcomSource)
        #if newGCcomSource==None:
        #    print "<><>"
        #    print (newVarianceIntra, -1, -1, newGCcomCible, newVarcomCible , classesNumber) 
        #    return (newVarianceIntra, -1, -1, newGCcomCible, newVarcomCible , classesNumber) # y aura des choses a rajouter
        #print "!!!!!!"
        #print "retour ",newVarianceIntra, newGCcomSource, newVarcomSource, newGCcomCible, newVarcomCible , classesNumber
        #pprint(newVarianceIntra)
        
        #assert(newVarianceIntra>=0)
        
        #print "self.varianceUnionTwoSets",moyGroup, varGroup, moyComTarget, varComTarget, nbNodesGroup, nbNodesTarget
        #print self.varianceUnionTwoSets(moyGroup, varGroup, moyComTarget, varComTarget, nbNodesGroup, nbNodesTarget)
        
        #print "resultats de moveGroup"
        #print (newVarianceIntra, newGCcomSource, newVarcomSource, newGCcomCible, newVarcomCible , classesNumber)
        #for i in self.comVar:
        #    print i, self.comVar[i], self.comGC[i], self.nbOfNodesInComs[i]
        
        #pprint(self.comVar[comSource])
        #pprint(self.comVar[comTarget])
        #print comSource, comTarget
        #pprint(newVarianceIntra)
        #if (self.comGC[comSource] != 0):
            #raw_input("resultats de moveGroup")
        #if (self.comGC[comTarget] != 0):
            #raw_input("resultats de moveGroup")
        
        # plutot inter pour la premiere
        
        
        
        
        
        valeurs=np.zeros(len(self.firstGlobalTfIdfTab))
        
        for i in range(len(self.firstGlobalTfIdfTab)):
            valeurs[i] = self.firstGlobalTfIdfTab[i][0]
        
        #pprint(valeurs)
        #exit()
        
        #pprint(self.currentPartition.node2com)
        
        #pprint(infererLaPartitionDeTousLesSommets(currentPartition=currentPartition))
        dic = self.hierarchy.infererLaPartitionDeTousLesSommets(currentPartition=self.currentPartition.node2com)
        #pprint(dic)
        #dic = self.currentPartition.node2com.copy()
        
        #pprint(dic)
        
        #assert(dic[group]==comSource)
        #, comSource, comTarget
        dic[group]=comTarget
        #pprint(dic)
        popVar=np.var(valeurs)
        moyPop=np.mean(valeurs)
        
        
        #exit()
        modAttr = -2 #ma.modulAttr(dic, valeurs, popVar, moyPop)
		
        
        return (modAttr, newVarianceIntra, newGCcomSource,newGCcomCible, classesNumber)
        #return (newVarianceIntra, newGCcomSource, newVarcomSource, newGCcomCible, newVarcomCible , classesNumber) # y aura des choses a rajouter
        #self.verifCoherence()
        
        
        
    def moveGravityCentersVersionGroups(self,newGroup, effectifGroup, initialCom2GravityCenter, initialAmountOfNodesInCom2, initialCom1GravityCenter,  initialAmountOfNodesInCom1):
        """
        On donne le nouveau group, le GC1 de destination et le poids de la com 1 et on retourne le nouveau centre de gravite
        Idem avec la communaute source
        """
        #pprint(initialCom1GravityCenter)
        
        newGC1=zeros([len(initialCom1GravityCenter)])
        newGC2=zeros([len(initialCom1GravityCenter)])
        
        #pprint(newGroup)
        
        #pprint(initialCom1GravityCenter)
        #pprint(self.groupGC[newGroup])
        #print "parameters",newGroup, effectifGroup, initialCom2GravityCenter, initialAmountOfNodesInCom2, initialCom1GravityCenter,  initialAmountOfNodesInCom1
        
        for i in range(len(initialCom1GravityCenter)):
            #print "initialAmountOfNodesInCom1",initialAmountOfNodesInCom1
            #print "initialCom1GravityCenter",initialCom1GravityCenter
            #pprint(newGroup)
            #pprint(self.groupGC[newGroup][i]) 
            #print "initialAmountOfNodesInCom1",initialAmountOfNodesInCom1
            
            
            #! corriger
            
            #print initialCom1GravityCenter[i]
            #pprint(self.groupGC)
            #print self.groupGC[newGroup]
            #[i]
            
            
            newGC1[i]=(initialAmountOfNodesInCom1*initialCom1GravityCenter[i]+self.groupGC[newGroup][i]*effectifGroup)/(initialAmountOfNodesInCom1+effectifGroup)
            assert(newGC1[i]!=float("inf"))
            assert(newGC1[i]!=float("-inf"))
            
            if((initialAmountOfNodesInCom2-effectifGroup)!=0):
                newGC2[i]=((initialAmountOfNodesInCom2)*initialCom2GravityCenter[i]-self.groupGC[newGroup][i]*effectifGroup)/(initialAmountOfNodesInCom2-effectifGroup)
                assert(newGC2[i]!=float("inf"))
                assert(newGC2[i]!=float("-inf"))
            else:#si la com2 n'existe plus:
                newGC2=None
        return (newGC1, newGC2)
    
    
    
    
    
    
        
    
    def varianceUnionTwoSets(self, moyA, varA, moyB, varB, nA, nB):
        '''
        http://math.stackexchange.com/a/80373/27820
        Retourne la valeur de la variance quand on unit 2 ensembles.
        moyA moyenne des element de A
        varA variance des element de A
        nA nombre d'element dans l'ensemble A
        
        Fonction verifiee
        '''
        out = zeros([len(moyA)])
        for i in range(len(moyA)):
            #pprint(varA)
            #pprint(varB)
            #print varA[i]
            #print varB[i]
            #print (nA*varA[i]+nB*varB[i])/(nA+nB)
            #print ((nA*nB)/((nA+nB)**2))*((moyA[i]-moyB[i])**2)
            #print ((float(nA)*float(nB))/((float(nA)+float(nB))**2))
            #print ((moyA[i]-moyB[i])**2)
            out[i]=(float(nA)*varA[i]+float(nB)*varB[i])/(float(nA)+float(nB)) + ((float(nA)*float(nB))/((float(nA)+float(nB))**2))*((moyA[i]-moyB[i])**2)
        return out #(nA*varA+nB*varB)/(nA+nB) + ((nA*nB)/((nA+nB)**2))*((moyA-moyB)**2)
        
        
    def __varianceDifferenceTwoSets(self, moyA, varA, moyB, varB, nA, nB):
        """
        http://math.stackexchange.com/a/80373/27820
        Retourne la valeur de la variance quand on enleve un ensemble B d'un autre plus grand A.
        """
        if nA==nB: # Plus de noeuds
            return None
        #print moyA, varA, moyB, varB, nA, nB
        #pprint((nA*varA+nB*varB)/(nA+nB) + ((nA*nB)/((nA+nB)**2))*((moyA-moyB)**2))
        
        return (nA*varA+nB*varB)/(nA+nB) + ((nA*nB)/((nA+nB)**2))*((moyA-moyB)**2)
        
        #var(a U b) = (nA*varA+nB*varB)/(nA+nB) + ((nA*nB)/((nA+nB)**2))*((moyA-moyB)**2)
        
        #formule a verifier:
        #var (a prive de  b) = ( (nA*varA)/(nA+nB)  - varAUnionB       + ((nA*nB)/((nA+nB)**2))*((moyA-moyB)**2)                      ) * ((nA+nB)/(-nB))
        
        
    def moveNode(self, node, comSource, comTarget   ): #, varInterSubCountComSource,   varInterSubCountComTarget    ):
        """
        On imagine que l'on deplace un noeud.
        On renvoie la varInter qu'il y aurait dans cette configuration, le cdg source et le cdg cible, et le nombre de classes
        Deprecated
        """
        exit()
        #on enleve le carre de la distance avec les anciens centres de gravite qu'on soustrait
        if (comSource == comTarget ):
            return self.currentVarInter,self.varInterSubGravityCenter[comSource] ,self.varInterSubGravityCenter[comTarget] ,self.classesNumber
        #pprint(self.currentPartition.node2com)    
        assert(self.currentPartition.node2com[node]==comSource)
        
        
        # on refait un calcul de style centre de gravite:
        # on prend la variance inter courante (une moyenne)
        res=self.currentVarInter
        
        # on enleve la variance inter d'avant de la communaute source * truc/n
        res=res-  self.__partialVarianceInterOfClass(comSource)* (len(self.currentPartition.dicInverse()[comSource]))   /len(self.firstGlobalTfIdfTab)   
        
        
        # on enleve la variance inter d'avant de la communaute cible * truc2/n
        res=res-  self.__partialVarianceInterOfClass(comTarget)* (len(self.currentPartition.dicInverse()[comTarget]))   /len(self.firstGlobalTfIdfTab)   
        
        
        # calcul des 2 nouveaux gc
        gc1dest, gc2source = self.moveGravityCenters( self.varInterSubGravityCenter[comTarget]      ,      self.firstGlobalTfIdfTab[self.globAuthorIndex[node]]    , (len(self.currentPartition.dicInverse()[comTarget]))  , self.varInterSubGravityCenter[comSource], (len(self.currentPartition.dicInverse()[comSource])))
        
        if (gc2source!=None): 
            classesNumber = self.classesNumber
            
        else: #si une communaute a disparu
            classesNumber = self.classesNumber -1
        
        
        # on ajoute la variance inter de la communaute source * truc-1/n
        if (gc2source!=None): # si gc2source n'a pas disparu
            res=res+     self.calculatePartOfVarInterWithGivenGC(gc2source)  * (len(self.currentPartition.dicInverse()[comSource])-1) / float(len(self.globalTfIdfTab))
        
        # on ajoute la variance inter de la communaute cible * truc2+1/n
        res=res+     self.calculatePartOfVarInterWithGivenGC(gc1dest)  *   (len(self.currentPartition.dicInverse()[comTarget])+1) / float(len(self.globalTfIdfTab))
        
        
        
        
        
        
        
        #on prend les nouveau centre et on calcule les nouveaux carres qu'on ajoute            
        
        assert(self.currentPartition.node2com[node]==comSource)
        
        return res, gc2source, gc1dest, classesNumber
    
    
    
    def __partialVarianceInterOfClass(self, cla): # a revoir
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
        print "cla"
        pprint(cla)
        
        dis=euclidean(self.varInterSubGravityCenter[cla],self.centreDeGraviteTotal)
        #print "distance "+ str(dis)
        #print "coef: "+str(len(self.currentPartition.dicInverse()[cla]))+"/"+str(len(self.globalTfIdfTab))+"  -----   "+str((len(self.currentPartition.dicInverse()[cla])/float(len(self.globalTfIdfTab))))
        #print (len(self.currentPartition.dicInverse()[cla])/float(len(self.globalTfIdfTab)))*(dis)**2
        return dis**2
        #dic=currentPartition.dicInverse()
        #listIndices=dic[cla]
    
    
    
    
    
    
    
    def reallyMoveGroup(self, group, comSource, comTarget, damax=-1,cnewVarianceIntra=None,cx1=None, cx2=None, cclassesNumber=None):
        '''
        appel de movegroup et maj de:
        CG pour coms
        var pour coms
        nb nodes dans les coms
        nombre de classes
        
        evidemment la varinter actuelle
        '''
        self.verifCoherence()
        
        if comSource==comTarget:
            return self.currentModAttr, self.currentVarInter
        
        
        #if comSource != comTarget:
        #    print group, comSource, comTarget
        #    raw_input()
            
        # y aura des choses a rajouter
        
        
        #if(cnewVarianceIntra!=None):
        
        #    newVarianceIntra, newGCcomSource, newVarcomSource, newGCcomCible, newVarcomCible, classesNumber = cnewVarianceIntra, cx1, cx2, cx3, cx4, cclassesNumber
        #else:
        #print newVarianceIntra
        
        #pprint(self.moveGroup( group, comSource, comTarget, damax=damax))
        (modAttr, newVarianceIntra, newGCcomSource,  newGCcomCible,  classesNumber) = self.moveGroup( group, comSource, comTarget, damax=damax)
        #if classesNumber==2000:
        #    exit()
        #print newVarianceIntra, cnewVarianceIntra
        
        #print newVarianceIntra
        
        #!
        #if(self.comGC[comSource][0]   != self.comGC[comTarget][0]):
        
        #    print newVarianceIntra, newGCcomSource, newVarcomSource, newGCcomCible, newVarcomCible , classesNumber
            #assert(newVarcomCible[0]!=0)
            #raw_input("important")
        #    pass
        #if(newVarcomCible[0] != 0):
        #    pprint(newVarcomCible)
        #    print "Bug"
        #    exit()
        
        
        
        self.currentVarInter=newVarianceIntra
        self.currentModAttr=modAttr
        
        
        
        
        self.comGC[comSource]=newGCcomSource
        self.comGC[comTarget]=newGCcomCible
        #self.comVar[comSource]=newVarcomSource
        #self.comVar[comTarget]=newVarcomCible
        
        
        #print "WWWvar"
        #pprint(self.comVar)
        #print "MMMvar"
        #print "WWWGC"
        #pprint(self.comGC)
        #print "MMMGC"
        #print "WWWnbnodes"
        #pprint(self.nbOfNodesInComs)
        #print "MMMnbnodes"
        
        
        self.verifCoherence()
        
        
        nbNodesGroup = self.hierarchy.nbNodesInAGroup(group)
        
        nbNodesSource= self.nbInitialNodes(comSource)
        #print "ajout: ",nbNodesGroup, nbNodesSource
        
        if nbNodesGroup== nbNodesSource:
            # supprimer la com source
            self.nbOfNodesInComs[comSource] = 0
            #del self.comVar[comSource]
            #pprint(self.currentPartition.node2com)
            #print self.classesNumber
            #raw_input("classesNumber diminue")
        else:
            self.nbOfNodesInComs[comSource] = self.nbOfNodesInComs[comSource]-nbNodesGroup
        self.nbOfNodesInComs[comTarget] = self.nbOfNodesInComs[comTarget] + nbNodesGroup 
        
        self.verifCoherence()
        
        
        
        #maj de la partition
        self.currentPartition.update(group,comTarget)
        #cp=deepcopy(self.currentPartition)
        #cp.node2com[group]=comTarget
        #self.currentPartition=APartition(cp.node2com)
        self.classesNumber=classesNumber
        self.verifCoherence()
        
        
        #nbElements=[]
        
        #centre de gravite et variance de chaque groupe
        #groupVariance=[]
        
        
        
        #centre de gravite et variance de chaque communaute
        # plutot cdg
        varInterSubGravityCenter=[]
        
        centreDeGraviteTotal=[]
        
        #la partition courante
        currentNode2Com={}
        
        euclideanDistanceOfComsToGlobGC=[]
        
        
        #nbInitialNodesInGroups={}
        
        # preciser ce truc
        #nbOfNodesInComs=[]
        
        
        #self.nbOfNodesInComs[comSource] = self.nbInitialNodes(comSource)-nbNodesGroup
        #self.nbOfNodesInComs[comTarget] = self.nbInitialNodes(comTarget)+nbNodesGroup
        #nbNodesSource= self.nbInitialNodes(comSource) 
        
        
        self.verifCoherence()
        
        
        #################################MODATTR
        # valeurs=np.zeros(len(self.firstGlobalTfIdfTab))
        
        # for i in range(len(self.firstGlobalTfIdfTab)):
            # valeurs[i] = self.firstGlobalTfIdfTab[i][0]
        
        # pprint(valeurs)
        # #exit()
        
        # dic = self.currentPartition.node2com.copy()
        
        # #pprint(dic)
        
        # #assert(dic[group]==comSource)
        # #, comSource, comTarget
        # #dic[group]=comTarget
        # #pprint(dic)
        # popVar=np.var(valeurs)
        # moyPop=np.mean(valeurs)
        
        
        # #exit()
        # modAttr = ma.modulAttr(dic, valeurs, popVar, moyPop)
        
        
        return modAttr, newVarianceIntra
        
        
        
        
        
    def inducedVarStatus(self,partition):
        '''
        Passage au meta, les com deviennent des groupes.
        '''
        
        #pprint(partition)
        
        #print(str(self))
        
        #print "avant induced"
        
        #pprint(self.currentPartition.node2com) 
        #pprint(self.nbInitialNodesInGroups)
        
        
        #pprint(self.comVar)
        #raw_input("etat de la variance")
        #pprint()
        
        
        self.verifCoherence()
        
        
        # a partir de la partition actuelle
        # on en fait une nouvelle, vierge
        
        
        #pprint(self.currentPartition.node2com)
        
        # on renumerote
        from communityThreshold import renumber
        renumberedPartition, nbClasses = renumber(self.currentPartition.node2com)
        
        renumberedPartition = partition
        
        #print "renumberedPartition variance"
        #pprint(renumberedPartition)
        
        effectiveInCom = zeros(nbClasses, dtype=int)
        
        #print "ancien compte de noeuds"
        #pprint(self.nbOfNodesInComs)
        
        #print "effectif des groups"
        #pprint(self.nbInitialNodesInGroups)
        
        
        
        
        # on compte le nouveau nombre de noeud dans chaque nouvelle communaute
        #for v in range(len(self.nbInitialNodesInGroups)):
        for v in range(len(self.nbInitialNodesInGroups)):
            #print v
            k=self.nbInitialNodesInGroups[v]
            #print k
        #for v,k in self.nbInitialNodesInGroups.items(): # pour chaque group
            effectiveInCom[renumberedPartition[v]] =  effectiveInCom[renumberedPartition[v]] + self.nbInitialNodesInGroups[v]
            #self.nbOfNodesInComs[k] = self.nbOfNodesInComs[k] + self.nbOfNodesInGroup[v]
        #pprint(effectiveInCom)
        #exit()
        # nombres de noeuds dans les nouveaux groupes
        self.nbInitialNodesInGroups = effectiveInCom
        
        #pprint(self.nbInitialNodesInGroups)
        #raw_input("nbInitialNodesInGroups")
        
        # nombres de groups dans les nouvelles communautes
        
        self.nbOfNodesInComs = zeros([nbClasses])
        
        

        
        
        #pprint(self.nbOfNodesInComs)
        
        #nouvelle node2com
        discPart = {} #zeros([nbClasses])
        j = 0
        for i in range(len(self.nbOfNodesInComs)):
            discPart[i] = i
            
        # la nouvelle partition (discrete)
        #pprint(discPart)
        
        

        
        
        
        # reaffectation des GC
        newComGC={}
        ap = APartition(renumberedPartition)
        inver=deepcopy(ap.dicInverse())
        
        

        
        for i in range(nbClasses):
            #!
            #print "i",i
            #print "discPart",discPart
            #print "node2com", self.currentPartition.node2com
            #print "com2node", self.currentPartition.com2node
            #print "renumberedPartition",renumberedPartition
            #pprint(self.comGC)
            #pprint(inver)
            #print "comGC"
            #pprint(self.comGC)
            #print "cg de i? i=",i
            #print "i contient",inver[i]
            #print "inver[i].pop()",inver[i].pop()
            #print "inver[i]",inver[i]
            #pprint(("node2com[inver[i]]",self.currentPartition.node2com[inver[i].pop()]))
            #print "node2com[inver[i].pop()]",self.currentPartition.node2com[inver[i].pop()]
            
            
            #pprint(self.comGC)
            #pprint(self.currentPartition.node2com)
            #pprint(inver)
            
            elt=inver[i].pop()
            inver[i].add(elt)
            
            
            
            newComGC[i] = self.comGC[self.currentPartition.node2com[elt]]
            #pprint(newComGC)
            #self.newComGC[i]=self.comGC[inver[self.currentPartition.com2node[i]]]
            #pprint(inver)
            #raw_input("Press ENTER to exit")
            
            
        #pprint(newComGC)
        
        self.comGC=newComGC
        
        # pour eviter les effets de bord parce qu'on a pope
        inver=deepcopy(ap.dicInverse())
        #pprint(inver)
        #pprint(self.comGC)
        #exit()
        
        
        
        
        self.groupGC={}
        
        # reaffectation des Variances
        #!
        newComVar = {}
        #pprint(self.comVar)
        for i in range(nbClasses):
        #var actuelle
            #print "i",i
            #newComVar[i] = self.comVar[self.currentPartition.node2com[inver[i].pop()]]
            
            self.groupGC[i] = self.comGC[i]
            #self.groupMoyenne[i] = self.comGC[i]
            #self.groupVariance[i] = 0
            
        #pprint(newComVar)
        
        #raw_input("pprint(newComVar)")
        
        #self.groupVariance = deepcopy(newComVar)
        #self.comVar = newComVar
        
        
        
        
        
        
        self.currentPartition = APartition(discPart)
        #pprint(self.currentPartition.node2com)
        
        
        # on ajoute la partition courante dans la hierarchie
        
        
        #print(renumberedPartition)
        
        self.hierarchy.ajouterUnePartition(renumberedPartition)
        
        
        
        
        
        euclideanDistanceOfComsToGlobGC=[]
        
        
        
        
        
        #print "apres induced"
        
        #print "self.groupVariance"
        #pprint(self.groupVariance)
        #pprint(self.currentPartition.node2com) 
        #pprint(self.hierarchy.numberOfNodes)
        
        
        #pprint(self.nbInitialNodesInGroups)
        #raw_input("nbInitialNodesInGroups")
        dpcp = deepcopy(self.nbInitialNodesInGroups)
        self.nbOfNodesInComs = dpcp
        #raw_input("finInducedVarStatus")
        self.verifCoherence()
        #print(str(self))
        
        #print str(self)
        #exit()
        
        
        
    def modulariteDesAttributs():
        
        # (la variance Inter de P) moins (la variance Inter attendue si la partition avait ete la meme (en particulier autant d'individus dans chaque classe) mais que les individus avaient ete repartis aleatoirement dans chaque classe)
        
        
        
        
        # calculer la variance intra de si les sommets avaient ete repartis aleatoirement:
        
        
        
        nbOfNodesInComs
        
        
        
        
        varInter - nbOfNodesInComs
        
        
        return -1
        
        
        
        
        
    def __str__(self):
        out = "-----VarStatus-----\n|hierarchie\n|"
        out=out+ str(self.hierarchy)
        out=out+ "\n|currentPartition\n|"
        out=out+str(self.currentPartition)
        out=out+ "\n|comGC\n|"
        out=out+str(self.comGC)
        out=out+ "\n|comVar\n|"
        out=out+str(self.comVar)
        
        
        
        out=out+ "\n|groupVariance\n|"
        out=out+str(self.groupVariance)
        
        
        #out=out+ "\ngroupMoyenne\n"
        #out=out+str(self.groupMoyenne)
        
        
        out=out+ "\n|groupGC\n|"
        out=out+str(self.groupGC)
        
        
        
        out=out+ "\n-----FinVarStatus-----\n"
        return out
# Les fonctions a produire sont
# - variInter( node, dest)
# - (variTotale)
# - hypothese de bougeage: moveNode(self, node, comSource, comTarget)
# - bougeage effectif: reallyMoveNode(self, node, comSource, comTarget)



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
    
    def __init__(self, graphDeprecated, globalTfIdfTab, partition, globAuthorIndex ) :
        self.globAuthorIndex = globAuthorIndex
        self.varInterSubCount = dict([])
        self.globalTfIdfTab=globalTfIdfTab
        self.centreDeGraviteTotal=self.moyenne_points(range(len(globalTfIdfTab)))
        
        #print "init: cdg total"
        #pprint(self.centreDeGraviteTotal)
        
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
        Deprecated
        """
        exit()
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
        
        self.currentPartition.update(node,comTarget)
        
        assert(self.currentPartition.node2com[node] == comTarget)
        
        
        
        #modAttr=
        
        return modAttr, self.currentVarInter

    




    def moveGravityCenters(self,initialCom1GravityCenter, newPoint, initialAmountOfNodesInCom1, initialCom2GravityCenter, initialAmountOfNodesInCom2):
        """
        On donne le nouveau point, le GC1 de destination et le poids de la com 1 et on retourne le nouveau centre de gravite
        Idem avec la communaute source
        """
        newGC1=zeros([len(initialCom1GravityCenter)])
        newGC2=zeros([len(initialCom1GravityCenter)])
        
        for i in range(len(initialCom1GravityCenter)):
            newGC1[i]=(initialAmountOfNodesInCom1*initialCom1GravityCenter[i]+newPoint[i])/(initialAmountOfNodesInCom1+1)
            #assert(newGC1[i]!=float("inf"))
            #assert(newGC1[i]!=float("-inf"))
            #si la com2 n'existe plus:
            if((initialAmountOfNodesInCom2-1)!=0):
                newGC2[i]=((initialAmountOfNodesInCom2)*initialCom2GravityCenter[i]-newPoint[i])/(initialAmountOfNodesInCom2-1)
                #assert(newGC2[i]!=float("inf"))
                #assert(newGC2[i]!=float("-inf"))
            else:
                newGC2=None
        return (newGC1, newGC2)
    
    
    def metaMoveGravityCenters(self,initialCom1GravityCenter, newPoint, initialAmountOfNodesInCom1, initialCom2GravityCenter, initialAmountOfNodesInCom2, nbNodesInGroup):
        """
        
        """
        newGC1=zeros([len(initialCom1GravityCenter)])
        newGC2=zeros([len(initialCom1GravityCenter)])
        for i in range(len(initialCom1GravityCenter)):
            newGC1[i]=(initialAmountOfNodesInCom1*initialCom1GravityCenter[i]+newPoint[i])/(initialAmountOfNodesInCom1+nbNodesInGroup)
            assert(newGC1[i]!=float("inf"))
            assert(newGC1[i]!=float("-inf"))
            if((initialAmountOfNodesInCom2-nbNodesInGroup)!=0):
                newGC2[i]=((initialAmountOfNodesInCom2)*initialCom2GravityCenter[i]-newPoint[i])/(initialAmountOfNodesInCom2-nbNodesInGroup)
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



def testMetaMoveGravityCenters():
    pass

def testMoveNode():
    part=APartition({"alpha":0,"beta":0,"delta":0,"gamma":1})
    global globTfIdfTab
    G=graphExample()
    vs = VarianceStatus(G, globTfIdfTab, part, globAuthorIndex)
    
    print(vs.__str__())
    
    #print "move node"
    var, rien1, rien2, classesNumber = vs.moveNode( "gamma", 1, 0)#, varInterSubCountComSource,   varInterSubCountComTarget    )
    #print "var inter avec noeud deplace: "+str(var)
    # la var inter de la partition unitaire est zero
    assert(var<0.01)
    assert(var>-0.01)
    print "Test de bougeage de noeud passe avec succes\n--------------------Fin du test--------------------"


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
    vs = VarianceStatus(G, globTfIdfTab, part, globAuthorIndex)
    
    
    
    
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
    vs = VarianceStatus(G, globTfIdfTab, part, globAuthorIndex)
    
    print str(vs)
    
    
    
    
    
    
    
    
    exit()
    
    global globTfIdfTab
    G=graphExample()
    print G.nodes()
        
    print globTfIdfTab
    
    
def looseAssertlt(v1,v2):    
    #pprint(v1)
    #pprint(v2)
    assert(math.floor(100000000.0*v1)/100000000.0<=v2)
    
    
    
    
def test():
    #assert(varianceUnionTwoSets(moyA, varA, moyB, varB, nA, nB) ==   )
    assert(varianceUnionTwoSets(1, 0, 2, 0, 1, 1) == 0.25  )

if __name__ == "__main__" :
    __main()




