totemClust
==========

ToTeM Complex networkx clustering


Appeler best_partition du fichier communityThreshold.py.

Des exemples d'appels et d'usage de networkx sont propos�s dans testsCommunityThreshold.py

best_partition(graph, 8, None, None, None)
Entr�es : 
- le graphe
- un seuil (plus utilis� je crois)
- un dictionnaire qui assigne un auteur � son "i" dans le tableau des attributs "tfIdfTab"
- une partition de d�part (pas de SAV sur ce point, � v�rifier)
- la partition de la v�rit� terrain, pour faire une �valuation des partitions interm�diaires en cours d'ex�cution de l'algorithme (v�rifier le comportement si on s'en sert �galement)
