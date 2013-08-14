totemClust
==========

ToTeM Complex networkx clustering


Appeler best_partition du fichier communityThreshold.py.

Des exemples d'appels et d'usage de networkx sont proposés dans testsCommunityThreshold.py

best_partition(graph, 8, None, None, None)
Entrées : 
- le graphe
- un seuil (plus utilisé je crois)
- un dictionnaire qui assigne un auteur à son "i" dans le tableau des attributs "tfIdfTab"
- une partition de départ (pas de SAV sur ce point, à vérifier)
- la partition de la vérité terrain, pour faire une évaluation des partitions intermédiaires en cours d'exécution de l'algorithme (vérifier le comportement si on s'en sert également)
