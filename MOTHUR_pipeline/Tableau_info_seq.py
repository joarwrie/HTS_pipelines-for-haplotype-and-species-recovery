#!/usr/bin/python2
# coding: utf-8

import re
import sys

# Ouverture du fichier resultat
fichierFinal=open(sys.argv[2], "w")

# Pour chaque ligne du fichier
with open(sys.argv[1]) as f:
	for ligne in f:
		seqid=ligne.split("\t")[0]
		liste_seq=(ligne.split("\t")[1]).split(",")
		for sequence in liste_seq:
			fichierFinal.write (seqid + "\t" + sequence + "\n")
		# ---- Done ----
	# ---- Done ----		
# ---- Fin with ----

fichierFinal.close()
