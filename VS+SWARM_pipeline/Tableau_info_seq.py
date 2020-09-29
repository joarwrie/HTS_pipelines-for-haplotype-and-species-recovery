#!/usr/bin/python2
# coding: utf-8

import re
import sys

# Ouverture du fichier resultat
fichierFinal=open(sys.argv[2], "w")

# Pour chaque ligne du fichier
with open(sys.argv[1]) as f:
	for ligne in f:
		seqid=(ligne.split(" ")[0]).split(";")[0]
		liste_seq=ligne.split(" ")
		for sequence in liste_seq:
			fichierFinal.write (seqid + "\t" + sequence.split(";")[0] + "\n")
		# ---- Done ----
	# ---- Done ----		
# ---- Fin with ----

fichierFinal.close()