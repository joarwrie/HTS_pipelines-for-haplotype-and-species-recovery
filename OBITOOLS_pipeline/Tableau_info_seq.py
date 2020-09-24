#!/usr/bin/python2
# coding: utf-8

import re
import sys

# Ouverture du fichier resultat
fichierFinal=open(sys.argv[2], "w")

# Pour chaque ligne du fichier
with open(sys.argv[1]) as f:
	for ligne in f:
		if re.search("^>", ligne) is not None:
			seqid=(ligne.split(" ")[0]).split(">")[1]
			liste_params=((ligne.split(" merged_Echantillon={")[1]).split("};")[0]).split(", ")
			for params in liste_params:
				fichierFinal.write (seqid + "\t" + params.split("'")[1] + "\t" + params.split(": ")[1] + "\n")
			# ---- Done ----
		# ---- Fi ----		
	# ---- Done ----		
# ---- Fin with ----

fichierFinal.close()

