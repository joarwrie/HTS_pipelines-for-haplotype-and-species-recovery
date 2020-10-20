#! /bin/bash

touch AllSamples_final_table_$3.shared

old_IFS=$IFS
IFS=$'\n'

compt=0
for ligne in $(cat $1);do
	(( compt++ ))
	echo $compt
	Char1=$(echo $ligne | head -c 1)
	if [ $Char1 == "#" ];then
		echo $ligne | sed 's/Representative_Sequence/Representative_Sequence\tOTU/g' >> AllSamples_final_table_$3.shared
	else
		Nom=$(echo $ligne | cut -f1)
		OTU=$(grep ${Nom} $2 | cut -f1)
		test $(echo $OTU | grep -c "M02439") -eq 0 && OTU="FAIL"
		echo $ligne | sed "s/${Nom}/${Nom}\t${OTU}/g" >> AllSamples_final_table_$3.shared
	fi
done

echo $compt

IFS=${old_IFS}