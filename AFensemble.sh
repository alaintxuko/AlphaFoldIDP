#!/bin/bash

#Define the number of structures to generate for the ensemble
Nens=100

###############################################################################

#uncompress the zip files
for i in $(ls $NAME*.zip); do unzip -o $i; done

#extract the root name of the files
RNAME=$(ls -1d *pdb | awk -F "_" '{print $1}' | head -n1)
NAME=$(ls -1d *pdb | awk -F "_" '{print $1}' | head -n1)_f

#set variables for the residue numbers
INI=$(ls -1d *_relaxed_rank_1_model_*.pdb | awk -F "$NAME" '{print $2}' | awk -F "_" '{print $1}' | sort -h | head -n1)
LAST=$(ls -1d *_relaxed_rank_1_model*.pdb | awk -F "$NAME" '{print $2}' | awk -F "_" '{print $1}' | sort -h | tail -n1)
INI1=$(echo $INI | awk -F "-" '{print $1}')
INI2=$(echo $INI | awk -F "-" '{print $2}')
LAST1=$(echo $LAST | awk -F "-" '{print $1}')
LAST2=$(echo $LAST | awk -F "-" '{print $2}')
Ndiff=$[$INI1-1]
N=$[$INI2-$INI1+1]
N1=$[$N-1]
N2=$[$N-2]
Npept=$[$LAST1-$INI1+1]
Size=$[$LAST2-$INI1+1]
Nfiles=$(ls -1d *_relaxed_rank_1_model_*.pdb | grep -c pdb)

#loop to simplify the filenames
n=0
for i in `seq $INI1 $LAST1`; do
n=$[$n+1]
#    cp "$NAME$i"-*_relaxed_rank_1_model_*.pdb "$NAME$n"_relaxed_rank_0.pdb
#    cp "$NAME$i"-*_relaxed_rank_2_model_*.pdb "$NAME$n"_relaxed_rank_1.pdb
#    cp "$NAME$i"-*_relaxed_rank_3_model_*.pdb "$NAME$n"_relaxed_rank_2.pdb
#    cp "$NAME$i"-*_relaxed_rank_4_model_*.pdb "$NAME$n"_relaxed_rank_3.pdb
#    cp "$NAME$i"-*_relaxed_rank_5_model_*.pdb "$NAME$n"_relaxed_rank_4.pdb
    mv "$NAME$i"-*_relaxed_rank_1_model_*.pdb "$NAME$n"_relaxed_rank_0.pdb
    mv "$NAME$i"-*_relaxed_rank_2_model_*.pdb "$NAME$n"_relaxed_rank_1.pdb
    mv "$NAME$i"-*_relaxed_rank_3_model_*.pdb "$NAME$n"_relaxed_rank_2.pdb
    mv "$NAME$i"-*_relaxed_rank_4_model_*.pdb "$NAME$n"_relaxed_rank_3.pdb
    mv "$NAME$i"-*_relaxed_rank_5_model_*.pdb "$NAME$n"_relaxed_rank_4.pdb
done

#check that the dataset is complete
if [ $Npept != $Nfiles ]; then
    echo The number of peptides is not correct
    exit 0
fi

#extract the sequence from the files
SEQ=""
for i in `seq $[$INI1-$Ndiff] $N $[$LAST1-$Ndiff]`; do
    SEQ2=$(cat "$NAME$i"_relaxed_rank_0.pdb | awk '/ATOM/ && $3 == "CA" && $5 == "A" {print $4}' | tr '\n' ' ' | sed 's/ALA/A/g;s/CYS/C/g;s/ASP/D/g;s/GLU/E/g;s/PHE/F/g;s/GLY/G/g;s/HIS/H/g;s/ILE/I/g;s/LYS/K/g;s/LEU/L/g;s/MET/M/g;s/ASN/N/g;s/PRO/P/g;s/GLN/Q/g;s/ARG/R/g;s/SER/S/g;s/THR/T/g;s/VAL/V/g;s/TRP/W/g;s/TYR/Y/g' | sed 's/ //g')
    SEQ=$SEQ$SEQ2
done
if (( $Size % $N != 0 )); then
end=$[$LAST1-$i-$Ndiff]
SEQ2=$(cat "$NAME$[$LAST1-$Ndiff]"*_relaxed_rank_0.pdb | awk '/ATOM/ && $3 == "CA" && $5 == "A" {print $4}' | tr '\n' ' ' | sed 's/ALA/A/g;s/CYS/C/g;s/ASP/D/g;s/GLU/E/g;s/PHE/F/g;s/GLY/G/g;s/HIS/H/g;s/ILE/I/g;s/LYS/K/g;s/LEU/L/g;s/MET/M/g;s/ASN/N/g;s/PRO/P/g;s/GLN/Q/g;s/ARG/R/g;s/SER/S/g;s/THR/T/g;s/VAL/V/g;s/TRP/W/g;s/TYR/Y/g' | sed 's/ //g')
SEQ=$SEQ${SEQ2: -$end}
fi
ResN=$(echo $SEQ | awk '{ print length }')
FragN=$[$ResN-$N+1]

#clean a previous run
rm scr_pp.pml
rm -r 000ensemble*
rm -r 00000END

#create output folder
mkdir 00000END

ensF=0 #final structures
ens=0 #counter of analyzed structures
e=0 #counter of rank combinations

#loop to generate the ensemble
while [ $ensF -le $[$Nens-1] ]; do
    e=$[$e+1]
    #extraction of the rank combination from the random list
    lineA=$(cat random.list | sed -n ''$e'p')
    array[0]=$(echo $lineA | awk '{print $1}')
    array[1]=$(echo $lineA | awk '{print $2}')
    array[2]=$(echo $lineA | awk '{print $3}')
    array[3]=$(echo $lineA | awk '{print $4}')
    array[4]=$(echo $lineA | awk '{print $5}')
    array[5]=$(echo $lineA | awk '{print $6}')
    #loop over the different starting points
    for i in `seq 1 $N2`; do
	ens=$[$ens+1]
	mkdir 000ensemble$ens
	c=0
	#add pdbs to load in pymol script
	for j in `seq $i $N2 $FragN`; do
	    k=${array[$c]}
	    echo "load "$NAME$j"_relaxed_rank_$k.pdb, p$j" >> scr_pp.pml
	    lastF=$j
	    if [ $c -ne 5 ]; then
		c=$[$c+1]
	    else
		c=0
	    fi
	done
	#add alignment commands to pymol script
	for j in `seq $[$i+$N2] $N2 $FragN`; do
	    echo "align (p$j and resi 1 and name ca+c) or (p$j and resi 2 and name n), (p$[$j-$N2] and resi $[1+$N2] and name ca+c) or (p$[$j-$N2] and resi $[2+$N2] and name n)" >> scr_pp.pml
	done
	lastR=$[$lastF+$N2]
	#add N-terminal residues if needed
	if [ $i -ne 1 ]; then
	    echo "load "$NAME"1_relaxed_rank_$k.pdb, p1" >> scr_pp.pml
	    echo "align (p1 and resi $i and name ca+c) or (p1 and resi $[$i+1] and name n), (p$i and resi 1 and name ca+c) or (p$i and resi 2 and name n)" >> scr_pp.pml
	    echo "save 000ensemble$ens/p1.pdb, p1 and resi 1-$i" >> scr_pp.pml
	fi
	#renumber residues
	for j in `seq $i $N2 $FragN`; do
	    j1=$[$j-1]
	    echo "alter p$j, resi=str(int(resi)+$j1)" >> scr_pp.pml
	done
	echo "sort" >> scr_pp.pml
	#save pdbs
	for j in `seq $i $N2 $FragN`; do
	    if  [ $j -ne 1 ]; then
	        if  [ $j -ne $FragN ]; then
		    echo "save 000ensemble$ens/p$j.pdb, p$j and resi $[$j+1]-$[$j+$N2]" >> scr_pp.pml
	        else
		    echo "save 000ensemble$ens/p$j.pdb, p$j and resi $[$j+1]-$[$j+$N2+1]" >> scr_pp.pml
	        fi
	    else
		echo "save 000ensemble$ens/p1.pdb, p1 and resi 1-$[1+$N2]" >> scr_pp.pml
	    fi
	done
	#add C-terminal residues if needed
	if  [ $lastF -ne $FragN ]; then
	    echo "load "$NAME$FragN"_relaxed_rank_$k.pdb, p$FragN" >> scr_pp.pml
	    j1=$[$FragN-1]
	    echo "alter p$FragN, resi=str(int(resi)+$j1)" >> scr_pp.pml
	    echo "sort" >> scr_pp.pml
	    echo "align (p$FragN and resi $lastR and name ca+c) or (p$FragN and resi $[$lastR+1] and name n), (p$lastF and resi $lastR and name ca+c) or (p$lastF and resi $[$lastR+1] and name n)" >> scr_pp.pml
	    echo "save 000ensemble$ens/p$[$lastR+1].pdb, p$FragN and resi $[$lastR+1]-$ResN" >> scr_pp.pml
	fi
	#run pymol script
	pymol -c scr_pp.pml
	#combine all fragment pdbs in one final pdb
	cat 000ensemble$ens/p*.pdb | grep -v END > 000ensemble$ens/000ensemble$ens.pdb
	echo "load 000ensemble$ens/000ensemble$ens.pdb" > scr_pp.pml
	echo "save 000ensemble$ens/000ensemble$ens.pdb" >> scr_pp.pml
	pymol -c scr_pp.pml
	rm scr_pp.pml
	#check for clashes and save final pdb
	echo "run get_raw_distances.py" > scr_checkclashes.pml
	echo "load 000ensemble$ens/000ensemble$ens.pdb" >> scr_checkclashes.pml
	echo "remove hydro" >> scr_checkclashes.pml
	echo "distance dist, 000ensemble$ens, 000ensemble$ens, 1.2, 0" >> scr_checkclashes.pml
	echo "get_raw_distances dist" >> scr_checkclashes.pml
	pymol -c scr_checkclashes.pml > distances.txt
	CLASHES=$(cat distances.txt | grep get_raw_distances | grep -c 000ensemble$ens)
	rm distances.txt
	rm scr_checkclashes.pml
	if [ $CLASHES -eq 0 ] ; then
	    ensF=$[$ensF+1]
	    echo "load 000ensemble$ens/000ensemble$ens.pdb" > scr_save.pml
	    echo "remove hydro" >> scr_save.pml
	    echo "save 000ensembles_calc/Str$ensF.pdb, 000ensemble$ens" >> scr_save.pml
	    echo "alter 000ensemble$ens, resi=str(int(resi)+$INI1-1)" >> scr_save.pml
	    echo "sort" >> scr_save.pml
	    echo "save 000ensemble$ens/Str$ensF.pdb, 000ensemble$ens" >> scr_save.pml
	    pymol -c scr_save.pml
	    cp 000ensemble$ens/Str$ensF.pdb 00000END/Str$ensF.pdb
	    rm scr_save.pml
	fi
	#finish the loop when the number of structures of the ensemble is reached
	if [ $ensF -eq $Nens ] ; then
	    break
	fi
    done #starting point loop
done #while loop
#save the sequence file
echo $SEQ > 00000END/$RNAME.seq
