#!/bin/bash
#Define the path to SSP software (read/write access needed)
SSPDIR=/path/to/SSP/folder
#Define the number of structures to analyze
Nens=100
#Define SSP threshold
SSPmin=0.20

###############################################################################

echo Usage: ./analysisEnsemble

#Define initial variables
Pf=Str
PWD=$(pwd)
PWD2=$PWD
SEQ=$(cat "$Pf"1.pdb | awk '/ATOM/ && $3 == "CA" {print $4}' | tr '\n' ' ' | sed 's/ALA/A/g;s/CYS/C/g;s/ASP/D/g;s/GLU/E/g;s/PHE/F/g;s/GLY/G/g;s/HIS/H/g;s/ILE/I/g;s/LYS/K/g;s/LEU/L/g;s/MET/M/g;s/ASN/N/g;s/PRO/P/g;s/GLN/Q/g;s/ARG/R/g;s/SER/S/g;s/THR/T/g;s/VAL/V/g;s/TRP/W/g;s/TYR/Y/g' | sed 's/ //g')
sl=$(echo $SEQ | awk '{ print length }')
echo $SEQ > $Pf.seq

#clean a previous run
rm ssp/res*
rm SSP.out

#loop the 100 structures of the e
mkdir ssp
for i in `seq 1 $Nens`; do
    if test -f "$Pf$i.pdb"; then
	#Run SPARTA+
	sparta+ -in $Pf$i.pdb
	#Prepare files from SPARTA+ for SSP
	Fres=$[$[$(cat pred.tab | grep FIRST_RESID | awk '{print $3}')]-1]
	cat pred.tab | grep CA | grep -v WARNING | grep -v SEQUENCE | awk -v Fres=$Fres '{print $1-Fres" "$5}' > SPA$i.ca
	cat pred.tab | grep CB | grep -v WARNING | awk -v Fres=$Fres '{print $1-Fres" "$5}' > SPA$i.cb
	cp SPA$i.ca $SSPDIR/
	cp SPA$i.cb $SSPDIR/
	cp $Pf.seq $SSPDIR/
	#Run SSP
	cd $SSPDIR
	./ssp -s $Pf.seq -ca SPA$i.ca -cb SPA$i.cb > out
	cp $SSPDIR/out $PWD2/ssp/$i.out
	#clean temporary files
	rm SPA$i.ca
	rm SPA$i.cb
	rm $Pf.seq
	rm out
	cd $PWD2
	rm SPA$i.ca
	rm SPA$i.cb
	#Save the SSP data per residue for further analysis
	NLINES=$(awk '{nlines++} END {print nlines}' ssp/$i.out)
	for pp in `seq 1 $NLINES`; do
	    res=$(sed -n ''$pp'p' ssp/$i.out | awk '{print $1}')
	    sspV=$(sed -n ''$pp'p' ssp/$i.out | awk '{print $2}')
	    echo $sspV >> ssp/res$res
	done
    fi
done
#clean temporary files
rm $Pf.seq
rm pred.tab
rm struct.tab

#Calculate average and std values per residue
for i in `seq 1 $sl`; do
    echo $i $(cat ssp/res$i | awk -F ' '  '{   sum=sum+$1 ; sumX2+=(($1)^2)} END { avg=sum/NR; printf "%f %f\n", avg, sqrt(sumX2/(NR-1) - 2*avg*(sum/(NR-1)) + ((NR*(avg^2))/(NR-1)))}') >> SSP.out
done

#Extract the filtered data with the defined threshold (SSPmin)
#Define initial variables
SSPf=SSP.out

#clean previous runs
rm SS3finder.list

#loop to find groups of 3 residues with values over/below the thresholds
c=0
for i in `seq 1 $sl`; do
    sspL=$(sed -n ''$i'p' $SSPf)
    ssp=$(echo $sspL | awk '{print $2}')
    if [ "$ssp" == "" ]; then
	c=0
    else
	ssp2=${ssp#-}
	sspC=`echo "$ssp2  > $SSPmin" | bc`
	if [ $sspC == 1 ]; then
	    c=$[$c+1]
	    if [ $c == 1 ]; then
	        echo $sspL > SS3finder.temp
	    elif [ $c == 2 ]; then
	        echo $sspL >> SS3finder.temp
	    elif [ $c == 3 ]; then
	        echo $sspL >> SS3finder.temp
	        cat SS3finder.temp >> SS3finder.list
	    else
	        echo $sspL >> SS3finder.list
	    fi
	else
	    c=0
	fi
    fi
done
rm SS3finder.temp

#loop to fill empty spaces
NLINES=$(awk '{nlines++} END {print nlines}' SS3finder.list)
ini=1
for k in `seq 1 $sl`
do
    found=0
    for i in `seq $ini $NLINES`
    do
    line=$(sed -n ''$i'p' SS3finder.list)
    res=$(echo $line | awk '{print $1}' | grep -o '[0-9]\+')
    if [ $[$res] = $[$k] ]; then
	found=1
	ini=$i
        echo $line >> SSPfilt.out
        break
    fi
    if [ $[$res] -gt $[$k] ]; then
        break
    fi
    done
    if [ $[$found] = 0 ]; then
	echo $k >> SSPfilt.out
    fi
done
rm SS3finder.list
cat SSP.out | awk -v Fres=$Fres '{print $1+Fres" "$2" "$3}' > SSPtemp.out
mv SSPtemp.out SSP.out
cat SSPfilt.out | awk -v Fres=$Fres '{print $1+Fres" "$2" "$3}' > SSPfilttemp.out
mv SSPfilttemp.out SSPfilt.out
