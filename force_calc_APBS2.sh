#!/bin/bash

if [ -z $1 ] ; then 
    echo "Usage: $0 <Calculation directory to analyze> " 
    exit 
    fi 

if [ ! -d $1 ] ; then 
    echo "ERROR: $1 directory does not exist. Have you run this MD yet?" 
    exit 
    fi 

printf "\n\t*** Program Beginning ***\n\n" 

MOLEC=$1
TOP=$(PWD)
XTC=$TOP/$MOLEC/nopbc/$MOLEC.nopbc.xtc 
TPR=$TOP/$MOLEC/Production/$MOLEC.tpr
tot=100000 #ps
numFrames=1000 

dirName="APBS_fixed"

pdb2pqr="/Users/jfirst/pdb2pqr-osx-bin64-2.1.0/pdb2pqr"

if [[ ${MOLEC:0:8} == 'SNase_WT' || ${MOLEC} == "SNase_V104X" ]] ; then 
    echo "No nitrile on this molecuele. Nothing to do." 
    exit 
    fi 

if (( $tot % 4 )) ; then 
    echo "WARNING: $tot % 4 != 0. There may not be frames for requested time steps" 
    fi 
timeStep=$(echo "$tot / $numFrames" | bc) #ps 
if [[ $timeStep -lt 20 ]] ; then 
    echo "ERROR: Requested time step is $timeStep. Frames only printed every 20 ps." 
    exit ; fi 
logFile=$TOP/$MOLEC/forces.log
errFile=$TOP/$MOLEC/forces.err

if [ ! -f $XTC ] ; then 
    echo "ERROR: $XTC not found. Did the calculation finish?" 
    exit 
    fi 

if [ ! -f $TPR ] ; then 
    echo "ERROR: $TPR not found " 
    exit 
    fi 

cd $MOLEC

check(){
    for arg in $@ ; do 
        if [ ! -s $arg ] ; then 
            echo ; echo "ERROR: $arg does not exist" 
            exit 
            fi 
        done 
}

if [ ! -d $dirName ] ; then mkdir $dirName ; fi 
cd $dirName
cp ../../free_energy_files/AMBER.DAT . 
cp ../../free_energy_files/AMBER.names . 
cp ../../GMXFF/*.dat . 
cp -r ../../GMXFF/amber03.ff . 

check AMBER.DAT AMBER.names

printf "\n\t*** Pre-Compress      ***\n\n" 
printf "\tExtracting relevant frames.............." 
rm compress.xtc index.ndx 
if [ ! -f compress.xtc ] ; then 
    gmx select -s $TPR -on index.ndx -select "group \"Protein\" or resname \"CA\"" >> $logFile 2>> $errFile 
    #gmx select -s $TPR -on index.ndx -select "group \"Protein\"" >> $logFile 2>> $errFile 
    check index.ndx

    gmx trjconv -f $XTC -s $TPR -o compress.xtc -n index.ndx -b 0 -e $tot -dt $timeStep -tu ps >> $logFile 2>> $errFile 
    check compress.xtc
    printf "Complete\n" 
else 
    printf "Skipped\n" 
    fi 

XTC=compress.xtc

printf "\n\t*** pdb2pqr and APBS    ***\n\n" 
if [ -f rxn_field.out ] ; then rm rxn_field.out ; fi 
if [ -f coloumb_field.out ] ; then rm coloumb_field.out ; fi 
for frame in $(seq 0 $timeStep $tot) ; do 
    printf "\tReading %5i of %5i..." $frame $tot

    if [ ! -f time_${frame}.pdb ] ; then 
        gmx trjconv -f $XTC -s $TPR -n index.ndx -o time_${frame}.pdb -dump $frame -tu ps >> $logFile 2>> $errFile 
        check time_${frame}.pdb 

        ####Replace C-Termin GLN with resname CGLN
        sed "s/OC1 GLN/O   GLN/" time_${frame}.pdb > temp.pdb 
        sed "s/OC2 GLN/OXT GLN/" temp.pdb > time_${frame}.pdb
        fi 



    if [ ! -f time_${frame}.pqr ] ; then 
        $pdb2pqr time_${frame}.pdb time_${frame}.pqr --userff AMBER.DAT --usernames AMBER.names --assign-only >> $logFile 2>> $errFile
        check time_${frame}.pqr 
        fi 

    if [ ! -f time_${frame}+DUM.pqr ] ; then 
        python ../../free_energy_files/add_dummy_pqr.py time_${frame}.pqr 0.1 CNC CD NE > time_${frame}+DUM.pqr 
        fi 

    if [[ ! -f time_${frame}_78.txt || ! -f time_${frame}_1.txt ]] ; then 
        sed "s/SDIE/78/" ../../free_energy_files/force_temp.in | sed "s/FRAME/${frame}/" >> time_${frame}_78.in 
        check time_${frame}_78.in 
        apbs time_${frame}_78.in >> $logFile 2>> $errFile 
        check time_${frame}_78.txt

        sed "s/SDIE/1/" ../../free_energy_files/force_temp.in | sed "s/FRAME/${frame}/" >> time_${frame}_1.in 
        check time_${frame}_1.in 

        apbs time_${frame}_1.in >> $logFile 2>> $errFile 
        check time_${frame}_1.txt
        fi 

    ../../free_energy_files/read_apbs_rxn_field.py time_${frame}+DUM.pqr time_${frame}_78.txt time_${frame}_1.txt >> rxn_field.out 

    python ../../free_energy_files/analytic_coloumb.py time_${frame}.pqr CNC CD NE >> coloumb_field.out 2>> /dev/null
    printf "Complete\n" 
    
    done 

printf "\n\t *** Program Complete  *** \n\n" 




