#!/bin/bash

usage(){
    echo "USAGE: $0 <PDB file {molec.pdb} > "
    exit 
}

if [ -z $1 ] ; then 
    usage 
    fi 

fileName=$1 
if [ ! -f $fileName ] ; then 
    echo "ERROR: $fileName not found " 
    exit 
    fi 
if [[ $fileName == *.pdb ]] ; then 
    MOLEC=$(basename $fileName)
    MOLEC=${MOLEC%.*}
else 
    echo "ERROR: Input file must be PDB file (*.pdb)" 
    exit 
    fi 
if [ ! -d $MOLEC ] ; then mkdir $MOLEC ; fi 
if [ ! -f $MOLEC/$fileName ] ; then cp $fileName $MOLEC/. ; fi 

TOP=${PWD}
MDP=$TOP/mdp_files
logFile=$TOP/$MOLEC/$MOLEC.log 
errFile=$TOP/$MOLEC/$MOLEC.err 
FF=$TOP/GMXFF
forceField=amber03
if [ ! -d $FF/$forceField.ff ] ; then 
    echo ; echo "ERROR: FF not found" 
    exit
    fi  

check(){
    for arg in $@ ; do  
         if [ ! -s $arg ] ; then 
             echo ; echo "ERROR: $arg missing. Exitting" 
             exit 
             fi  
         done 
}

clean(){
    if [ -d $forceField.ff ] ; then rm -r $forceField.ff *.dat ; fi  
}

create_dir(){
    if [ -z $1 ] ; then 
        echo "ERROR: create_dir requires argument. " ; exit ; fi  

    dirName=$1 
    if [ ! -d $dirName ] ; then mkdir $dirName ; fi  
                                                        
    if [ ! -d $dirName/$forceField.ff ] ; then 
        if [ -d $FF/$forceField.ff ] ; then 
            cp -r $FF/$forceField.ff $dirName
            cp $FF/*.dat $dirName/. 
        else 
            echo "FF not found" 
            exit 
            fi  
        fi  
}

protein_steep(){
    printf "\t\tProtein steep............................." 
    if [ ! -f Protein_steep/protein_steep.gro ] ; then 
        create_dir Protein_steep
        
        cp $MOLEC.pdb Protein_steep/.
        cd Protein_steep

        gmx pdb2gmx -f $MOLEC.pdb \
            -p $MOLEC.top \
            -ff $forceField \
            -water tip3p \
            -o $MOLEC.gro >> $logFile 2>> $errFile 
        check $MOLEC.gro 

        echo 'Backbone' | gmx editconf -f $MOLEC.gro \
            -d 1.5 \
            -bt dodecahedron \
            -o boxed.gro >> $logFile 2>> $errFile
        check boxed.gro 

        gmx grompp -f $MDP/protein_steep.mdp \
            -c boxed.gro \
            -p $MOLEC.top \
            -o protein_steep.tpr >> $logFile 2>> $errFile 
        check protein_steep.tpr 

        gmx mdrun -deffnm protein_steep >> $logFile 2>> $errFile 
        check protein_steep.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi
} 

solvate(){
    printf "\t\tSolvating protein........................." 
    if [ ! -f Solvate/neutral.top ] ; then 
        create_dir Solvate
        
        cp Protein_steep/protein_steep.gro Solvate/. 
        cp Protein_steep/$MOLEC.top Solvate/. 
        cp Protein_steep/$MOLEC*.itp Solvate/. 
        cd Solvate

        gmx solvate -cp protein_steep.gro \
            -p $MOLEC.top \
            -o solvated.gro >> $logFile 2>> $errFile 
        check solvated.gro

        gmx grompp -f $MDP/vac_md.mdp \
            -p $MOLEC.top \
            -c solvated.gro \
            -o genion.tpr >> $logFile 2>> $errFile 
        check genion.tpr
        
        echo 'SOL' | gmx genion -s genion.tpr \
            -neutral \
            -nname 'CL' \
            -pname 'NA' \
            -o neutral.gro >> $logFile 2>> $errFile 
        check neutral.gro 

        gmx pdb2gmx -f neutral.gro \
            -ff $forceField \
            -water tip3p \
            -p neutral.top \
            -o neutral.gro >> $logFile 2>> $errFile 
        check neutral.top 

        sed 's/POSRES/POSRES_IONS/' neutral_Ion2.itp > temp.itp 
        mv temp.itp neutral_Ion2.itp 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

solvent_steep(){
    printf "\t\tSolvent steep............................." 
    if [ ! -f Solvent_steep/solvent_steep.gro ] ; then 
        create_dir Solvent_steep
        
        cp Solvate/neutral.gro Solvent_steep/. 
        cp Solvate/neutral.top Solvent_steep/. 
        cp Solvate/neutral*.itp Solvent_steep/. 
        cp Solvate/posre*.itp Solvent_steep/. 
        cd Solvent_steep

        gmx grompp -f $MDP/solvent_steep.mdp \
            -p neutral.top \
            -c neutral.gro \
            -o solvent_steep.tpr >> $logFile 2>> $errFile 
        check solvent_steep.tpr 

        gmx mdrun -deffnm solvent_steep >> $logFile 2>> $errFile 
        check solvent_steep.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

solvent_nvt(){
    printf "\t\tSolvent NVT relaxation...................." 
    if [ ! -f Solvent_nvt/solvent_nvt.gro ] ; then 
        create_dir Solvent_nvt
        
        cp Solvent_steep/solvent_steep.gro Solvent_nvt/. 
        cp Solvent_steep/neutral.top Solvent_nvt/. 
        cp Solvent_steep/*.itp Solvent_nvt/. 
        cd Solvent_nvt

        gmx grompp -f $MDP/solvent_nvt_relax.mdp \
            -c solvent_steep.gro \
            -p neutral.top \
            -o solvent_nvt.tpr >> $logFile 2>> $errFile 
        check solvent_nvt.tpr 

        gmx mdrun -deffnm solvent_nvt >> $logFile 2>> $errFile 
        check solvent_nvt.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

solvent_npt(){
    printf "\t\tSolvent NPT relaxation...................." 
    if [ ! -f Solvent_npt/solvent_npt.gro ] ; then 
        create_dir Solvent_npt
        
        cp Solvent_nvt/solvent_nvt.gro Solvent_npt/. 
        cp Solvent_nvt/neutral.top Solvent_npt/. 
        cp Solvent_nvt/*.itp Solvent_npt/. 
        cd Solvent_npt

        gmx grompp -f $MDP/solvent_npt_relax.mdp \
            -c solvent_nvt.gro \
            -p neutral.top \
            -o solvent_npt.tpr >> $logFile 2>> $errFile 
        check solvent_npt.tpr 

        gmx mdrun -deffnm solvent_npt >> $logFile 2>> $errFile 
        check solvent_npt.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

production(){
    printf "\t\tProduction run............................" 
    if [ ! -f Production/$MOLEC.nopbc.gro ] ; then 
        create_dir Production
        
        cp Solvent_npt/neutral.top Production/.
        cp Solvent_npt/solvent_npt.gro Production/.
        cp Solvent_npt/*.itp Production/. 
        cd Production

        if [ ! -f $MOLEC.gro ] ; then 
            if [ ! -f $MOLEC.tpr ] ; then 
                gmx grompp -f $MDP/production.mdp \
                    -p neutral.top \
                    -c solvent_npt.gro \
                    -o $MOLEC.tpr >> $logFile 2>> $errFile 
                fi 
                check $MOLEC.tpr 

            if [ -f $MOLEC.cpt ] ; then 
                gmx mdrun -deffnm $MOLEC -cpi $MOLEC.cpt >> $logFile 2>> $errFile  
            else 
                gmx mdrun -deffnm $MOLEC >> $logFile 2>> $errFile 
                fi 
            fi 
        check $MOLEC.gro 

        if [ ! -f $MOLEC.nopbc.xtc ] ; then 
            echo 'Protein System' | gmx trjconv -f $MOLEC.xtc \
                -center \
                -s $MOLEC.tpr \
                -ur compact \
                -pbc mol \
                -o $MOLEC.nopbc.xtc >> $logFile 2>> $errFile 
            fi 
        check $MOLEC.nopbc.xtc 

        if [ ! -f $MOLEC.nopbc.gro ] ; then 
            echo 'Protein System' | gmx trjconv -f $MOLEC.gro \
                -center \
                -s $MOLEC.tpr \
                -ur compact \
                -pbc mol \
                -o $MOLEC.nopbc.gro >> $logFile 2>> $errFile 
            fi 
        check $MOLEC.nopbc.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

hbond(){
    printf "\t\tAnalyzing hydrogen bonds.................." 
    if [ ! -f hbond/all_hbnum.xvg ] ; then 
        create_dir hbond
        cp Production/solvent_npt.gro hbond/. 
        cp Production/neutral.top hbond/. 
        cp Production/*.itp hbond/. 
        cd hbond
        
        ##Version 4 tpr required for Andrew's nitrile code
        grompp -f $MDP/vac_md.mdp \
            -c solvent_npt.gro \
            -p neutral.top \
            -maxwarn 3 \
            -o v4.tpr >> $logFile 2>> $errFile 
        check v4.tpr 

        CD=`grep CNC solvent_npt.gro | grep CD | awk '{print $3}'`
        NE=`grep CNC solvent_npt.gro | grep NE | awk '{print $3}'`
        ~/andrews_gmx/g_nitrile_hbond/g_nitrile_hbond -f ../Production/$MOLEC.xtc \
            -s v4.tpr \
            -select 'resname SOL and same residue as within 0.5 of resname CNC and name NE' \
            -a1 $CD \
            -a2 $NE \
            -op persistent.xvg \
            -or geometry.xvg \
            -o frame_hb.xvg >> $logFile 2>> $errFile 
        check frame_hb.xvg geometry.xvg persistent.xvg 

        clean
        rm solvent_npt.gro neutral.top *.itp 
   
        echo "r CNC & a NE" > selection.dat 
        echo "r SOL" >> selection.dat 
        echo "!r CNC" >> selection.dat  
        echo "q" >> selection.dat  

        touch empty.ndx 
        cat selection.dat | gmx make_ndx -f ../Production/solvent_npt.gro \
            -n empty.ndx \
            -o index.ndx >> $logFile 2>> $errFile 
        check index.ndx 

        echo '0 1 0' | gmx hbond -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -shell 1.0 \
            -ac wat_hbac.xvg \
            -dist wat_hbdist.xvg \
            -ang wat_hbang.xvg \
            -life wat_hblife.xvg \
            -num wat_hbnum.xvg >> $logFile 2>> $errFile 
        check wat_hbnum.xvg 

        echo '0 2 0' | gmx hbond -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -shell 1.0 \
            -ac all_hbac.xvg \
            -dist all_hbdist.xvg \
            -ang all_hbang.xvg \
            -life all_hblife.xvg \
            -num all_hbnum.xvg >> $logFile 2>> $errFile 
        check all_hbnum.xvg 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

minimage(){
    printf "\t\tCalculating minimum image................." 
    if [ ! -f minimage/mindist.xvg ] ; then 
        create_dir minimage
        cd minimage
        clean 
        
        echo 'Protein' | gmx mindist -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -pi \
            -od mindist.xvg >> $logFile 2>> $errFile 
            check mindist.xvg 

        printf "Success\n" 
        cd ../ 
    else
        printf "Skipped\n"
        fi  
} 

rmsd(){
    printf "\t\tCalculating RMSD.........................." 
    if [ ! -f rmsd/without_ter.xvg ] ; then 
        create_dir rmsd
        cd rmsd
        clean 

        echo 'Backbone Backbone' | gmx rms -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -o backbone.xvg >> $logFile 2>> $errFile 
        check backbone.xvg 

        echo '"Backbone" & ri 6-130' > selection.dat 
        echo "q" >> selection.dat 

        cat selection.dat | gmx make_ndx -f ../Production/solvent_npt.gro \
            -o index.ndx >> $logFile 2>> $errFile 
        check index.ndx 
        
        echo "Backbone_&_r_6_130 Backbone_&_r_6_130" | gmx rms -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -o without_ter.xvg >> $logFile 2>> $errFile 
        check without_ter.xvg

        printf "Success\n" 
        cd ../ 
    else
        printf "Skipped\n"
        fi  
}

chi(){
    printf "\t\tCalculating phi psi chi dihedrals........." 
    if [ ! -f chi/order.xvg ] ; then 
        create_dir chi
        cd chi
        clean 

        gmx chi -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -norad \
            -rama >> $logFile 2>> $errFile 
        check order.xvg 

        printf "Success\n" 
        cd ../ 
    else
        printf "Skipped\n"
        fi  
} 

rgyr(){
    printf "\t\tCalculating radius of gyration............" 
    if [ ! -f rgyr/gyrate.xvg ] ; then 
        create_dir rgyr
        cd rgyr
        clean 

        echo 'Backbone' | gmx gyrate -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -o backbone.xvg >> $logFile 2>> $errFile 
        check backbone.xvg 

        echo '"Backbone" & ri 6-130' > selection.dat 
        echo "q" >> selection.dat 

        cat selection.dat | gmx make_ndx -f ../Production/solvent_npt.gro \
            -o index.ndx >> $logFile 2>> $errFile 
        check index.ndx 
        
        echo "Backbone_&_r_6_130" | gmx gyrate -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -o without_ter.xvg >> $logFile 2>> $errFile 
        check without_ter.xvg

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

printf "\n\t\t*** Program Beginning ***\n\n" 
cd $MOLEC
protein_steep
solvate
solvent_steep
solvent_nvt
solvent_npt
production 
if grep -sq CNC $MOLEC.pdb ; then 
    hbond 
    fi 
minimage
rmsd 
chi
rgyr
cd ../

printf "\n\n\t\t*** Program Ending    ***\n\n" 
