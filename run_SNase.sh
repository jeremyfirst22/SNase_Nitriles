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
if [ ! -f $MOLEC/$fileName ] ; then cp $fileName $MOLEC/ ; fi 

TOP=${PWD}
MDP=$TOP/mdp_files
logFile=$TOP/$MOLEC/$MOLEC.log 
errFile=$TOP/$MOLEC/$MOLEC.err 
FF=$TOP/GMXFF
FORCE_TOOLS=/Users/jeremyfirst/force_calc_tools

check(){
   for var in $@ ; do 
        if [ ! -s $var ] ; then 
            echo ; echo $var missing, exitting... 
            exit 
            fi 
        done 
}

clean(){
    if [ -d amber03.ff ] ; then rm -r amber03.ff *.dat ; fi 
}

create_dir(){
    if [ -z $1 ] ; then 
        echo "ERROR: create_dir requires argument. " ; exit ; fi 

    dirName=$1 
    if [ ! -d $dirName ] ; then mkdir $dirName ; fi 
    
    if [ ! -d $dirName/amber03.ff ] ; then 
        if [ -d $FF/amber03.ff ] ; then 
            cp -r $FF/amber03.ff $dirName/amber03.ff 
            cp $FF/*.dat $dirName/. 
        else 
            echo "FF not found. Standard ff not yet combined with CRO augmented force field." 
            exit 
            fi 
        fi 
}

protein_min(){
    printf "\t\tProtein STEEP................." 
    if [ ! -f Protein_min/$MOLEC.solute_min.gro ] ; then 
        create_dir Protein_min

        cp $MOLEC.pdb Protein_min
        cd Protein_min
        
        #Reads protein .pdb file as input and outputs .pdb as gromacs file and topology file
        #Applies force field amber03 and uses tip3p water molecule
        gmx pdb2gmx -f $MOLEC.pdb -o $MOLEC.gro -p $MOLEC.top -ff amber03 -water tip3p >> $logFile 2>> $errFile
        check $MOLEC.top $MOLEC.gro 
        
        #Reads input with MD parameters, gromacs file, and topology file and creates binary file to speed computation process.
        gmx grompp -f $MDP/solute_steep.mdp -c $MOLEC.gro -p $MOLEC.top -o $MOLEC.solute_min.tpr >> $logFile 2>> $errFile
        check $MOLEC.solute_min.tpr 

        #Runs molecular dynamics simulation 
        gmx mdrun -deffnm $MOLEC.solute_min >> $logFile 2>> $errFile
        check $MOLEC.solute_min.gro 
        
        #Converts minimized .gro file to minimized .pdb file
        gmx editconf -f $MOLEC.solute_min.gro -o $MOLEC.min.pdb >> $logFile 2>> $errFile
        check $MOLEC.min.pdb 
        cp $MOLEC.min.pdb ../

        clean
        printf "Success\n" 
        cd ../

    else 
        printf "Skipped\n" 
        fi 
    check Protein_min/$MOLEC.solute_min.gro 
}

solvate(){
    printf "\t\tSolvating....................." 
    if [ ! -f Solvate/$MOLEC.neutral.gro ] ; then 
        create_dir Solvate

        cp Protein_min/$MOLEC.solute_min.gro Solvate
        cd Solvate
        # Generate a topology from $MOLEC.solute_min.gro with Amber03 force field and tip3p water molecule
        gmx pdb2gmx -f $MOLEC.solute_min.gro -o $MOLEC.gro -p $MOLEC.top -water tip3p -ff amber03 >> $logFile 2>> $errFile 
        check $MOLEC.gro $MOLEC.top 
        
        #Uses octahedron box for solvation and centers molecule
        gmx editconf -f $MOLEC.gro -bt octahedron -box 8 -o temp.centered.gro >> $logFile 2>> $errFile 
        check temp.centered.gro 
        
        #Solvates the molecule with octahedron
        gmx solvate -cp temp.centered.gro -o temp.solvated.gro >> $logFile 2>> $errFile  
        check temp.solvated.gro 
        
        #Generate topology from temp using water tip3p and Amber03s
        gmx pdb2gmx -f temp.solvated.gro -water tip3p -ff amber03 -o temp.$MOLEC.solvated.gro -p temp.$MOLEC.solvated.top >> $logFile 2>> $errFile 
        check temp.$MOLEC.solvated.gro temp.$MOLEC.solvated.top

        #Creates binary file for MD
        gmx grompp -f $MDP/vac_md.mdp -p temp.$MOLEC.solvated.top -c temp.$MOLEC.solvated.gro -o temp.genion.tpr >> $logFile 2>> $errFile 
        check temp.genion.tpr

        #Randomly replaces solvent ions with Cl- and Na+ to neutralize system
        echo SOL | gmx genion -s temp.genion.tpr -neutral -nname 'CL' -pname 'NA' -o temp.neutral.gro >> $logFile 2>> $errFile 
        check temp.neutral.gro 

        #Generate topology of netural molecule
        gmx pdb2gmx -f temp.neutral.gro -water tip3p -ff amber03 -p $MOLEC.neutral.top -o $MOLEC.neutral.gro >> $logFile 2>> $errFile 
        check $MOLEC.neutral.gro 

        clean
        printf "Success\n" 
        cd ../

   else 
       printf "Skipped\n"      
       fi 
   check Solvate/$MOLEC.neutral.gro Solvate/$MOLEC.neutral.top 
}

solvent_min(){
    printf "\t\tSolvent minimization.........." 
    if [ ! -f Solvent_min/$MOLEC.npt_relax.gro ] ; then 
        create_dir Solvent_min

        cp Solvate/$MOLEC.neutral.gro Solvent_min
        cp Solvate/$MOLEC.neutral.top Solvent_min
        cp Solvate/*.itp Solvent_min
        cd Solvent_min
        if [ ! -f $MOLEC.minimize.gro ] ; then
            #Creates binary file for MD 
            gmx grompp -f $MDP/solvent_min.mdp -c $MOLEC.neutral.gro -p $MOLEC.neutral.top -o $MOLEC.minimize.tpr >> $logFile 2>> $errFile 
            check $MOLEC.minimize.tpr

            #Runs MD simulation
            gmx mdrun -deffnm $MOLEC.minimize >> $logFile 2>> $errFile 
            check $MOLEC.minimize.gro
            fi 
        if [ ! -f $MOLEC.nvt_relax.gro ] ; then 
            #constant volume, temperature, protein is frozen
            gmx grompp -f $MDP/solvent_nvt_relax.mdp -c $MOLEC.minimize.gro -p $MOLEC.neutral.top -o $MOLEC.nvt_relax.tpr >> $logFile 2>> $errFile 
            check $MOLEC.nvt_relax.tpr
            
            gmx mdrun -deffnm $MOLEC.nvt_relax >> $logFile 2>> $errFile 
            check $MOLEC.nvt_relax.gro 
            fi 

        if [ ! -f $MOLEC.npt_relax.gro ] ; then
             #constant pressure and temperature 
            gmx grompp -f $MDP/solvent_npt_relax.mdp -c $MOLEC.nvt_relax.gro -p $MOLEC.neutral.top -o $MOLEC.npt_relax.tpr >> $logFile 2>> $errFile 
            check $MOLEC.npt_relax.tpr 

            gmx mdrun -deffnm $MOLEC.npt_relax >> $logFile 2>> $errFile 
            check $MOLEC.npt_relax.gro 
            fi 

        clean
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n" 
        fi 
    check Solvent_min/$MOLEC.npt_relax.gro Solvent_min/$MOLEC.neutral.top 
} 

production_run(){
    printf "\t\tProduction run ..............." 
    if [ ! -f Production/$MOLEC.production.nopbc.gro ] ; then 
        create_dir Production
        cp Solvent_min/$MOLEC.npt_relax.gro Production
        cp Solvent_min/$MOLEC.neutral.top Production
        cp Solvent_min/*.itp Production
        cd Production
        
        if [ ! -f $MOLEC.production.gro ] ; then 
            if [ ! -f $MOLEC.production.tpr ] ; then 
                gmx grompp -f $MDP/production_gfp.mdp -o $MOLEC.production.tpr -p $MOLEC.neutral.top -c $MOLEC.npt_relax.gro >> $logFile 2>> $errFile 
                fi 
            check $MOLEC.production.tpr 
            if [ -f $MOLEC.production.cpt ] ; then 
                gmx mdrun -deffnm $MOLEC.production -cpi $MOLEC.production.cpt >> $logFile 2>> $errFile 
            else 
                gmx mdrun -deffnm $MOLEC.production >> $logFile 2>> $errFile 
                fi 
            fi 
        check $MOLEC.production.gro 
        
        if [ ! -f $MOLEC.production.nopbc.gro ] ; then 
            echo '1 0' | gmx trjconv -f $MOLEC.production.xtc -s $MOLEC.production.tpr -pbc mol -ur compact -center -o $MOLEC.production.nopbc.xtc >> $logFile 2>> $errFile 
            echo '1 0' | gmx trjconv -f $MOLEC.production.gro -s $MOLEC.production.tpr -pbc mol -ur compact -center -o $MOLEC.production.nopbc.gro >> $logFile 2>> $errFile 
            fi 
        check $MOLEC.production.nopbc.gro 

        clean
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n" 
        fi 
}

analyze_hbond(){
    printf "\t\tAnalyzing Hbond content ......" 
    create_dir HBond
    cd HBond
    
    if grep -sq CNF ../Production/$MOLEC.production.nopbc.gro ; then 
        if [ ! -f cnf_num.xvg ] ; then 
            echo "r CNF & a NH" > selection.dat 
            echo "q" >> selection.dat 
            cat selection.dat | gmx make_ndx -f ../Production/$MOLEC.production.gro -o cnf_nh.ndx >> $logFile 2>> $errFile 
            check cnf_nh.ndx 

            echo '18 14 18' | gmx hbond -f ../Production/$MOLEC.production.xtc -s ../Production/$MOLEC.production.tpr -n cnf_nh.ndx -shell 1 -r 0.3 -a 20 -num cnf_num.xvg >> $logFile 2>> $errFile 
            fi 
        check cnf_num.xvg 
        fi 
        
    if [ ! -f cro_num.xvg ] ; then 
        echo "r CROn & a OH or a HO" > selection.dat 
        echo "r CROn & a OH" >> selection.dat
        echo "q" >> selection.dat 
        cat selection.dat | gmx make_ndx -f ../Production/$MOLEC.production.gro -o cro_o_h.ndx >> $logFile 2>> $errFile
        check cro_o_h.ndx

        echo '18 14 19' | gmx hbond -f ../Production/$MOLEC.production.xtc -s ../Production/$MOLEC.production.tpr -n cro_o_h.ndx -shell 1 -r 0.3 -a 20 -num cro_num.xvg >> $logFile 2>> $errFile 
        fi 
    check cro_num.xvg

    clean
    printf "Success\n" 
    cd ../
} 

analyze_hbond_nit(){
    printf "\t\tAnaylzing HBond_nit content..."
    if [ ! -s HBond_nit/$MOLEC.hb_count.xvg ] ; then 
        create_dir HBond_nit
        cd HBond_nit

    if [ ! -d ../Production/amber03.ff ] ; then 
        cp $FF/*.dat ../Production/. 
        cp -r $FF/amber03.ff ../Production/.
        fi 
    check ../Production/amber03.ff/forcefield.itp 

        ## We use veriosn 4.6 of Gromacs for this grompp command, because g_insert_dummy is written for version 4.6
        ## We allow for warnings, since we are generated .tpr from a gromacs 5 mdp file. We are only inserting
        ## atoms this should not matter. 
        if [ ! -f $MOLEC.production.v4.tpr ] ; then 
            grompp -f $MDP/production_gfp.mdp -o $MOLEC.production.v4.tpr -p ../Production/$MOLEC.neutral.top -c ../Production/$MOLEC.npt_relax.gro -maxwarn 3 >> $logFile 2>> $errFile 
            fi
        check $MOLEC.production.v4.tpr 
    
        CT=`grep CNF ../Production/$MOLEC.production.nopbc.gro | grep CT | awk '{print $3}'`
        NH=`grep CNF ../Production/$MOLEC.production.nopbc.gro | grep NH | awk '{print $3}'`
        #echo $CT $NH
        
        if [ ! -s $MOLEC.hb_count.xvg ] ; then  
            $HOME/andrews_gmx/g_nitrile_hbond/g_nitrile_hbond \
                -s $MOLEC.production.v4.tpr \
                -f ../Production/$MOLEC.production.nopbc.xtc \
                -a1 $CT -a2 $NH \
                -select "resname SOL and same residue as within 0.5 of resname CNF and name NH" \
                -o $MOLEC.frame_hb.xvg -op $MOLEC.persistent.xvg \
                -oa $MOLEC.hb_count.xvg -or $MOLEC.geometry.xvg  >> $logFile 2>> $errFile
            check $MOLEC.hb_count.xvg $MOLEC.geometry.xvg $MOLEC.persistent.xvg $MOLEC.frame_hb.xvg 
        fi 
        
        check $MOLEC.hb_count.xvg 
        printf "Success\n" 
        cd ../
    else 
        printf "Skipped\n" 
        fi 
    check HBond_nit/$MOLEC.hb_count.xvg 
}

force_calc(){
    printf "\n\t\tCalculating force:\n" 
    if [[ ! -f force_calc/$MOLEC.solvent_rxn_field.projected.xvg || ! -f force_calc/$MOLEC.external_field.projected.xvg || ! -f force_calc/$MOLEC.total_field.xvg ]] ; then 

    if [ ! -f $FORCE_TOOLS/g_insert_dummy_atom ] ; then 
        printf "\t\t\tERROR: Force tools not found. Skipping force calc\n" 
        return  
        fi 

    create_dir force_calc
    cd force_calc 

    if [ ! -d ../Production/amber03.ff ] ; then 
        cp $FF/*.dat ../Production/. 
        cp -r $FF/amber03.ff ../Production/.
        fi 
    check ../Production/amber03.ff/forcefield.itp 

    ## We use veriosn 4.6 of Gromacs for this grompp command, because g_insert_dummy is written for version 4.6
    ## We allow for warnings, since we are generated .tpr from a gromacs 5 mdp file. We are only inserting
    ## atoms this should not matter. 
    if [ ! -f $MOLEC.production.v4.tpr ] ; then 
        grompp -f $MDP/production_gfp.mdp -o $MOLEC.production.v4.tpr -p ../Production/$MOLEC.neutral.top -c ../Production/$MOLEC.npt_relax.gro -maxwarn 3 >> $logFile 2>> $errFile 
        fi
    check $MOLEC.production.v4.tpr 
    
    CT=`grep CNF ../Production/$MOLEC.production.nopbc.gro | grep CT | awk '{print $3}'`
    NH=`grep CNF ../Production/$MOLEC.production.nopbc.gro | grep NH | awk '{print $3}'`
    #echo $CT $NH
    
    printf "\t\t\tInserting dummy atoms............................" 
    if [ ! -s $MOLEC.with_dummy.xtc ] ; then 
        $FORCE_TOOLS/g_insert_dummy_atom -s $MOLEC.production.v4.tpr -f ../Production/$MOLEC.production.nopbc.xtc -o $MOLEC.with_dummy.xtc -a1 $CT -a2 $NH >> $logFile 2>> $errFile 
        printf "Done\n" 
    else 
        printf "Skipped\n" 
        fi 
    check $MOLEC.with_dummy.xtc

    ## We use the initial configuration so that titration states are conserved (ie, at the end of the production run, pdb2gmx might assign a different titration state to a histidine, which causes it to fail. 
    if [ ! -s $MOLEC.with_dummy.gro ] ; then 
        echo '1 0' | gmx trjconv -s ../Solvent_min/$MOLEC.minimize.tpr -f ../Solvate/$MOLEC.neutral.gro -center -ur compact -pbc mol -o $MOLEC.nopbc.gro >> $logFile 2>> $errFile 
        check $MOLEC.nopbc.gro 

        $FORCE_TOOLS/g_insert_dummy_atom -s $MOLEC.production.v4.tpr -f $MOLEC.nopbc.gro -o $MOLEC.with_dummy.gro -a1 $CT -a2 $NH >> $logFile 2>> $errFile 
        fi 
    check $MOLEC.with_dummy.gro 

    if [ ! -s $MOLEC.with_dummy.top ] ; then 
        gmx pdb2gmx -f $MOLEC.with_dummy.gro -p $MOLEC.with_dummy.top -water tip3p -ff amber03 >> $logFile 2>> $errFile 
        fi 
    check $MOLEC.with_dummy.top 
    
    ##Find new atom numbers 
    CT=`grep CNF $MOLEC.with_dummy.gro | grep CT | awk '{print $3}'`
    NH=`grep CNF $MOLEC.with_dummy.gro | grep NH | awk '{print $3}'`

    echo "[ probe ]" > probe.ndx 
    echo "$CT $NH" >> probe.ndx 

    echo "[ protein ]" > protein.ndx 
    grep -v TCHG $MOLEC.with_dummy.gro | grep -v SOL | grep -v Na | grep -v Cl | tail -n+3 | sed '$d' | awk '{print $3}' >> protein.ndx 

    cp $MOLEC.with_dummy.top $MOLEC.total_field.top 

    if [ ! -s $MOLEC.solvent_rxn_field.top ] ; then 
        $FORCE_TOOLS/zero_charges.py $MOLEC.with_dummy.top protein.ndx $MOLEC.solvent_rxn_field.top >> $logFile 2>> $errFile 
        fi 

    if [ ! -s $MOLEC.external_field.top ] ; then 
        $FORCE_TOOLS/zero_charges.py $MOLEC.with_dummy.top probe.ndx $MOLEC.external_field.top >> $logFile 2>> $errFile 
        fi 
    check $MOLEC.total_field.top $MOLEC.external_field.top $MOLEC.solvent_rxn_field.top 

    for field in total_field external_field solvent_rxn_field ; do 
        printf "\t\t%20s..." "$field" 

        ##Extract forces 
        if [ ! -s $MOLEC.$field.projected.xvg ] ; then 
            printf "forces..." 
            if [ ! -s $MOLEC.$field.xvg ] ; then 
                if [ ! -s $MOLEC.$field.tpr ] ; then 
                    gmx grompp -f $MDP/rerun.mdp -p $MOLEC.$field.top -c $MOLEC.with_dummy.gro -o $MOLEC.$field.tpr  >> $logFile 2>> $errFile 
                    fi 
                check $MOLEC.$field.tpr 
 
                if [ ! -s $MOLEC.$field.trr ] ; then 
                    gmx mdrun -rerun $MOLEC.with_dummy.xtc -s $MOLEC.$field.tpr -deffnm $MOLEC.$field >> $logFile 2>> $errFile 
                    fi 
                check $MOLEC.$field.trr 

                echo 2 | gmx traj -f $MOLEC.$field.trr -s $MOLEC.$field.tpr -of $MOLEC.$field.xvg -xvg none  >> $logFile 2>> $errFile 
                rm $MOLEC.$field.trr 
            fi 
            check $MOLEC.$field.xvg 

            ##extract postions for bond vector
            printf "positions..." 
            if [ ! -s $MOLEC.positions.xvg ] ; then 
                gmx traj -f $MOLEC.with_dummy.xtc -s $MOLEC.$field.tpr -n probe.ndx -ox $MOLEC.positions.xvg -xvg none  >> $logFile 2>> $errFile 
                fi 
            check $MOLEC.positions.xvg 

            ##project force along bond vector 
            printf "Projecting..." 
            $FORCE_TOOLS/get_force.py $MOLEC.positions.xvg $MOLEC.$field.xvg $MOLEC.$field.projected.xvg 
            check $MOLEC.$field.projected.xvg 
            printf "Done\n"  
        else 
            printf "..................................Skipped\n" 
            fi 
        done 
    check $MOLEC.total_field.projected.xvg $MOLEC.external_field.projected.xvg $MOLEC.solvent_rxn_field.projected.xvg 
    clean 
    cd ../Production/
    clean
    cd ../
    
    else 
        printf "\t\t\t\t  ............Skipped\n" 
        fi 
    printf "\n" 
}

chi1_his148(){
    printf "\t\tCalculation chi1 of H148......" 
    if [[ ! -f chi1_his148/$MOLEC.angaver.xvg || ! -f chi1_his148/$MOLEC.angdist.xvg ]] ; then 
        create_dir chi1_his148 
        cd chi1_his148 
    
        N=`grep " N " ../Production/$MOLEC.production.nopbc.gro | grep " 147HIS" | awk '{print$3}'`
        CA=`grep " CA" ../Production/$MOLEC.production.nopbc.gro | grep " 147HIS" | awk '{print$3}'`
        CB=`grep " CB" ../Production/$MOLEC.production.nopbc.gro | grep " 147HIS" | awk '{print$3}'`
        CG=`grep " CG" ../Production/$MOLEC.production.nopbc.gro | grep " 147HIS" | awk '{print$3}'`
        echo "[ chi1_his148 ] " > chi1_his148.ndx 
        echo "$N   $CA   $CB   $CG   " >> chi1_his148.ndx 
        echo " " >> chi1_his148.ndx 
        check chi1_his148.ndx 

        gmx angle -f ../Production/$MOLEC.production.nopbc.xtc -type dihedral -n chi1_his148.ndx -od $MOLEC.angdist.xvg -ov $MOLEC.angaver.xvg >> $logFile 2>> $errFile 

        check $MOLEC.angaver.xvg $MOLEC.angdist.xvg
        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n" 
        fi 
} 

printf "\n\t\t*** Program Beginning $MOLEC***\n\n" 
cd $MOLEC
protein_min
solvate
solvent_min
production_run 
#if grep -sq CNF Production/$MOLEC.production.nopbc.gro ; then 
#    force_calc
#    analyze_hbond_nit
#    fi 
#analyze_hbond
#chi1_his148
#chi1_cnf
#sasa_cro
cd ../

printf "\n\n\t\t*** Program Ending  $MOLEC***\n\n" 
