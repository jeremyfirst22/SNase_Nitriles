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
FORCE_TOOLS=/Users/jfirst/force_calc_tools
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
    if [ ! -f Production/$MOLEC.gro ] ; then
        create_dir Production

        cp Solvent_npt/neutral.top Production/.
        cp Solvent_npt/solvent_npt.gro Production/.
        cp Solvent_npt/*.itp Production/.
        cd Production

        maxTime=20  

        if [ ! -f $MOLEC.gro ] ; then
            if [ ! -f $MOLEC.tpr ] ; then
                gmx grompp -f $MDP/production.mdp \
                    -p neutral.top \
                    -c solvent_npt.gro \
                    -o $MOLEC.tpr >> $logFile 2>> $errFile
                fi  
                check $MOLEC.tpr

                if [ ! -f $maxTime.tpr ] ; then
                    gmx convert-tpr -s $MOLEC.tpr \
                        -until $((maxTime*1000)) \
                        -o $maxTime.tpr >> $logFile 2>> $errFile
                    fi  
                check $maxTime.tpr

            if [ -f $MOLEC.cpt ] ; then 
#                ibrun mdrun_mpi -s $maxTime.tpr \
                gmx mdrun -s $maxTime.tpr \
                    -deffnm $MOLEC \
                    -cpi $MOLEC.cpt >> $logFile 2>> $errFile
            else 
                #ibrun mdrun_mpi -s $maxTime.tpr \
                gmx mdrun -s $maxTime.tpr \
                    -deffnm $MOLEC >> $logFile 2>> $errFile
                fi  
            fi  
        check $MOLEC.gro 

        clean
        printf "Success\n" 
        cd ../ 
    else
        printf "Skipped\n"
        fi  
} 

nopbc(){
    printf "\t\tCreating nopbc copy......................."
    if [ ! -f nopbc/small.xtc ] ; then
        create_dir nopbc
        cd nopbc

        echo "[ protein ]" > protein.ndx 
        grep -v -e HOH -e SOL -e CL -e NA ../Production/$MOLEC.gro | tail -n+3 | sed '$d' | awk '{print $3}' >> protein.ndx 

        if [ ! -f all.ndx ] ; then 
            gmx select -f ../Production/$MOLEC.gro \
                -select "all" \
                -on all.ndx >> $logFile 2>> $errFile 
        fi 
        check all.ndx 

        cat protein.ndx > both.ndx 
        cat all.ndx >> both.ndx 


        if [ ! -f interfacial.ndx ] ; then 
            echo "ri 39-42" > selection.dat 
            echo "q" >> selection.dat 
            cat selection.dat | gmx make_ndx -f ../Production/$MOLEC.gro \
                -o interfacial.ndx >> $logFile 2>> $errFile 
        fi 
        check interfacial.ndx 

        if [ ! -f whole.xtc ] ; then 
            echo 'System' | gmx trjconv -s ../Production/$MOLEC.tpr \
                -f ../Production/$MOLEC.xtc \
                -pbc nojump \
                -ur compact \
                -o whole.xtc >> $logFile 2>> $errFile 
        fi 
        check whole.xtc 

        if [ ! -f centered.xtc ] ; then 
            echo 'r_39-42 System' | gmx trjconv -s ../Production/$MOLEC.tpr \
                -f whole.xtc \
                -n interfacial.ndx \
                -pbc mol \
                -ur compact \
                -center \
                -o centered.xtc >> $logFile 2>> $errFile 
        fi 
        check centered.xtc 

        if [ ! -f $MOLEC.nopbc.xtc ] ; then  
            echo 'protein all' | gmx trjconv -s ../Production/$MOLEC.tpr \
                -f centered.xtc \
                -n both.ndx \
                -pbc mol \
                -ur compact \
                -center \
                -o $MOLEC.nopbc.xtc >> $logFile 2>> $errFile 
        fi 
        check $MOLEC.nopbc.xtc 

        if [ ! -f $MOLEC.nopbc.gro ] ; then 
            echo 'r_39-42 System' | gmx trjconv -s ../Production/$MOLEC.tpr \
                -f ../Production/solvent_npt.gro \
                -n interfacial.ndx \
                -center \
                -pbc mol \
                -ur compact \
                -o $MOLEC.nopbc.gro >> $logFile 2>> $errFile 
        fi
        check $MOLEC.nopbc.gro 

        ## Check to see if Ca ion in same peridodic image as protein
        echo "CA Protein" | gmx mindist -s $MOLEC.nopbc.gro \
            -f $MOLEC.nopbc.xtc \
            -pbc no \
            -od mindist.xvg >> $logFile 2>> $errFile 
        check mindist.xvg 

        rm whole.xtc centered.xtc 

        if [ ! -f small.xtc ] ; then 
            echo 'System' | gmx trjconv -s ../Production/$MOLEC.tpr \
                -f $MOLEC.nopbc.xtc \
                -dt 100 \
                -o small.xtc >> $logFile 2>> $errFile 
        fi 
        check small.xtc 

        clean
        printf "Success\n"
        cd ../
    else
        printf "Skipped\n"
        fi
}

hbond(){
    printf "\t\tAnalyzing hydrogen bonds.................." 
    if [ ! -f hbond/geometry.xvg ] ; then 
        create_dir hbond
        cp Production/solvent_npt.gro hbond/. 
        cp Production/neutral.top hbond/. 
        cp Production/*.itp hbond/. 
        cd hbond
        
        ##Version 4 tpr required for custom built nitrile code
        grompp -f $MDP/vac_md.mdp \
            -c solvent_npt.gro \
            -p neutral.top \
            -maxwarn 3 \
            -o v4.tpr >> $logFile 2>> $errFile 
        check v4.tpr 

        CT=`grep CNC solvent_npt.gro | grep CD | awk '{print $3}'`
        NH=`grep CNC solvent_npt.gro | grep NE | awk '{print $3}'`
        source /usr/local/gromacs/bin/GMXRC
        g_nitrile_hbond -f ../Production/$MOLEC.xtc \
            -s v4.tpr \
            -select 'not resname CNC and (same residue as within 0.5 of resname CNC and name NE)' \
            -a1 $CT \
            -a2 $NH \
            -op persistent.xvg \
            -or geometry.xvg \
            -oa hb_count.xvg \
            -onwr nw_geometry.xvg \
            -o frame_hb.xvg >> $logFile 2>> $errFile 
        check frame_hb.xvg geometry.xvg persistent.xvg 

        clean
        rm solvent_npt.gro neutral.top *.itp 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

hbond_1(){
    printf "\t\tAnalyzing residence time with forgive 1..." 
    if [ ! -f hbond_1/persistent.xvg ] ; then 
        create_dir hbond_1
        cp Production/solvent_npt.gro hbond_1/. 
        cp Production/neutral.top hbond_1/. 
        cp Production/*.itp hbond_1/. 
        cd hbond_1
        
        ##Version 4 tpr required for custom built nitrile code
        grompp -f $MDP/vac_md.mdp \
            -c solvent_npt.gro \
            -p neutral.top \
            -maxwarn 3 \
            -o v4.tpr >> $logFile 2>> $errFile 
        check v4.tpr 

        CT=`grep CNC solvent_npt.gro | grep CD | awk '{print $3}'`
        NH=`grep CNC solvent_npt.gro | grep NE | awk '{print $3}'`
        source /usr/local/gromacs/bin/GMXRC
        g_nitrile_hbond -f ../Production/$MOLEC.xtc \
            -s v4.tpr \
            -select 'not resname CNC and (same residue as within 0.5 of resname CNC and name NE)' \
            -a1 $CT \
            -a2 $NH \
            -forgiveness 1 \
            -op persistent.xvg \
            -o frame_hb.xvg >> $logFile 2>> $errFile 
        check persistent.xvg 

        clean
        rm solvent_npt.gro neutral.top *.itp 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

hbond_2(){
    printf "\t\tAnalyzing residence time with forgive 2..." 
    if [ ! -f hbond_2/persistent.xvg ] ; then 
        create_dir hbond_2
        cp Production/solvent_npt.gro hbond_2/. 
        cp Production/neutral.top hbond_2/. 
        cp Production/*.itp hbond_2/. 
        cd hbond_2
        
        ##Version 4 tpr required for custom built nitrile code
        grompp -f $MDP/vac_md.mdp \
            -c solvent_npt.gro \
            -p neutral.top \
            -maxwarn 3 \
            -o v4.tpr >> $logFile 2>> $errFile 
        check v4.tpr 

        CT=`grep CNC solvent_npt.gro | grep CD | awk '{print $3}'`
        NH=`grep CNC solvent_npt.gro | grep NE | awk '{print $3}'`
        source /usr/local/gromacs/bin/GMXRC
        g_nitrile_hbond -f ../Production/$MOLEC.xtc \
            -s v4.tpr \
            -select 'not resname CNC and (same residue as within 0.5 of resname CNC and name NE)' \
            -a1 $CT \
            -a2 $NH \
            -forgiveness 2 \
            -op persistent.xvg \
            -o frame_hb.xvg >> $logFile 2>> $errFile 
        check persistent.xvg 

        clean
        rm solvent_npt.gro neutral.top *.itp 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

fit_hbond(){
    printf "\t\tFitting hbond analysia...................."
    if [ ! -f hbond/geometry.xvg ] ; then
        printf "Skipping\n"
        break
        fi

    if [ ! -f ~/normal_distribution/tiltAngle ] ; then
        printf "Skipping\n"
        break
        fi

    if [[ ! -f fit_hbond/dist.poly || ! -f fit_hbond/theta1.poly || ! -f fit_hbond/theta2.poly ]] ; then
        create_dir fit_hbond
        cp hbond/*geometry.xvg fit_hbond/.
        cd fit_hbond
        clean

        check geometry.xvg nw_geometry.xvg
        cat geometry.xvg | grep -v "^[#@]" | awk '{print $3}' > clean_dist.dat
        cat nw_geometry.xvg | grep -v "^[#@]" | awk '{print $3}' >> clean_dist.dat
        check clean_dist.dat

        /Users/jfirst/normal_distribution/tiltAngle -f clean_dist.dat -o dist.his -g dist.gaus -p dist.poly
        check dist.poly

        check geometry.xvg nw_geometry.xvg
        cat geometry.xvg | grep -v "^[#@]" | awk '{print $4}' > clean_theta1.dat
        cat nw_geometry.xvg | grep -v "^[#@]" | awk '{print $4}' >> clean_theta1.dat
        check clean_theta1.dat

        /Users/jfirst/normal_distribution/tiltAngle -f clean_theta1.dat -o theta1.his -g theta1.gaus -p theta1.poly -n 81
        check theta1.poly

        check geometry.xvg nw_geometry.xvg
        cat geometry.xvg | grep -v "^[#@]" | awk '{print $5}' > clean_theta2.dat
        cat nw_geometry.xvg | grep -v "^[#@]" | awk '{print $5}' >> clean_theta2.dat
        check clean_theta2.dat

        /Users/jfirst/normal_distribution/tiltAngle -f clean_theta2.dat -o theta2.his -g theta2.gaus -p theta2.poly
        check theta2.poly

        printf "Success\n"
        cd ../
    else
        printf "Skipped\n"
        fi
}

sasa(){
    printf "\t\tAnalyzing SASA of CNF....................."
    if [ ! -f sasa/nitrile_area.xvg ] ; then
        create_dir sasa
        cd sasa

        if [ ! -s cnc_area.xvg ] ; then
            gmx sasa -f ../Production/$MOLEC.xtc \
                -s ../Production/$MOLEC.tpr \
                -surface 'Protein' \
                -output 'resname CNC' \
                -ndots 240 \
                -or cnc_resarea.xvg \
                -o cnc_area.xvg >> $logFile 2>> $errFile
            fi
        check cnc_area.xvg

        if [ ! -s sidechain_area.xvg ] ; then
            gmx sasa -f ../Production/$MOLEC.xtc \
                -s ../Production/$MOLEC.tpr \
                -surface 'Protein' \
                -output 'resname CNC and (name NE CD SG CB)' \
                -ndots 240 \
                -o sidechain_area.xvg >> $logFile 2>> $errFile
            fi
        check sidechain_area.xvg

        if [ ! -s thiocyanate.xvg ] ; then
            gmx sasa -f ../Production/$MOLEC.xtc \
                -s ../Production/$MOLEC.tpr \
                -surface 'Protein' \
                -output 'resname CNC and (name NE CD SG)' \
                -ndots 240 \
                -o thiocyanate.xvg >> $logFile 2>> $errFile
            fi
        check thiocyanate.xvg

        if [ ! -s nitrile_area.xvg ] ; then
            gmx sasa -f ../Production/$MOLEC.xtc \
                -s ../Production/$MOLEC.tpr \
                -surface 'Protein' \
                -output 'resname CNC and (name NE CD)' \
                -ndots 240 \
                -o nitrile_area.xvg >> $logFile 2>> $errFile
            fi
        check nitrile_area.xvg

        clean
        printf "Success\n"
        cd ../
    else
        printf "Skipped\n"
        fi
}

force_calc(){
    printf "\n\t\tCalculating force:\n" 
    if [[ ! -f force_calc_ca/$MOLEC.solvent_rxn_field.projected.xvg || ! -f force_calc_ca/$MOLEC.external_field.projected.xvg || ! -f force_calc_ca/$MOLEC.total_field.projected.xvg || ! -f force_calc_ca/$MOLEC.protein_field.projected.xvg ]] ; then 

        if [ ! -f $FORCE_TOOLS/g_insert_dummy_atom ] ; then 
            printf "\t\t\tERROR: Force tools not found. Skipping force calc\n" 
            return  
            fi 

        create_dir force_calc_ca
        cp Production/*.itp force_calc_ca/. 
        cp Production/neutral.top force_calc_ca/. 
        cp Production/solvent_npt.gro force_calc_ca/. 
        cp Solvate/neutral.gro force_calc_ca/. 
        cd force_calc_ca 

        ## We use version 4.6 of Gromacs for this grompp command, because g_insert_dummy is written for version 4.6
        ## We allow for warnings, since we are generated .tpr from a gromacs 5 mdp file. We are only inserting
        ## atoms this should not matter. 
        grompp -f $MDP/vac_md.mdp \
            -c solvent_npt.gro \
            -p neutral.top \
            -maxwarn 3 \
            -o v4.tpr >> $logFile 2>> $errFile 
        check v4.tpr 
        
        CD=`grep CNC ../Production/solvent_npt.gro | grep CD | awk '{print $3}'`
        NE=`grep CNC ../Production/solvent_npt.gro | grep NE | awk '{print $3}'`
        #echo $CD $NE
        
        printf "\t\t\tInserting dummy atoms............................" 
        if [[ ! -s $MOLEC.with_dummy.top || ! -s $MOLEC.with_dummy.xtc ]] ; then 
            if [ ! -s $MOLEC.with_dummy.xtc ] ; then 
                source /usr/local/gromacs/bin/GMXRC
                $FORCE_TOOLS/g_insert_dummy_atom -f ../nopbc/$MOLEC.nopbc.xtc \
                    -s v4.tpr \
                    -a1 $CD \
                    -a2 $NE \
                    -o $MOLEC.with_dummy.xtc >> $logFile 2>> $errFile 
            fi 
            check $MOLEC.with_dummy.xtc

            if [ ! -s $MOLEC.with_dummy.gro ] ; then 
                source /usr/local/gromacs/bin/GMXRC
                $FORCE_TOOLS/g_insert_dummy_atom -f ../nopbc/$MOLEC.nopbc.gro \
                    -s v4.tpr \
                    -a1 $CD \
                    -a2 $NE \
                    -o $MOLEC.with_dummy.gro >> $logFile 2>> $errFile 
                fi 
            check $MOLEC.with_dummy.gro 

            if [ ! -s $MOLEC.with_dummy.top ] ; then 
                gmx pdb2gmx -f $MOLEC.with_dummy.gro \
                    -water tip3p \
                    -ff amber03 \
                    -p $MOLEC.with_dummy.top \
                    -o $MOLEC.with_dummy.gro >> $logFile 2>> $errFile 
                fi 
            check $MOLEC.with_dummy.top 
            printf "Done\n" 
        else 
            printf "Skipped\n" 
            fi 
        
        ##Find new atom numbers 
        CD=`grep CNC $MOLEC.with_dummy.gro | grep CD | awk '{print $3}'`
        NE=`grep CNC $MOLEC.with_dummy.gro | grep NE | awk '{print $3}'`

        echo "[ probe ]" > probe.ndx 
        echo "$CD $NE" >> probe.ndx 

        echo "[ protein ]" > protein.ndx 
        #grep -v -e TCHG -e SOL -e HOH -e NA -e CL -e 144CA $MOLEC.with_dummy.gro | tail -n+3 | sed '$d' | awk '{print $3}' >> protein.ndx 
        grep -v -e TCHG -e SOL -e HOH -e NA -e CL $MOLEC.with_dummy.gro | tail -n+3 | sed '$d' | awk '{print $3}' >> protein.ndx 

        echo "[ solvent ]" > solvent.ndx 
        #grep $MOLEC.with_dummy.gro -e SOL -e HOH -e NA -e CL -e 144CA | cut -c 16-20 >> solvent.ndx 
        grep $MOLEC.with_dummy.gro -e SOL -e HOH -e NA -e CL | cut -c 16-20 >> solvent.ndx 


        ###
        ### Create topologies
        ###
        cp $MOLEC.with_dummy.top $MOLEC.total_field.top 

        if [ ! -s $MOLEC.solvent_rxn_field.top ] ; then 
            $FORCE_TOOLS/zero_charges.py $MOLEC.with_dummy.top protein.ndx $MOLEC.solvent_rxn_field.top >> $logFile 2>> $errFile 
            fi 

        if [ ! -s $MOLEC.external_field.top ] ; then 
            $FORCE_TOOLS/zero_charges.py $MOLEC.with_dummy.top probe.ndx $MOLEC.external_field.top >> $logFile 2>> $errFile 
            fi 

        if [ ! -s $MOLEC.protein_field.top ] ; then 
            $FORCE_TOOLS/zero_charges.py $MOLEC.external_field.top solvent.ndx $MOLEC.protein_field.top >> $logFile 2>> $errFile 

            sed "s/tip3p.itp/zerocharge.itp/" $MOLEC.protein_field.top > temp.top 
            mv temp.top $MOLEC.protein_field.top 

            head -n 6 amber03.ff/tip3p.itp > zerocharge.itp 
            echo "  1   OW          1       SOL       OW       1       0.000    16.00000" >> zerocharge.itp 
            echo "  2   HW          1       SOL       HW1      1       0.000     1.00800" >> zerocharge.itp 
            echo "  3   HW          1       SOL       HW2      1       0.000     1.00800" >> zerocharge.itp 
            tail -n 25 amber03.ff/tip3p.itp >> zerocharge.itp 

            mv zerocharge.itp amber03.ff/.
            fi 
        check $MOLEC.total_field.top $MOLEC.external_field.top $MOLEC.solvent_rxn_field.top $MOLEC.protein_field.top

        for field in total_field external_field solvent_rxn_field protein_field ; do 
            printf "\t\t%20s..." $field 

            ##Extract forces 
            if [ ! -s $MOLEC.$field.projected.xvg ] ; then 
                printf "forces..." 
                if [ ! -s $MOLEC.$field.xvg ] ; then 
                    if [ ! -s $MOLEC.$field.tpr ] ; then 
                        gmx grompp -f $MDP/rerun.mdp \
                            -p $MOLEC.$field.top \
                            -c $MOLEC.with_dummy.gro \
                            -o $MOLEC.$field.tpr  >> $logFile 2>> $errFile 
                        fi 
                    check $MOLEC.$field.tpr 
 
                    if [ ! -s $MOLEC.$field.trr ] ; then 
                        gmx mdrun -s $MOLEC.$field.tpr \
                            -rerun $MOLEC.with_dummy.xtc \
                            -deffnm $MOLEC.$field >> $logFile 2>> $errFile 
                        fi 
                    check $MOLEC.$field.trr 

                    echo 2 | gmx traj -f $MOLEC.$field.trr \
                        -s $MOLEC.$field.tpr \
                        -xvg none \
                        -of $MOLEC.$field.xvg >> $logFile 2>> $errFile 
                    rm $MOLEC.$field.trr 
                fi 
                check $MOLEC.$field.xvg 

                ##extract postions for bond vector
                printf "positions..." 
                if [ ! -s $MOLEC.positions.xvg ] ; then 
                    gmx traj -f $MOLEC.with_dummy.xtc \
                        -s $MOLEC.$field.tpr \
                        -n probe.ndx \
                        -xvg none \
                        -ox $MOLEC.positions.xvg >> $logFile 2>> $errFile 
                    fi 
                check $MOLEC.positions.xvg 

                ##project force along bond vector 
                printf "Projecting..." 
                $FORCE_TOOLS/get_force.py $MOLEC.positions.xvg $MOLEC.$field.xvg $MOLEC.$field.projected.xvg 
                check $MOLEC.$field.projected.xvg 

                ~/normal_distribution/tiltAngle -f $MOLEC.$field.projected.xvg -o $field.out -g $field.gaus -p $field.poly -t 25 >> $logFile 2>> $errFile 
                check $field.gaus $field.poly

                printf "Done\n"  
            else 
                printf "..................................Skipped\n" 
                fi 
            done 

        for field in external-srf ; do 
            printf "\t\t%20s..." $field 

            ##Extract forces 
            if [ ! -s $MOLEC.$field.projected.xvg ] ; then 
                echo "import numpy as np" > $field.py 
                echo "srf = np.genfromtxt('$MOLEC.solvent_rxn_field.projected.xvg')" >> $field.py 
                echo "ext = np.genfromtxt('$MOLEC.external_field.projected.xvg')" >> $field.py 
                echo "np.savetxt('$MOLEC.$field.projected.xvg',ext - srf)" >> $field.py 
                python $field.py 
                check $MOLEC.$field.projected.xvg
                

                ~/normal_distribution/tiltAngle -f $MOLEC.$field.projected.xvg -o $field.out -g $field.gaus -p $field.poly -t 25 >> $logFile 2>> $errFile 
                check $field.gaus $field.poly

                printf "..................................Done\n"  
            else 
                printf "..................................Skipped\n" 
                fi 
            done 
        check $MOLEC.total_field.projected.xvg $MOLEC.external_field.projected.xvg $MOLEC.solvent_rxn_field.projected.xvg 

        #rm $MOLEC.with_dummy.xtc  ##Extra 4 GB per trajectory that isn't needed anymore. 

        clean 
        cd ../
    else 
        printf "\t\t\t\t  ........................Skipped\n" 
        fi 
    printf "\n" 
}

rmsd(){
    printf "\t\tCalculating RMSD.........................."
    if [ ! -f rmsd/crystal.xvg ] ; then
        create_dir rmsd
        cd rmsd

        if  [ ! -f backbone.xvg ] ; then 
            echo 'Backbone Backbone' | gmx rms -f ../Production/$MOLEC.xtc \
                -s ../Production/$MOLEC.tpr \
                -b 10000 \
                -o backbone.xvg >> $logFile 2>> $errFile
        fi 
        check backbone.xvg

        if [ ! -f crystal.xvg ] ; then 
            echo "4 & ri 7-141" > selection.dat 
            echo "q" >> selection.dat 

            cat selection.dat | gmx make_ndx -f ../Production/$MOLEC.gro \
                -o crystal.ndx >> $logFile 2>> $errFile 
            check crystal.ndx 

            echo 'Backbone_&_r_7-141 Backbone_&_r_7-141' | gmx rms -f ../Production/$MOLEC.xtc \
                -s ../Production/$MOLEC.tpr \
                -n crystal.ndx \
                -o crystal.xvg >> $logFile 2>> $errFile 
        fi 
        check crystal.xvg 

        clean
        printf "Success\n"
        cd ../
    else
        printf "Skipped\n"
        fi
}

chi(){
    printf "\t\tCalculating Chi1 and Chi2................."
    if [ ! -f chi/chi2.xvg ] ; then
        create_dir chi
        cd chi

        indexN=`cat ../Production/$MOLEC.gro | grep CNC | grep " N " | awk '{print $3}'`
        indexCA=`cat ../Production/$MOLEC.gro | grep CNC | grep " CA " | awk '{print $3}'`
        indexCB=`cat ../Production/$MOLEC.gro | grep CNC | grep " CB " | awk '{print $3}'`
        indexSG=`cat ../Production/$MOLEC.gro | grep CNC | grep " SG " | awk '{print $3}'`
        indexCD=`cat ../Production/$MOLEC.gro | grep CNC | grep " CD " | awk '{print $3}'`

        echo [ chi1 ] > chi1.ndx 
        echo "$indexN $indexCA $indexCB $indexSG" >> chi1.ndx 

        echo [ chi2 ] > chi2.ndx 
        echo "$indexCA $indexCB $indexSG $indexCD" >> chi2.ndx 

        if [ ! -f chi1.xvg ] ; then 
            gmx angle -f ../Production/$MOLEC.xtc \
                -n chi1.ndx \
                -type dihedral \
                -ov chi1.xvg >> $logFile 2>> $errFile  
        fi 
        check chi1.xvg 

        if [ ! -f chi2.xvg ] ; then 
            gmx angle -f ../Production/$MOLEC.xtc \
                -n chi2.ndx \
                -type dihedral \
                -ov chi2.xvg >> $logFile 2>> $errFile  
        fi 
        check chi2.xvg 

        clean
        printf "Success\n"
        cd ../
    else
        printf "Skipped\n"
        fi
}

printf "\n\t\t*** Program Beginning ***\n\n" 
cd $MOLEC
#protein_steep
#solvate
#solvent_steep
#solvent_nvt
#solvent_npt
#production 
nopbc
if grep -sq CNC $MOLEC.pdb ; then 
    hbond 
    hbond_1
    hbond_2
    fit_hbond
#    sasa
#    force_calc
#    chi
    fi 
#minimage
#rmsd 
#rgyr
cd ../

printf "\n\n\t\t*** Program Ending    ***\n\n" 
