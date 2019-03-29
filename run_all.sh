#!/bin/bash

molecList="
SNase_WT 
SNase_V23X
SNase_L25X
SNase_L38X
SNase_A58X 
SNase_T62X
SNase_V66X
SNase_A90X
SNase_I92X
SNase_A109X
SNase_V104X
SNase_N118X
"

for molec in $molecList ; do 
    printf "\n\t$molec\n"
    if [ ! -f StartingStructures/$molec.pdb ] ; then 
        echo "ERROR: $molec.pdb not found in starting structures. Continuing" 
        continue
        fi 

    if [ ! -f submit_$molec ] ; then 
        echo "#!/bin/bash" > submit_$molec
        echo >> submit_$molec
        echo "#SBATCH -J $molec " >> submit_$molec
        echo "#SBATCH -o $molec.o%j" >> submit_$molec 
        echo "#SBATCH -N 1 "  >> submit_$molec
        echo "#SBATCH -n 48 "  >> submit_$molec
        echo "#SBATCH -p skx-normal" >> submit_$molec 
        echo "#SBATCH -t 48:00:00" >> submit_$molec 
        echo "#SBATCH -A Ras "  >> submit_$molec
        echo "#SBATCH --mail-user=Jeremy_first@utexas.edu"  >> submit_$molec
        echo "#SBATCH --mail-type=all"  >> submit_$molec

        echo >> submit_$molec
        echo "module load gromacs" >> submit_$molec 
        
        echo >> submit_$molec 
        echo "bash run_SNase.sh StartingStructures/$molec.pdb "  >> submit_$molec

        fi 
     
     #sbatch submit_$molec
     bash run_SNase.sh StartingStructures/$molec.pdb
     #bash force_calc_APBS2.sh $molec
     done 

        



