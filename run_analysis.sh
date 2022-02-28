#!/bin/bash

SHORT=i:,T:,r:,s:,h

OPTS=$(getopt --options $SHORT -- "$@")

PANEL="Untargeted"
RES=1.4
STEMNESS="Stem"
while getopts $SHORT opt; do 
    case ${opt} in
        i )
            ID=$OPTARG
            ;;
        T )
            PANEL=$OPTARG
            ;;
        r )
            RES=$OPTARG
            ;;
        s )
            STEMNESS=$OPTARG
            ;;
        \? )
            echo "Invalid option: $OPTARG" 1>&2
            ;;
        : )
            echo "Invalid option: $OPTARG requires an argument" 1>&2
            ;;
    esac
done 

( export ID=$ID PANEL=$PANEL RES=$RES STEMNESS=$STEMNESS; sbatch slurm/start_subset.sh )
