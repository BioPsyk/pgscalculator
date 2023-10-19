#!/usr/bin/env bash

# This script is a wrapper to run the singularity image of this pipeline, where folders are parsed and mounted in the right place.

################################################################################
# Help page
################################################################################

function general_usage(){
 echo "Usage:"
 echo " ./pgscalculator.sh -i <file> -o <dir>"
 echo ""
 echo ""
 echo "options:"
 echo "-h		 Display help message for pgscalculator"
 echo "-s <file> 	 path to infile"
 echo "-l <dir> 	 LD map dir, absolute paths"
 echo "-g <dir> 	 target genotypes"
 echo "-f <file> 	 target genotypes files in genotype folder"
 echo "-m <value> 	 method (default: prscs)"
 echo "-c <file> 	 run specific config file"
 echo "-o <dir> 	 path to output directory"
 echo "-b <dir> 	 path to system tmp or scratch (default: /tmp)"
 echo "-w <dir> 	 path to workdir/intermediate files (default: work)"
 echo "-d  	 	 dev mode, saving intermediate files, no cleanup of workdir(default: not active)"
 echo "-v  	 	 get the version number"
}

################################################################################
# Prepare path parsing
################################################################################
# All paths we see will start from the project root, even if the command is called from somewhere else
present_dir="${PWD}"
project_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

################################################################################
# Parameter parsing
################################################################################
# whatever the input make it array
paramarray=($@)

# starting getops with :, puts the checking in silent mode for errors.
getoptsstring=":hvs:o:b:w:l:g:f:m:c:d"

infile=""
ldfile=""
genodir=""
genofile=""
conffile=""
outdir="out"
method=""

# some logical defaults
infile_given=false
lddir_given=false
genodir_given=false
genofile_given=false
conffile_given=false
outdir_given=false
tmpdir_given=false
devmode_given=false

# default system tmp
tmpdir="/tmp"
workdir="${present_dir}/work"
devmode=""

while getopts "${getoptsstring}" opt "${paramarray[@]}"; do
  case ${opt} in
    h )
      general_usage 1>&2
      exit 0
      ;;
    v )
      #write a something that parses the actual version number
      cat ${project_dir}/VERSION 1>&2
      exit 0
      ;;
    s )
      infile="$OPTARG"
      infile_given=true
      ;;
    l )
      lddir="$OPTARG"
      lddir_given=true
      ;;
    g )
      genodir="$OPTARG"
      genodir_given=true
      ;;
    f )
      genofile="$OPTARG"
      genofile_given=true
      ;;
    c )
      conffile="$OPTARG"
      conffile_given=true
      ;;
    m )
      method="$OPTARG"
      method_given=true
      ;;
    o )
      outdir="$OPTARG"
      outdir_given=true
      ;;
    b )
      tmpdir="$OPTARG"
      tmpdir_given=true
      ;;
    w )
      workdir="$OPTARG"
      workdir_given=true
      ;;
    d )
      devmode="--dev"
      devmode_given=true
      ;;
    \? )
      echo "Invalid Option: -$OPTARG" 1>&2
      exit 1
      ;;
    : )
      echo "Invalid Option: -$OPTARG requires an argument" 1>&2
      exit 1
      ;;
  esac
done

################################################################################
# Check if the provided paths exist
################################################################################

# make outdir if it doesn't already exist
mkdir -p ${outdir}

# make workdir if it doesn't already exist
mkdir -p ${workdir}

# make workdir if it doesn't already exist
mkdir -p ${tmpdir}

infile_host=$(realpath "${infile}")
lddir_host=$(realpath "${lddir}")
genodir_host=$(realpath "${genodir}")
genofile_host=$(realpath "${genofile}")
conffile_host=$(realpath "${conffile}")
outdir_host=$(realpath "${outdir}")
tmpdir_host=$(realpath "${tmpdir}")
workdir_host=$(realpath "${workdir}")

# Test that file and folder exists, all of these will always get mounted
if [ ! -f $infile_host ]; then
  >&2 echo "infile doesn't exist"
  exit 1
fi
if [ ! -d $lddir_host ]; then
  >&2 echo "lddir doesn't exist"
  exit 1
fi
if [ ! -d $genodir_host ]; then
  >&2 echo "genodir doesn't exist"
  exit 1
fi
if [ ! -f $genofile_host ]; then
  >&2 echo "genofile doesn't exist"
  >&2 echo "path tried: $genofile_host"
  exit 1
fi
if [ ! -f $conffile_host ]; then
  >&2 echo "conffile doesn't exist"
  >&2 echo "path tried: $conffile_host"
  exit 1
fi
if [ ! -d $outdir_host ]; then
  >&2 echo "outdir doesn't exist"
  exit 1
fi
if [ ! -d $tmpdir_host ]; then
  >&2 echo "tmpdir doesn't exist"
  exit 1
fi
if [ ! -d $workdir_host ]; then
  >&2 echo "workdir doesn't exist"
  exit 1
fi

################################################################################
# Prepare container variables
################################################################################

# All paths we see will start from the project root, even if the command is called from somewhere else
project_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source "${project_dir}/scripts/init-containerization.sh"

mount_flags=$(format_mount_flags "-B")

# indir
indir_host=$(dirname "${infile_host}")
infile_name=$(basename "${infile_host}")
indir_container="/pgscalculator/input"
infile_container="${indir_container}/${infile_name}"

# lddir
lddir_container="/pgscalculator/lddir"

# genodir
genodir_container="/pgscalculator/genodir"

# genofile
genodir2_host=$(dirname "${genofile_host}")
genofile_name=$(basename "${genofile_host}")
genodir2_container="/pgscalculator/genodir2"
genofile_container="${genodir2_container}/${genofile_name}"

# genofile
confdir_host=$(dirname "${conffile_host}")
conffile_name=$(basename "${conffile_host}")
confdir_container="/pgscalculator/confdir"
conffile_container="${confdir_container}/${conffile_name}"

# outdir
outdir_container="/pgscalculator/outdir"

# tmpdir
tmpdir_container="/tmp"

# workdir
workdir_container="/pgscalculator/work"

# Use outdir as fake home to avoid lock issues for the hidden .nextflow/history file
FAKE_HOME="${outdir_container}"
export SINGULARITY_HOME="${FAKE_HOME}"

singularity run \
   --contain \
   --cleanenv \
   ${mount_flags} \
   -B "${indir_host}:${indir_container}" \
   -B "${outdir_host}:${outdir_container}" \
   -B "${lddir_host}:${lddir_container}" \
   -B "${genodir_host}:${genodir_container}" \
   -B "${genodir2_host}:${genodir2_container}" \
   -B "${confdir_host}:${confdir_container}" \
   -B "${tmpdir_host}:${tmpdir_container}" \
   -B "${workdir_host}:${workdir_container}" \
   "tmp/${singularity_image_tag}" \
   nextflow \
     -log "${outdir_container}/.nextflow.log" \
     run /pgscalculator ${runtype} \
     ${devmode} \
     --method ${method} \
     --input "${infile_container}" \
     --lddir "${lddir_container}" \
     --genodir "${genodir_container}" \
     --genofile "${genofile_container}" \
     --conffile "${conffile_container}" \
     --outdir "${outdir_container}" 

#Set correct permissions to pipeline_info files
chmod -R ugo+rwX ${outdir_host}/pipeline_info

#remove .nextflow directory by default
if ${devmode_given} ;
then
  :
else
  function cleanup {
    echo ">> Cleaning up (disable with -l) "
    echo ">> Removing ${outdir_host}/.nextflow"
    rm -r ${outdir_host}/.nextflow
    echo ">> Done"
  }
  trap cleanup EXIT
fi
 
echo "pgscalculator.sh reached the end: $(date)"
