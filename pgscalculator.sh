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
 echo "-i <dir> 	 path to sumstats folder (or posteriors folder if -1)"
 echo "-l <dir> 	 LD map dir, absolute paths"
 echo "-g <dir> 	 target genotypes"
 echo "-f <file> 	 target genotypes files in genotype folder"
 echo "-s <file> 	 target genotypes snp list filtering(default: none)"
 echo "-c <file> 	 run specific config file"
 echo "-o <dir> 	 path to output directory"
 echo "-b <dir> 	 path to system tmp or scratch (default: /tmp)"
 echo "-w <dir> 	 path to workdir/intermediate files (default: work)"
 echo "-j  	 	 image mode, run docker or singularity (default: singularity)"
 echo "-d  	 	 dev mode, no cleanup of intermediates(default: not active)"
 echo "-v  	 	 get the version number"
 echo "-1  	 	 disable step1, calc posteriors "
 echo "-2  	 	 disable step2, calc score "
 echo ""
 echo "NOTES: 	 By default all steps are run (step1 and step2). "
 echo "          Use -1 or -2 to disable one of the steps. "
 echo "          Disabling all steps will only return formatted sumstats."
 
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

# check for modifiers
if [ ${paramarray[0]} == "placeholder" ] ; then
  runtype="placeholder"
  # remove modifier, 1st element
  paramarray=("${paramarray[@]:1}")
elif [ ${paramarray[0]} == "test" ] ; then
  runtype="test"
  paramarray=("${paramarray[@]:1}")
elif [ ${paramarray[0]} == "utest" ] ; then
  runtype="utest"
  paramarray=("${paramarray[@]:1}")
elif [ ${paramarray[0]} == "etest" ] ; then
  runtype="etest"
  paramarray=("${paramarray[@]:1}")
else
  runtype="default"
fi

# starting getops with :, puts the checking in silent mode for errors.
getoptsstring=":hvi:o:b:w:l:g:f:s:m:c:db:12j:"

infold=""
lddir=""
genodir=""
genofile=""
snpfile=""
conffile=""
outdir="out"
container_image=""
calc_posterior=true
calc_score=true

# some logical defaults
infold_given=false
lddir_given=false
genodir_given=false
genofile_given=false
snpfile_given=false
conffile_given=false
outdir_given=false
tmpdir_given=false
devmode_given=false
calc_posterior_given=false
calc_score_given=false
container_image_given=false

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
    i )
      infold="$OPTARG"
      infold_given=true
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
    s )
      snpfile="$OPTARG"
      snpfile_given=true
      ;;
    c )
      conffile="$OPTARG"
      conffile_given=true
      ;;
    o )
      outdir="$OPTARG"
      outdir_given=true
      ;;
    j )
      container_image="$OPTARG"
      container_image_given=true
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
    1 )
      calc_posterior=false
      calc_posterior_given=true
      ;;
    2 )
      calc_score=false
      calc_score_given=true
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
# Setup quick-run example options
################################################################################
if [ "${runtype}" == "test" ] || [ "${runtype}" == "utest" ] || [ "${runtype}" == "etest" ]; then
  # All are placeholders and not used
  infold="${project_dir}/tests"
  genodir="${project_dir}/tests"
  genofile="${project_dir}/VERSION"
  conffile="${project_dir}/VERSION"
  lddir="${project_dir}/tests"
fi

################################################################################
# Check if the provided paths exist
################################################################################

# make outdir if it doesn't already exist
mkdir -p ${outdir}

# make workdir if it doesn't already exist
mkdir -p ${workdir}

# make workdir if it doesn't already exist
mkdir -p ${tmpdir}

infold_host=$(realpath "${infold}")
conffile_host=$(realpath "${conffile}")
outdir_host=$(realpath "${outdir}")
tmpdir_host=$(realpath "${tmpdir}")
workdir_host=$(realpath "${workdir}")

# Test that file and folder exists, all of these will always get mounted
if [ ! -d $infold_host ]; then
  >&2 echo "infold doesn't exist"
  >&2 echo "path tried: $infold_host"
  exit 1
fi

# if calc posterior active
if $calc_posterior; then
  lddir_host=$(realpath "${lddir}")
  if [ ! -d $lddir_host ]; then
    >&2 echo "lddir doesn't exist"
    >&2 echo "path tried: $lddir_host"
    exit 1
  fi
else
  lddir_host="$(realpath ${project_dir}/tests/example_data/ldref)"
fi

# if calc score active
#if $calc_score; then
#  if [ "$build" != "37" ] && [ "$build" != "38" ]; then
#    >&2 echo "build not available"
#    >&2 echo "build tried: $build"
#    exit 1
#  fi
#else
#  build="not_defined"
#fi

#if $calc_score; then

# Always use genodir and genofile. Because we need at least the .bim for the posterior calculation
  genodir_host=$(realpath "${genodir}")
  if [ ! -d $genodir_host ]; then
    >&2 echo "genodir doesn't exist"
    >&2 echo "path tried: $genodir_host"
    exit 1
  fi

  genofile_host=$(realpath "${genofile}")
  if [ ! -f $genofile_host ]; then
    >&2 echo "genofile doesn't exist"
    >&2 echo "path tried: $genofile_host"
    exit 1
  fi

  if ${snpfile_given}; then
    snpfile_host=$(realpath "${snpfile}")
    if [ ! -f $snpfile_host ]; then
      >&2 echo "snpfile doesn't exist"
      >&2 echo "path tried: $snpfile_host"
      exit 1
    fi
  fi

#else
#  genodir_host="$(realpath ${project_dir}/tests/example_data/genotypes)"
#  genofile_host="$(realpath ${project_dir}/tests/example_data/genotypes/genofiles_placehoder.txt)"
#fi

# Always check these
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

# which mount symbol to use
if [ "${container_image}" == "docker" ]; then
  mountflag="-v"
elif [ "${container_image}" == "dockerhub_biopsyk" ]; then
  mountflag="-v"
else
  mountflag="-B"
fi

# indir
#indir_host=$(dirname "${infile_host}")
#infile_name=$(basename "${infile_host}")
indir_container="/pgscalculator/input"
#infile_container="${indir_container}/${infile_name}"

# lddir
foldername=$(basename "$lddir_host")
lddir_container="/pgscalculator/$foldername"
#lddir_container="/pgscalculator/lddir"

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

if ${snpfile_given}; then
  snpdir_host=$(dirname "${snpfile_host}")
  snpfile_name=$(basename "${snpfile_host}")
  snpdir_container="/pgscalculator/snpdir"
  snpfile_container="${snpdir_container}/${snpfile_name}"
  snplist_host_container="${mountflag} ${snpdir_host}:${snpdir_container}"
  snplist_container="--snplist ${snpfile_container}"
else
  snplist_host_container=""
  snplist_container=""
fi

# outdir
outdir_container="/pgscalculator/outdir"

# tmpdir
tmpdir_container="/tmp"

# workdir
workdir_container="/pgscalculator/work"

# Use outdir as fake home to avoid lock issues for the hidden .nextflow/history file
FAKE_HOME="${outdir_container}"
export SINGULARITY_HOME="${FAKE_HOME}"
export APPTAINER_HOME="${FAKE_HOME}"


# Which runscript to use
if [ "${runtype}" == "default" ]; then
  run_script="/pgscalculator/main.nf"
elif [ "${runtype}" == "test" ]; then
  run_script="/pgscalculator/tests/run-tests.sh"
elif [ "${runtype}" == "utest" ]; then
  run_script="/pgscalculator/tests/run-unit-tests.sh"
elif [ "${runtype}" == "etest" ]; then
  mkdir -p tmp
  run_script="/pgscalculator/tests/run-e2e-tests.sh"
else
  echo "option not available"
  exit 1
fi

# set image tags and make mount function available
source "${project_dir}/scripts/init-containerization.sh"

# Which image is to be used
if [ "${container_image}" == "docker" ]; then
  runimage="${image_tag}" 
elif [ "${container_image}" == "dockerhub_biopsyk" ]; then
  runimage="${deploy_image_tag_docker_hub}" 
elif [ "${container_image}" == "" ]; then
  #if not set, assume image is in sif folder
  runimage="sif/${singularity_image_tag}" 
else
  runimage="${container_image}" 
fi

# Source Docker configuration from the conf folder
source "${project_dir}/conf/init-docker-config.sh"

if [ "${runtype}" == "test" ] || [ "${runtype}" == "utest" ] || [ "${runtype}" == "etest" ]; then
  if [ "${container_image}" == "dockerhub_biopsyk" ]; then
    echo "container: $runimage"
    mount_flags=$(format_mount_flags "${mountflag}")
    #exec docker run --rm "${deploy_image_tag_docker_hub}" "${run_script}"
    exec docker run --rm ${mount_flags} "${runimage}" ${run_script}
  elif [ "${container_image}" == "docker" ]; then
    echo "container: $runimage"
    mount_flags=$(format_mount_flags "${mountflag}")
    exec docker run --rm ${mount_flags} "${runimage}" ${run_script}
  else
    echo "container: $runimage"
    mount_flags=$(format_mount_flags "${mountflag}")
    singularity run --contain --cleanenv ${mount_flags} "${runimage}" ${run_script}
  fi
elif [ "${container_image}" == "docker" ] || [ "${container_image}" == "dockerhub_biopsyk" ]; then
  echo "container: $runimage"
  mount_flags=$(format_mount_flags "${mountflag}")
  exec docker run \
     --rm \
     ${docker_run_args} \
     ${mount_flags} \
     ${mountflag} "${infold_host}:${indir_container}" \
     ${mountflag} "${outdir_host}:${outdir_container}" \
     ${mountflag} "${lddir_host}:${lddir_container}" \
     ${mountflag} "${genodir_host}:${genodir_container}" \
     ${mountflag} "${genodir2_host}:${genodir2_container}" \
     ${mountflag} "${confdir_host}:${confdir_container}" \
     ${mountflag} "${tmpdir_host}:${tmpdir_container}" \
     ${mountflag} "${workdir_host}:${workdir_container}" \
     ${snplist_host_container} \
     "${runimage}" \
     nextflow \
       -log "${outdir_container}/.nextflow.log" \
       -c ${conffile_container} \
       run ${run_script} \
       ${devmode} \
       --calc_posterior ${calc_posterior} \
       --calc_score ${calc_score} \
       --input "${indir_container}" \
       --lddir "${lddir_container}" \
       --genodir "${genodir_container}" \
       --genofile "${genofile_container}" \
       ${snplist_container} \
       --conffile "${conffile_container}" \
       --outdir "${outdir_container}"
else
  echo "container: $runimage"
  mount_flags=$(format_mount_flags "${mountflag}")
  singularity run \
     --contain \
     --cleanenv \
     ${mount_flags} \
     ${mountflag} "${infold_host}:${indir_container}" \
     ${mountflag} "${outdir_host}:${outdir_container}" \
     ${mountflag} "${lddir_host}:${lddir_container}" \
     ${mountflag} "${genodir_host}:${genodir_container}" \
     ${mountflag} "${genodir2_host}:${genodir2_container}" \
     ${mountflag} "${confdir_host}:${confdir_container}" \
     ${mountflag} "${tmpdir_host}:${tmpdir_container}" \
     ${mountflag} "${workdir_host}:${workdir_container}" \
     ${snplist_host_container} \
     "${runimage}" \
     nextflow \
       -log "${outdir_container}/.nextflow.log" \
       -c ${conffile_container} \
       run ${run_script} \
       ${devmode} \
       --calc_posterior ${calc_posterior} \
       --calc_score ${calc_score} \
       --input "${indir_container}" \
       --lddir "${lddir_container}" \
       --genodir "${genodir_container}" \
       --genofile "${genofile_container}" \
       ${snplist_container} \
       --conffile "${conffile_container}" \
       --outdir "${outdir_container}" 
fi       



#remove .nextflow directory by default
if ${devmode_given} ;
then
  #Set correct permissions to pipeline_info files
  chmod -R ugo+rwX ${outdir_host}/pipeline_info
elif [ "${runtype}" == "test" ] || [ "${runtype}" == "utest" ] || [ "${runtype}" == "etest" ]; then
  :
else
  #Set correct permissions to pipeline_info files
  chmod -R ugo+rwX ${outdir_host}/pipeline_info

  function cleanup {
    echo ">> Cleaning up (disable with -d) "
    echo ">> Removing ${outdir_host}/.nextflow"
    rm -r ${outdir_host}/.nextflow
    echo ">> Done"
  }
  trap cleanup EXIT
fi
 
echo "pgscalculator.sh reached the end: $(date)"
