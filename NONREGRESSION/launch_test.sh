#!/bin/bash
# Launch the test and compares obtained results to expected ones
set -euo pipefail 
#set -x

# This function launch the computation by calling the executable with arguments
# taken from args.txt file
function launch_computation {
  local readonly exe_path=$1
  local readonly data_dir=$2
  local readonly mpi_launcher=$3
  local return_code=0
  ${mpi_launcher} -n 1 $1 $data_dir/Donnees.arc
  if [[ $? -ne 0 ]]; then
    echo "A problem occured during test execution."
    return_code=1
  fi
  return ${return_code}
}
# This function launch the computation by calling the executable with arguments
# taken from args.txt file
function launch_computation_seq_pr {
  local readonly exe_path=$1
  local readonly data_dir=$2
  local readonly mpi_launcher=$3
  local return_code=0
  ${mpi_launcher} -n 1 $1 -arcane_opt max_iteration 10 $data_dir/Donnees.arc
  ${mpi_launcher} -n 1 $1 -arcane_opt continue $data_dir/Donnees.arc
  if [[ $? -ne 0 ]]; then
    echo "A problem occured during test execution."
    return_code=1
  fi
  return ${return_code}
}
# This function launch the computation by calling the executable with arguments
# taken from args.txt file
function launch_computation_para_4 {
  local readonly exe_path=$1
  local readonly data_dir=$2
  local readonly mpi_launcher=$3
  local return_code=0
  ${mpi_launcher} -n 4 $1 $data_dir/Donnees.arc
  if [[ $? -ne 0 ]]; then
    echo "A problem occured during test execution."
    return_code=1
  fi
  return ${return_code}
}
# This function launch the computation by calling the executable with arguments
# taken from args.txt file
function launch_computation_para_8 {
  local readonly exe_path=$1
  local readonly data_dir=$2
  local readonly mpi_launcher=$3
  local return_code=0
  ${mpi_launcher} -n 8 $1 $data_dir/Donnees.arc
  if [[ $? -ne 0 ]]; then
    echo "A problem occured during test execution."
    return_code=1
  fi
  return ${return_code}
}
# This function launch the computation by calling the executable with arguments
# taken from args.txt file
function launch_computation_cuda_1 {
  local readonly exe_path=$1
  local readonly data_dir=$2
  local readonly mpi_launcher=$3
  local return_code=0
  ${mpi_launcher} -n 1 $1 -A,AcceleratorRuntime=cuda $data_dir/Donnees.arc
  if [[ $? -ne 0 ]]; then
    echo "A problem occured during test execution."
    return_code=1
  fi
  return ${return_code}
}
# This function compares the results obtained in the output directory with those
# expected in the reference directory
function compare_results {
  local readonly reference_dir=$1
  local return_code=0
  ls -l "$reference_dir/output"
  echo "Differences in : ${PWD}/DIFF.txt"
  diff -r output/depouillement "$reference_dir/output/depouillement" > DIFF.txt 2>&1
  if [[ $? -ne 0 ]]; then
      echo "Test failure! Output is different from reference"
      echo ${test_dir} >>  list_of_cases_to_change
    return_code=1
  fi
  return ${return_code}
}

# Main function. Calls launch_computation and compare_results
# For each test a directory is created inside /tmp
function main {
  local readonly exe_path=$1
  local readonly test_dir=$2
  local readonly type=$3
  local readonly mpi_launcher=$4
  local readonly basetmpdir_nr=$5

  local readonly test_name=$(basename ${test_dir})
  echo "Executable path : ${exe_path}"
  echo "Executing test ${test_name}"
  echo "Executing test mode ${type}"
  if [ ${type} = "para" ]; then
    echo "Executing test parallele mode ${para}"
  fi
  
  tmpdir_base="${basetmpdir_nr}/MAHYCO/NONREGRESSION"
  mkdir -p ${tmpdir_base}
  if [ ${?} != 0 ]
  then
    echo "Unable to create ${$tmpdir_base}"
    echo "Aborting!"
    exit 4
  fi
  tmp_dir=$(mktemp -d -p ${tmpdir_base} mahyco-ci-${type}-${test_name}-XXXXXXXXXX)

  if [[ ! -d ${tmp_dir} ]]; then
    echo "Unable to create a temporary directory!"
    echo "Aborting!"
    exit 3
  fi

  cd ${tmp_dir}
  echo "This directory contains the output of the test under ${test_dir}" > README.txt
  pwd
  echo ${type}
  if [ ${type}  = "para_8" ]
  then
      echo " lancement parallele sur 8 coeurs" 
      launch_computation_para_8 ${exe_path} ${test_dir} ${mpi_launcher}
  elif [ ${type}  = "para_4" ]
  then
      echo " lancement parallele sur 4 coeurs" 
      launch_computation_para_4 ${exe_path} ${test_dir} ${mpi_launcher}
  elif [ ${type}  = "seq_pr" ]
  then
      echo " lancement sequentiel protection-reprise" 
      launch_computation_seq_pr ${exe_path} ${test_dir} ${mpi_launcher}
  elif [ ${type}  = "cuda_1" ]
  then
      echo " lancement sequentiel 1 GPU" 
      launch_computation_cuda_1 ${exe_path} ${test_dir} ${mpi_launcher}
  else
      echo " lancement sequentiel" 
      launch_computation ${exe_path} ${test_dir} ${mpi_launcher}
  fi    
  if [[ $? -ne 0 ]]; then
    echo "Aborting!"
    exit 1
  fi

  compare_results ${test_dir}
  if [[ $? -ne 0 ]]; then
    echo "Aborting!"
    exit 2
  else
    echo "Success!"
  fi
  exit 0
}

main $@
