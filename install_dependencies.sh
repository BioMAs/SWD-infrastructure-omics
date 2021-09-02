#!/bin/bash

# Stop on error
set -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

MAIN_ENV_NAME=srp

ENVS=$(conda env list | awk '{print $1}' )

FOUND=1

for ENV in ${ENVS}
do
	if [ "${ENV}" == "${MAIN_ENV_NAME}" ]; then
		FOUND=0
	fi
done

# Creation of main conda environment.
if [ ${FOUND} -eq 0 ]; then
	echo "${MAIN_ENV_NAME} already created"
else
	echo "Creating env ${MAIN_ENV_NAME}"
	conda env create -n ${MAIN_ENV_NAME} -f ${DIR}/CONDA/srp.yml
	source activate ${MAIN_ENV_NAME}
	Rscript ${DIR}/CONDA/installDeEnv.R
	source deactivate
fi