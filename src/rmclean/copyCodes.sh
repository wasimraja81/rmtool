#!/bin/bash -l 
SOURCE_PATH="/home/wasim/softwares/CURR_DEVEL/RM_CLEAN_TESTS/SOURCE/SUBROUTINES"
DEST_PATH="/home/wasim/softwares/rmtool/src/rmclean"

codeList=("rm_clean.f" \
	"compute_dirty_rmbeam.f" \
	"rm_restore.f" \
	"peak_interp.f" \
	"quad_interp.f" \
	"fourier_interp.f" \
	"fourier_interp_re.f" \
	"sinc_interp.f" \
	"index_absmax.f" \
	"bw_depol_correct.f" \
	"bw_depol_correct_setup.f")

nCodes=${#codeList[@]}
echo "Number of codes to copy: ${nCodes}"

for (( iCode=0; iCode<${nCodes}; iCode++ ))
do
	codeNow="${codeList[${iCode}]}"
	convertedCode="${codeList[${iCode}]%%.f}.f90"
	echo "Converting code: ${SOURCE_PATH}/${codeNow} to ${DEST_PATH}/${convertedCode}"
	findent < "${SOURCE_PATH}/${codeNow}" > "${DEST_PATH}/${convertedCode}"
done
