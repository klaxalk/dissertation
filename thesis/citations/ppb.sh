#!/bin/bash

# get me vim, we will be using it alot to postprocess the generated json files
if [ -x "$(whereis nvim | awk '{print $2}')" ]; then
  VIM_BIN="$(whereis nvim | awk '{print $2}')"
  HEADLESS="--headless"
elif [ -x "$(whereis vim | awk '{print $2}')" ]; then
  VIM_BIN="$(whereis vim | awk '{print $2}')"
  HEADLESS=""
fi

FULL_DIR=full_list
NO_AUTOCIT_DIR=no_autocit
ORIGINAL_DIR=original_exports
AUTHORS_DIR=author_list

FILES=`ls $ORIGINAL_DIR | grep .bib`

for file in $FILES; do

  INPUT_FILE=${file%.bib}

  # copy the original file
  filename=$FULL_DIR/$file
  original=$ORIGINAL_DIR/$file
  cp "$original" "$filename"

  echo Pre-processing $filename

  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^unique-id/norm ct=keywords ^f{C{$INPUT_FILE, mine}," -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^doi/norm maf{%mbd'a" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^isbn/norm maf{%mbd'a" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^issn/norm maf{%mbd'a" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^eissn/norm maf{%mbd'a" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^note/norm maf{%mbd'a" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^editor/norm maf{%mbd'a" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^series/norm maf{%mbd'a" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^organization/norm maf{%mbd'a" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^orcid-numbers/norm maf{%mbd'a" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^researcherid-numbers/norm maf{%mbd'a" -c "wqa" -- "$filename"

  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^pages/norm :s/--/-/g" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^pages/norm :s/-/--/g" -c "wqa" -- "$filename"

  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^booktitle/norm f{v%:s/\<\(\w\)\(\S*\)/\u\1\L\2/g" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^journal/norm f{v%:s/\<\(\w\)\(\S*\)/\u\1\L\2/g" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^title/norm f{v%:s/\<\(\w\)\(\S*\)/\u\1\L\2/g" -c "wqa" -- "$filename"

  $VIM_BIN $HEADLESS -nEs -c "%s/ieee/IEEE/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/rsj/RSJ/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/case/CASE/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/icuas/ICUAS/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/iros/IROS/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/iccps/ICCPS/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/aamas/AAMAS/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/iecon/IECON/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/xxi/XXI/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/icra/ICRA/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/etfa/ETFA/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/(irc/(IRC/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/ecmr/ECMR/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/spie/SPIE/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/ccdc/CCDC/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/ijcnn/IJCNN/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/icarcv/ICARCV/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/(ae)/(AE)/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/mbzirc/MBZIRC/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/mbz /MBZ /gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/Cern@school/CERN School/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/cern/CERN/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/vzlusat/VZLUSAT/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/euv/EUV/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/mesas/MESAS/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/rcar/RCAR/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/wrx-r/WRX-R/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/med-hoc-net/Med-Hoc-Net/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/eth-mav/ETH-MAV/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/uav/UAV/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/mavs/MAVs/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/X-ray/X-Ray/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/\<Uv\>/UV/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/Rtdp/RTDP/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/Hierarchy of Ai/Hierarchy of AI/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/gnss/GNSS/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/gps/GPS/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/embe dded/Embedded/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/wf-iot/WF-IoT/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/cnn/CNN/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/\svr\s/ VR /gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/space vi/Space VI/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/dcad/DCAD/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/nmpc/NMPC/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/\<nips\>/NIPS/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/px4/PX4/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/icsse/ICSSE/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/2-dof/2-DOF/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/ofCFA(2)FBcontrol/of CFA2FB Control/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/a CASE study/A Case Study/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/auv/AUV/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/msalc/MSALC/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/icaaid/ICAAID/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/a-accelerators/Accelerators/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/xvi/XVI/gi" -c "wqa" -- "$filename"

  $VIM_BIN $HEADLESS -nEs -c "%s/\sOf\s/ of /g" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/\sThe\s/ the /g" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/\sOn\s/ on /g" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/\sFor\s/ for /g" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/\sAnd\s/ and /g" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/Risepix-a/Risepix --- A/g" -c "wqa" -- "$filename"

  $VIM_BIN $HEADLESS -nEs -c "%g/\cmonth.*=/norm f=lC 1," -c "wqa" -- "$filename"

  filename2=$NO_AUTOCIT_DIR/$file
  cp "$filename" "$filename2"

  echo Removing second-level autocitations from $filename2

  for author in `cat $AUTHORS_DIR/$INPUT_FILE.txt`; do

    $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^author.*\(\(},\)\@<!\n.*\)\{0,5\}$author/norm dap" -c "wqa" -- "$filename2"

  done

done
