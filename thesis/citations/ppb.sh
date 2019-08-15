#!/bin/bash

# get me vim, we will be using it alot to postprocess the generated json files
if [ -x "$(whereis nvim | awk '{print $2}')" ]; then
  VIM_BIN="$(whereis nvim | awk '{print $2}')"
  HEADLESS="--headless"
elif [ -x "$(whereis vim | awk '{print $2}')" ]; then
  VIM_BIN="$(whereis vim | awk '{print $2}')"
  HEADLESS=""
fi

FILES=`ls | grep .bib`

for file in $FILES; do

  INPUT_FILE=${file%.bib}

  echo $INPUT_FILE

  $VIM_BIN $HEADLESS -nEs -c "%g/Booktitle/norm f{v%:s/\<\(\w\)\(\S*\)/\u\1\L\2/g" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "%g/Journal/norm f{v%:s/\<\(\w\)\(\S*\)/\u\1\L\2/g" -c "wqa" -- "$file"

  $VIM_BIN $HEADLESS -nEs -c "%g/^Unique-ID/norm ct=keywords ^f{C{$INPUT_FILE}," -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^doi/norm dd" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^isbn/norm dd" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^issn/norm dd" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^note/norm dd" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^editor/norm dd" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^series/norm dd" -c "wqa" -- "$file"

  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^pages/norm :s/--/-/g" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^pages/norm :s/-/--/g" -c "wqa" -- "$file"

  $VIM_BIN $HEADLESS -nEs -c "%s/ieee/IEEE/gi" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "%s/rsj/RSJ/gi" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "%s/case/CASE/gi" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "%s/icuas/ICUAS/gi" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "%s/iros/IROS/gi" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "%s/icra/ICRA/gi" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "%s/etfa/ETFA/gi" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "%s/(irc/(IRC/gi" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "%s/ecmr/ECMR/gi" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "%s/spie/SPIE/gi" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "%s/ijcnn/IJCNN/gi" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "%s/icarcv/ICARCV/gi" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "%s/(ae)/(AE)/gi" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "%s/mbzirc/MBZIRC/gi" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "%s/Cern@school/CERN school/gi" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "%s/cern/CERN/gi" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "%s/vzlusat/VZLUSAT/gi" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "%s/euv/EUV/gi" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "%s/mesas/MESAS/gi" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "%s/rcar/RCAR/gi" -c "wqa" -- "$file"
  $VIM_BIN $HEADLESS -nEs -c "%s/wrx-r/WRX-R/gi" -c "wqa" -- "$file"

done
