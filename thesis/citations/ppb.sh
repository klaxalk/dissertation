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
AUTHORS_DIR=author_list

FILES=`ls $FULL_DIR | grep .bib`

for file in $FILES; do

  INPUT_FILE=${file%.bib}

  filename=$FULL_DIR/$file

  echo Pre-processing $filename

  $VIM_BIN $HEADLESS -nEs -c "%g/Booktitle/norm f{v%:s/\<\(\w\)\(\S*\)/\u\1\L\2/g" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%g/Journal/norm f{v%:s/\<\(\w\)\(\S*\)/\u\1\L\2/g" -c "wqa" -- "$filename"

  $VIM_BIN $HEADLESS -nEs -c "%g/^Unique-ID/norm ct=keywords ^f{C{$INPUT_FILE}," -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^doi/norm dd" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^isbn/norm dd" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^issn/norm dd" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^note/norm dd" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^editor/norm dd" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^series/norm dd" -c "wqa" -- "$filename"

  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^pages/norm :s/--/-/g" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^pages/norm :s/-/--/g" -c "wqa" -- "$filename"

  $VIM_BIN $HEADLESS -nEs -c "%s/ieee/IEEE/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/rsj/RSJ/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/case/CASE/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/icuas/ICUAS/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/iros/IROS/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/icra/ICRA/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/etfa/ETFA/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/(irc/(IRC/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/ecmr/ECMR/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/spie/SPIE/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/ijcnn/IJCNN/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/icarcv/ICARCV/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/(ae)/(AE)/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/mbzirc/MBZIRC/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/Cern@school/CERN school/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/cern/CERN/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/vzlusat/VZLUSAT/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/euv/EUV/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/mesas/MESAS/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/rcar/RCAR/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/wrx-r/WRX-R/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/med-hoc-net/Med-Hoc-Net/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/eth-mav/ETH-MAV/gi" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/uav/UAV/gi" -c "wqa" -- "$filename"

  $VIM_BIN $HEADLESS -nEs -c "%s/ Of / of /g" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/ The / the /g" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/ On / on /g" -c "wqa" -- "$filename"
  $VIM_BIN $HEADLESS -nEs -c "%s/ For / for /g" -c "wqa" -- "$filename"

  filename2=$NO_AUTOCIT_DIR/$file
  cp "$filename" "$filename2"

  echo Removing second-level autocitations from $filename2

  for author in `cat $AUTHORS_DIR/$INPUT_FILE.txt`; do

    $VIM_BIN $HEADLESS -nEs -c "set ignorecase" -c "%g/^author.*$author/norm dap" -c "wqa" -- "$filename2"

  done

done
