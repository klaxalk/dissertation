#!/bin/bash

# get me vim, we will be using it alot to postprocess the generated json files
if [ -x "$(whereis nvim | awk '{print $2}')" ]; then
  VIM_BIN="$(whereis nvim | awk '{print $2}')"
  HEADLESS="--headless"
elif [ -x "$(whereis vim | awk '{print $2}')" ]; then
  VIM_BIN="$(whereis vim | awk '{print $2}')"
  HEADLESS=""
fi

INPUT_FILE=papis_export.bib
OUTPUT_FILE=main.bib

cp $INPUT_FILE $OUTPUT_FILE

$VIM_BIN $HEADLESS -nEs -c '%g/url =/norm dd' -c "wqa" -- "$OUTPUT_FILE"
$VIM_BIN $HEADLESS -nEs -c '%g/doi =/norm dd' -c "wqa" -- "$OUTPUT_FILE"
$VIM_BIN $HEADLESS -nEs -c '%s/Baca, Tomas/\\textbf\{T. Baca\}/g' -c "wqa" -- "$OUTPUT_FILE"
$VIM_BIN $HEADLESS -nEs -c '%s/Baca, T\./\\textbf\{T. Baca\}/g' -c "wqa" -- "$OUTPUT_FILE"
