#!/bin/bash

SINGLE_SCREEN="true"

SWAP_SCREENS="false"

PERSIST_CACHE="false"

DISABLE_CACHE="false"

# set the name for the monitors (use arandr to get them)
# MONITORS="--presenter-screen=eDP-1 --presentation-screen=DP-1"

# load special pdfpc notes
NOTES=""

# timer in minutes
TIMER="20"

# start from page #n
PAGE=""

#############################################################

ARG_SINGLE_SCREEN=""

if [ "$SINGLE_SCREEN" == "true" ]; then
  ARG_SINGLE_SCREEN="-S"
fi

ARG_SWAP_SCREENS=""

if [ "$SWAP_SCREENS" == "true" ]; then
  ARG_SWAP_SCREENS="-s"
fi

ARG_PERSIST_CACHE=""

if [ "$PERSIST_CACHE" == "true" ]; then
  ARG_PERSIST_CACHE="-s"
fi

ARG_DISABLE_CACHE=""

if [ "$DISABLE_CACHE" == "true" ]; then
  ARG_DISABLE_CACHE="--disable-cache"
fi

ARG_TIMER=""

if [ -n "$TIMER" ]; then
  ARG_TIMER="--duration=$TIMER"
fi

ARG_NOTES=""

if [ -n "$NOTES" ]; then
  ARG_NOTES="--notes=$NOTES"
fi

ARG_PAGE=""

if [ -n "$PAGE" ]; then
  ARG_PAGE="-P $PAGE"
fi

cmd="~/git/linux-setup/submodules/pdfpc/build/bin/pdfpc $MONITORS --note-format=markdown $ARG_PAGE $ARG_NOTES $ARG_SINGLE_SCREEN $ARG_SWAP_SCREENS $ARG_TIMER $ARG_PERSIST_CACHE $ARG_DISABLE_CACHE ./build/main.pdf"
echo $cmd
eval `echo $cmd`
