#!/bin/bash
# 'MakeFemtoDst.sh'

lvl=$1
typ=$2
if [ -z $lvl ]; then
  echo "HEY! Please specify \"particle\" or \"detector\"."
  exit
fi
if [ -z $typ ]; then
  echo "HEY! Please specify \"charged\" or \"full\"."
  exit
fi

# particle case
if [ $lvl = "particle" ]; then
  if [ $typ = "charged" ]; then
    root -b -q MakeFemtoDst.C'(250000,0,-1,true,false,0,"pp200py.effTest.et920pi0par.r03rm1chrg.d17m4y2018.root")'
  elif [ $typ = "full" ]; then
    root -b -q MakeFemtoDst.C'(250000,0,-1,true,false,1,"pp200py.effTest.et920pi0par.r03rm1full.d17m4y2018.root")'
  else
    echo "Hmmm... check what you typed..."
    exit
  fi
# detector case
elif [ $lvl = "detector" ]; then
  if [ $typ = "charged" ]; then
    root -b -q MakeFemtoDst.C'(250000,0,-1,false,true,0,"pp200py.effTest.eTtrg920pi0det.r03rm1chrg.d17m4y2018.root")'
  elif [ $typ = "full" ]; then
    root -b -q MakeFemtoDst.C'(250000,0,-1,false,true,1,"pp200py.effTest.eTtrg920pi0det.r03rm1full.d17m4y2018.root")'
  else
    echo "Hmmm... check what you typed..."
    exit
  fi
# other case
else
  echo "Hmmm... check what you typed..."
  exit
fi
