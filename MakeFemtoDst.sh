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
    root -b -q MakeFemtoDst.C'(100000,0,-1,true,false,0,"check.moreEvts.eTtrg920gam.r03rm1chrg.d4m12y2017.root")'
  elif [ $typ = "full" ]; then
    root -b -q MakeFemtoDst.C'(100000,0,-1,true,false,1,"check.moreEvts.eTtrg920gam.r03rm1full.d4m12y2017.root")'
  else
    echo "Hmmm... check what you typed..."
    exit
  fi
# detector case
elif [ $lvl = "detector" ]; then
  if [ $typ = "charged" ]; then
    root -b -q MakeFemtoDst.C'(100000,0,-1,false,true,0,"check.moreEvts.eTtrg920gamEff91.r03rm1chrg.d4m12y2017.root")'
  elif [ $typ = "full" ]; then
    root -b -q MakeFemtoDst.C'(100000,0,-1,false,true,1,"check.moreEvts.eTtrg920gamEff91.r03rm1full.d4m12y2017.root")'
  else
    echo "Hmmm... check what you typed..."
    exit
  fi
# other case
else
  echo "Hmmm... check what you typed..."
  exit
fi
