#!/bin/bash

usage='Usage: makeIasPredictionPlotsHtml.sh -d plotDirName -r rootFile'
  
#export SCRAM_ARCH=slc5_amd64_gcc434
#source /afs/cern.ch/cms/sw/cmsset_default.sh
# UMN
source /local/cms/sw/cmsset_default.sh


eval `scramv1 runtime -sh`

args=`getopt dr: -- "$@"`
if test $? != 0
     then
         echo $usage
         exit 1
fi

eval set -- "$args"
for i
  do
  case "$i" in
      -d) shift; plots_dir=$2;shift;;
      -r) shift; root_file=$2;shift;;
  esac
done

if [ "X"${plots_dir} == "X" ]; then
  echo "invalid plots dir! Please give a valid one."
  echo $usage
  exit 1
fi

if [ "X"${root_file} == "X" ]; then
  echo "invalid root file! Please give a valid one."
  echo $usage
  exit 1
fi

my_pwd=`pwd`
mkdir ${plots_dir}
cd ${plots_dir}

root -b $my_pwd/${root_file} <<!

.x /local/cms/user/cooper/cmssw/428/HSCPAnalysis/src/HSCP2011/ToyMC/bin/makeIasPredictionPngsFixedBins.C

.q
!

cat > index.html <<EOF

<HTML>

<HEAD><TITLE>Ias prediction plots</TITLE></HEAD>
 
<BODY link="Red">
<FONT color="Black">
 
<Center>
<h2><FONT color="Blue">Ias points and fit plots</FONT></h2>
</Center>

<table>
<tr>
<td></td>
<td>NoM 5-8</td>
<td>NoM 9-10</td>
<td>NoM 11-12</td>
<td>NoM 13-14</td>
<td>NoM 15-16</td>
<td>NoM 17+</td>
</tr>
<tr>
<td>0.0 <= |eta| < 0.2</td>
<td><a href="nom5to8eta0to2iasPred.png"><img height="300" src="nom5to8eta0to2iasPred.png"></a></td>
<td><a href="nom9to10eta0to2iasPred.png"><img height="300" src="nom9to10eta0to2iasPred.png"></a></td>
<td><a href="nom11to12eta0to2iasPred.png"><img height="300" src="nom11to12eta0to2iasPred.png"></a></td>
<td><a href="nom13to14eta0to2iasPred.png"><img height="300" src="nom13to14eta0to2iasPred.png"></a></td>
<td><a href="nom15to16eta0to2iasPred.png"><img height="300" src="nom15to16eta0to2iasPred.png"></a></td>
<td><a href="nom17to18eta0to2iasPred.png"><img height="300" src="nom17to18eta0to2iasPred.png"></a></td>
<td>0.0 <= |eta| < 0.2</td>
</tr>
<tr>
<td>0.2 < |eta| < 0.4</td>
<td><a href="nom5to8eta2to4iasPred.png"><img height="300" src="nom5to8eta2to4iasPred.png"></a></td>
<td><a href="nom9to10eta2to4iasPred.png"><img height="300" src="nom9to10eta2to4iasPred.png"></a></td>
<td><a href="nom11to12eta2to4iasPred.png"><img height="300" src="nom11to12eta2to4iasPred.png"></a></td>
<td><a href="nom13to14eta2to4iasPred.png"><img height="300" src="nom13to14eta2to4iasPred.png"></a></td>
<td><a href="nom15to16eta2to4iasPred.png"><img height="300" src="nom15to16eta2to4iasPred.png"></a></td>
<td><a href="nom17to18eta2to4iasPred.png"><img height="300" src="nom17to18eta2to4iasPred.png"></a></td>
<td>0.2 < |eta| < 0.4</td>
</tr>
<tr>
<td>0.4 < |eta| < 0.6</td>
<td><a href="nom5to8eta4to6iasPred.png"><img height="300" src="nom5to8eta4to6iasPred.png"></a></td>
<td><a href="nom9to10eta4to6iasPred.png"><img height="300" src="nom9to10eta4to6iasPred.png"></a></td>
<td><a href="nom11to12eta4to6iasPred.png"><img height="300" src="nom11to12eta4to6iasPred.png"></a></td>
<td><a href="nom13to14eta4to6iasPred.png"><img height="300" src="nom13to14eta4to6iasPred.png"></a></td>
<td><a href="nom15to16eta4to6iasPred.png"><img height="300" src="nom15to16eta4to6iasPred.png"></a></td>
<td><a href="nom17to18eta4to6iasPred.png"><img height="300" src="nom17to18eta4to6iasPred.png"></a></td>
<td>0.4 < |eta| < 0.6</td>
</tr>
<tr>
<td>0.6 < |eta| < 0.8</td>
<td><a href="nom5to8eta6to8iasPred.png"><img height="300" src="nom5to8eta6to8iasPred.png"></a></td>
<td><a href="nom9to10eta6to8iasPred.png"><img height="300" src="nom9to10eta6to8iasPred.png"></a></td>
<td><a href="nom11to12eta6to8iasPred.png"><img height="300" src="nom11to12eta6to8iasPred.png"></a></td>
<td><a href="nom13to14eta6to8iasPred.png"><img height="300" src="nom13to14eta6to8iasPred.png"></a></td>
<td><a href="nom15to16eta6to8iasPred.png"><img height="300" src="nom15to16eta6to8iasPred.png"></a></td>
<td><a href="nom17to18eta6to8iasPred.png"><img height="300" src="nom17to18eta6to8iasPred.png"></a></td>
<td>0.6 < |eta| < 0.8</td>
</tr>
<tr>
<td></td>
<td>NoM 5-8</td>
<td>NoM 9-10</td>
<td>NoM 11-12</td>
<td>NoM 13-14</td>
<td>NoM 15-16</td>
<td>NoM 17+</td>
</tr>
<tr>
<td>0.8 < |eta| < 1.0</td>
<td><a href="nom5to8eta8to10iasPred.png"><img height="300" src="nom5to8eta8to10iasPred.png"></a></td>
<td><a href="nom9to10eta8to10iasPred.png"><img height="300" src="nom9to10eta8to10iasPred.png"></a></td>
<td><a href="nom11to12eta8to10iasPred.png"><img height="300" src="nom11to12eta8to10iasPred.png"></a></td>
<td><a href="nom13to14eta8to10iasPred.png"><img height="300" src="nom13to14eta8to10iasPred.png"></a></td>
<td><a href="nom15to16eta8to10iasPred.png"><img height="300" src="nom15to16eta8to10iasPred.png"></a></td>
<td><a href="nom17to18eta8to10iasPred.png"><img height="300" src="nom17to18eta8to10iasPred.png"></a></td>
<td>0.8 < |eta| < 1.0</td>
</tr>
<tr>
<td>1.0 < |eta| < 1.2</td>
<td><a href="nom5to8eta10to12iasPred.png"><img height="300" src="nom5to8eta10to12iasPred.png"></a></td>
<td><a href="nom9to10eta10to12iasPred.png"><img height="300" src="nom9to10eta10to12iasPred.png"></a></td>
<td><a href="nom11to12eta10to12iasPred.png"><img height="300" src="nom11to12eta10to12iasPred.png"></a></td>
<td><a href="nom13to14eta10to12iasPred.png"><img height="300" src="nom13to14eta10to12iasPred.png"></a></td>
<td><a href="nom15to16eta10to12iasPred.png"><img height="300" src="nom15to16eta10to12iasPred.png"></a></td>
<td><a href="nom17to18eta10to12iasPred.png"><img height="300" src="nom17to18eta10to12iasPred.png"></a></td>
<td>1.0 < |eta| < 1.2</td>
</tr>
<tr>
<td>1.2 < |eta| < 1.4</td>
<td><a href="nom5to8eta12to14iasPred.png"><img height="300" src="nom5to8eta12to14iasPred.png"></a></td>
<td><a href="nom9to10eta12to14iasPred.png"><img height="300" src="nom9to10eta12to14iasPred.png"></a></td>
<td><a href="nom11to12eta12to14iasPred.png"><img height="300" src="nom11to12eta12to14iasPred.png"></a></td>
<td><a href="nom13to14eta12to14iasPred.png"><img height="300" src="nom13to14eta12to14iasPred.png"></a></td>
<td><a href="nom15to16eta12to14iasPred.png"><img height="300" src="nom15to16eta12to14iasPred.png"></a></td>
<td><a href="nom17to18eta12to14iasPred.png"><img height="300" src="nom17to18eta12to14iasPred.png"></a></td>
<td>1.2 < |eta| < 1.4</td>
</tr>
<tr>
<td>1.4 < |eta| < 1.6</td>
<td><a href="nom5to8eta14to16iasPred.png"><img height="300" src="nom5to8eta14to16iasPred.png"></a></td>
<td><a href="nom9to10eta14to16iasPred.png"><img height="300" src="nom9to10eta14to16iasPred.png"></a></td>
<td><a href="nom11to12eta14to16iasPred.png"><img height="300" src="nom11to12eta14to16iasPred.png"></a></td>
<td><a href="nom13to14eta14to16iasPred.png"><img height="300" src="nom13to14eta14to16iasPred.png"></a></td>
<td><a href="nom15to16eta14to16iasPred.png"><img height="300" src="nom15to16eta14to16iasPred.png"></a></td>
<td><a href="nom17to18eta14to16iasPred.png"><img height="300" src="nom17to18eta14to16iasPred.png"></a></td>
<td>1.4 < |eta| < 1.6</td>
</tr>
<tr>
<td></td>
<td>NoM 5-8</td>
<td>NoM 9-10</td>
<td>NoM 11-12</td>
<td>NoM 13-14</td>
<td>NoM 15-16</td>
<td>NoM 17+</td>
</tr>
</table>

</FONT>
</BODY>
</HTML>

EOF




exit

