#!/bin/bash

usage='Usage: makeIasPredictionDataPlotsHtml.sh -d plotDirName -r rootFile'
  
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

.x /local/cms/user/cooper/cmssw/428/HSCPAnalysis/src/HSCP2011/ToyMC/bin/makeIasPredictionDataPngsFixedBins.C

.q
!

cat > index.html <<EOF

<HTML>

<HEAD><TITLE>Ias prediction (limits) and Data</TITLE></HEAD>
 
<BODY link="Red">
<FONT color="Black">
 
<Center>
<h2><FONT color="Blue">Ias prediction nominal (limits)</FONT><FONT color="Red"> Ias prediction (-1 sigma limits)</FONT><FONT color="Green"> Ias prediction (+1 sigma limits)</FONT> and Data</h2>
</Center>

<table>
<tr>
<td></td>
<td>NoM 5-8</td>
<td>NoM 9-10</td>
<td>NoM 11-12</td>
<td>NoM 13-14</td>
<td>NoM 15-16</td>
<td>NoM 17-18</td>
<td>NoM 19-20</td>
<td>NoM 21-22</td>
</tr>
<tr>
<td>0.0 <= |eta| < 0.2</td>
<td><a href="nom5to8eta0to2iasDataPred.png"><img height="200" src="nom5to8eta0to2iasDataPred.png"></a></td>
<td><a href="nom9to10eta0to2iasDataPred.png"><img height="200" src="nom9to10eta0to2iasDataPred.png"></a></td>
<td><a href="nom11to12eta0to2iasDataPred.png"><img height="200" src="nom11to12eta0to2iasDataPred.png"></a></td>
<td><a href="nom13to14eta0to2iasDataPred.png"><img height="200" src="nom13to14eta0to2iasDataPred.png"></a></td>
<td><a href="nom15to16eta0to2iasDataPred.png"><img height="200" src="nom15to16eta0to2iasDataPred.png"></a></td>
<td><a href="nom17to18eta0to2iasDataPred.png"><img height="200" src="nom17to18eta0to2iasDataPred.png"></a></td>
<td><a href="nom19to20eta0to2iasDataPred.png"><img height="200" src="nom19to20eta0to2iasDataPred.png"></a></td>
<td><a href="nom21to22eta0to2iasDataPred.png"><img height="200" src="nom21to22eta0to2iasDataPred.png"></a></td>
<td>0.0 <= |eta| < 0.2</td>
</tr>
<tr>
<td>0.2 < |eta| < 0.4</td>
<td><a href="nom5to8eta2to4iasDataPred.png"><img height="200" src="nom5to8eta2to4iasDataPred.png"></a></td>
<td><a href="nom9to10eta2to4iasDataPred.png"><img height="200" src="nom9to10eta2to4iasDataPred.png"></a></td>
<td><a href="nom11to12eta2to4iasDataPred.png"><img height="200" src="nom11to12eta2to4iasDataPred.png"></a></td>
<td><a href="nom13to14eta2to4iasDataPred.png"><img height="200" src="nom13to14eta2to4iasDataPred.png"></a></td>
<td><a href="nom15to16eta2to4iasDataPred.png"><img height="200" src="nom15to16eta2to4iasDataPred.png"></a></td>
<td><a href="nom17to18eta2to4iasDataPred.png"><img height="200" src="nom17to18eta2to4iasDataPred.png"></a></td>
<td><a href="nom19to20eta2to4iasDataPred.png"><img height="200" src="nom19to20eta2to4iasDataPred.png"></a></td>
<td><a href="nom21to22eta2to4iasDataPred.png"><img height="200" src="nom21to22eta2to4iasDataPred.png"></a></td>
<td>0.2 < |eta| < 0.4</td>
</tr>
<tr>
<td>0.4 < |eta| < 0.6</td>
<td><a href="nom5to8eta4to6iasDataPred.png"><img height="200" src="nom5to8eta4to6iasDataPred.png"></a></td>
<td><a href="nom9to10eta4to6iasDataPred.png"><img height="200" src="nom9to10eta4to6iasDataPred.png"></a></td>
<td><a href="nom11to12eta4to6iasDataPred.png"><img height="200" src="nom11to12eta4to6iasDataPred.png"></a></td>
<td><a href="nom13to14eta4to6iasDataPred.png"><img height="200" src="nom13to14eta4to6iasDataPred.png"></a></td>
<td><a href="nom15to16eta4to6iasDataPred.png"><img height="200" src="nom15to16eta4to6iasDataPred.png"></a></td>
<td><a href="nom17to18eta4to6iasDataPred.png"><img height="200" src="nom17to18eta4to6iasDataPred.png"></a></td>
<td><a href="nom19to20eta4to6iasDataPred.png"><img height="200" src="nom19to20eta4to6iasDataPred.png"></a></td>
<td><a href="nom21to22eta4to6iasDataPred.png"><img height="200" src="nom21to22eta4to6iasDataPred.png"></a></td>
<td>0.4 < |eta| < 0.6</td>
</tr>
<tr>
<td>0.6 < |eta| < 0.8</td>
<td><a href="nom5to8eta6to8iasDataPred.png"><img height="200" src="nom5to8eta6to8iasDataPred.png"></a></td>
<td><a href="nom9to10eta6to8iasDataPred.png"><img height="200" src="nom9to10eta6to8iasDataPred.png"></a></td>
<td><a href="nom11to12eta6to8iasDataPred.png"><img height="200" src="nom11to12eta6to8iasDataPred.png"></a></td>
<td><a href="nom13to14eta6to8iasDataPred.png"><img height="200" src="nom13to14eta6to8iasDataPred.png"></a></td>
<td><a href="nom15to16eta6to8iasDataPred.png"><img height="200" src="nom15to16eta6to8iasDataPred.png"></a></td>
<td><a href="nom17to18eta6to8iasDataPred.png"><img height="200" src="nom17to18eta6to8iasDataPred.png"></a></td>
<td><a href="nom19to20eta6to8iasDataPred.png"><img height="200" src="nom19to20eta6to8iasDataPred.png"></a></td>
<td><a href="nom21to22eta6to8iasDataPred.png"><img height="200" src="nom21to22eta6to8iasDataPred.png"></a></td>
<td>0.6 < |eta| < 0.8</td>
</tr>
<tr>
<td></td>
<td>NoM 5-8</td>
<td>NoM 9-10</td>
<td>NoM 11-12</td>
<td>NoM 13-14</td>
<td>NoM 15-16</td>
<td>NoM 17-18</td>
<td>NoM 19-20</td>
<td>NoM 21-22</td>
</tr>
<tr>
<td>0.8 < |eta| < 1.0</td>
<td><a href="nom5to8eta8to10iasDataPred.png"><img height="200" src="nom5to8eta8to10iasDataPred.png"></a></td>
<td><a href="nom9to10eta8to10iasDataPred.png"><img height="200" src="nom9to10eta8to10iasDataPred.png"></a></td>
<td><a href="nom11to12eta8to10iasDataPred.png"><img height="200" src="nom11to12eta8to10iasDataPred.png"></a></td>
<td><a href="nom13to14eta8to10iasDataPred.png"><img height="200" src="nom13to14eta8to10iasDataPred.png"></a></td>
<td><a href="nom15to16eta8to10iasDataPred.png"><img height="200" src="nom15to16eta8to10iasDataPred.png"></a></td>
<td><a href="nom17to18eta8to10iasDataPred.png"><img height="200" src="nom17to18eta8to10iasDataPred.png"></a></td>
<td><a href="nom19to20eta8to10iasDataPred.png"><img height="200" src="nom19to20eta8to10iasDataPred.png"></a></td>
<td><a href="nom21to22eta8to10iasDataPred.png"><img height="200" src="nom21to22eta8to10iasDataPred.png"></a></td>
<td>0.8 < |eta| < 1.0</td>
</tr>
<tr>
<td>1.0 < |eta| < 1.2</td>
<td><a href="nom5to8eta10to12iasDataPred.png"><img height="200" src="nom5to8eta10to12iasDataPred.png"></a></td>
<td><a href="nom9to10eta10to12iasDataPred.png"><img height="200" src="nom9to10eta10to12iasDataPred.png"></a></td>
<td><a href="nom11to12eta10to12iasDataPred.png"><img height="200" src="nom11to12eta10to12iasDataPred.png"></a></td>
<td><a href="nom13to14eta10to12iasDataPred.png"><img height="200" src="nom13to14eta10to12iasDataPred.png"></a></td>
<td><a href="nom15to16eta10to12iasDataPred.png"><img height="200" src="nom15to16eta10to12iasDataPred.png"></a></td>
<td><a href="nom17to18eta10to12iasDataPred.png"><img height="200" src="nom17to18eta10to12iasDataPred.png"></a></td>
<td><a href="nom19to20eta10to12iasDataPred.png"><img height="200" src="nom19to20eta10to12iasDataPred.png"></a></td>
<td><a href="nom21to22eta10to12iasDataPred.png"><img height="200" src="nom21to22eta10to12iasDataPred.png"></a></td>
<td>1.0 < |eta| < 1.2</td>
</tr>
<tr>
<td>1.2 < |eta| < 1.4</td>
<td><a href="nom5to8eta12to14iasDataPred.png"><img height="200" src="nom5to8eta12to14iasDataPred.png"></a></td>
<td><a href="nom9to10eta12to14iasDataPred.png"><img height="200" src="nom9to10eta12to14iasDataPred.png"></a></td>
<td><a href="nom11to12eta12to14iasDataPred.png"><img height="200" src="nom11to12eta12to14iasDataPred.png"></a></td>
<td><a href="nom13to14eta12to14iasDataPred.png"><img height="200" src="nom13to14eta12to14iasDataPred.png"></a></td>
<td><a href="nom15to16eta12to14iasDataPred.png"><img height="200" src="nom15to16eta12to14iasDataPred.png"></a></td>
<td><a href="nom17to18eta12to14iasDataPred.png"><img height="200" src="nom17to18eta12to14iasDataPred.png"></a></td>
<td><a href="nom19to20eta12to14iasDataPred.png"><img height="200" src="nom19to20eta12to14iasDataPred.png"></a></td>
<td><a href="nom21to22eta12to14iasDataPred.png"><img height="200" src="nom21to22eta12to14iasDataPred.png"></a></td>
<td>1.2 < |eta| < 1.4</td>
</tr>
<tr>
<td>1.4 < |eta| < 1.6</td>
<td><a href="nom5to8eta14to16iasDataPred.png"><img height="200" src="nom5to8eta14to16iasDataPred.png"></a></td>
<td><a href="nom9to10eta14to16iasDataPred.png"><img height="200" src="nom9to10eta14to16iasDataPred.png"></a></td>
<td><a href="nom11to12eta14to16iasDataPred.png"><img height="200" src="nom11to12eta14to16iasDataPred.png"></a></td>
<td><a href="nom13to14eta14to16iasDataPred.png"><img height="200" src="nom13to14eta14to16iasDataPred.png"></a></td>
<td><a href="nom15to16eta14to16iasDataPred.png"><img height="200" src="nom15to16eta14to16iasDataPred.png"></a></td>
<td><a href="nom17to18eta14to16iasDataPred.png"><img height="200" src="nom17to18eta14to16iasDataPred.png"></a></td>
<td><a href="nom19to20eta14to16iasDataPred.png"><img height="200" src="nom19to20eta14to16iasDataPred.png"></a></td>
<td><a href="nom21to22eta14to16iasDataPred.png"><img height="200" src="nom21to22eta14to16iasDataPred.png"></a></td>
<td>1.4 < |eta| < 1.6</td>
</tr>
<tr>
<td>1.6 < |eta| < 1.8</td>
<td><a href="nom5to8eta16to18iasDataPred.png"><img height="200" src="nom5to8eta16to18iasDataPred.png"></a></td>
<td><a href="nom9to10eta16to18iasDataPred.png"><img height="200" src="nom9to10eta16to18iasDataPred.png"></a></td>
<td><a href="nom11to12eta16to18iasDataPred.png"><img height="200" src="nom11to12eta16to18iasDataPred.png"></a></td>
<td><a href="nom13to14eta16to18iasDataPred.png"><img height="200" src="nom13to14eta16to18iasDataPred.png"></a></td>
<td><a href="nom15to16eta16to18iasDataPred.png"><img height="200" src="nom15to16eta16to18iasDataPred.png"></a></td>
<td><a href="nom17to18eta16to18iasDataPred.png"><img height="200" src="nom17to18eta16to18iasDataPred.png"></a></td>
<td><a href="nom19to20eta16to18iasDataPred.png"><img height="200" src="nom19to20eta16to18iasDataPred.png"></a></td>
<td><a href="nom21to22eta16to18iasDataPred.png"><img height="200" src="nom21to22eta16to18iasDataPred.png"></a></td>
<td>1.6 < |eta| < 1.8</td>
</tr>
<tr>
<td>1.8 < |eta| < 2.0</td>
<td><a href="nom5to8eta18to20iasDataPred.png"><img height="200" src="nom5to8eta18to20iasDataPred.png"></a></td>
<td><a href="nom9to10eta18to20iasDataPred.png"><img height="200" src="nom9to10eta18to20iasDataPred.png"></a></td>
<td><a href="nom11to12eta18to20iasDataPred.png"><img height="200" src="nom11to12eta18to20iasDataPred.png"></a></td>
<td><a href="nom13to14eta18to20iasDataPred.png"><img height="200" src="nom13to14eta18to20iasDataPred.png"></a></td>
<td><a href="nom15to16eta18to20iasDataPred.png"><img height="200" src="nom15to16eta18to20iasDataPred.png"></a></td>
<td><a href="nom17to18eta18to20iasDataPred.png"><img height="200" src="nom17to18eta18to20iasDataPred.png"></a></td>
<td><a href="nom19to20eta18to20iasDataPred.png"><img height="200" src="nom19to20eta18to20iasDataPred.png"></a></td>
<td><a href="nom21to22eta18to20iasDataPred.png"><img height="200" src="nom21to22eta18to20iasDataPred.png"></a></td>
<td>1.8 < |eta| < 2.0</td>
</tr>
<tr>
<td>2.0 < |eta| < 2.2</td>
<td><a href="nom5to8eta20to22iasDataPred.png"><img height="200" src="nom5to8eta20to22iasDataPred.png"></a></td>
<td><a href="nom9to10eta20to22iasDataPred.png"><img height="200" src="nom9to10eta20to22iasDataPred.png"></a></td>
<td><a href="nom11to12eta20to22iasDataPred.png"><img height="200" src="nom11to12eta20to22iasDataPred.png"></a></td>
<td><a href="nom13to14eta20to22iasDataPred.png"><img height="200" src="nom13to14eta20to22iasDataPred.png"></a></td>
<td><a href="nom15to16eta20to22iasDataPred.png"><img height="200" src="nom15to16eta20to22iasDataPred.png"></a></td>
<td><a href="nom17to18eta20to22iasDataPred.png"><img height="200" src="nom17to18eta20to22iasDataPred.png"></a></td>
<td><a href="nom19to20eta20to22iasDataPred.png"><img height="200" src="nom19to20eta20to22iasDataPred.png"></a></td>
<td><a href="nom21to22eta20to22iasDataPred.png"><img height="200" src="nom21to22eta20to22iasDataPred.png"></a></td>
<td>2.0 < |eta| < 2.2</td>
</tr>
<tr>
<td>2.2 < |eta| < 2.4</td>
<td><a href="nom5to8eta22to24iasDataPred.png"><img height="200" src="nom5to8eta22to24iasDataPred.png"></a></td>
<td><a href="nom9to10eta22to24iasDataPred.png"><img height="200" src="nom9to10eta22to24iasDataPred.png"></a></td>
<td><a href="nom11to12eta22to24iasDataPred.png"><img height="200" src="nom11to12eta22to24iasDataPred.png"></a></td>
<td><a href="nom13to14eta22to24iasDataPred.png"><img height="200" src="nom13to14eta22to24iasDataPred.png"></a></td>
<td><a href="nom15to16eta22to24iasDataPred.png"><img height="200" src="nom15to16eta22to24iasDataPred.png"></a></td>
<td><a href="nom17to18eta22to24iasDataPred.png"><img height="200" src="nom17to18eta22to24iasDataPred.png"></a></td>
<td><a href="nom19to20eta22to24iasDataPred.png"><img height="200" src="nom19to20eta22to24iasDataPred.png"></a></td>
<td><a href="nom21to22eta22to24iasDataPred.png"><img height="200" src="nom21to22eta22to24iasDataPred.png"></a></td>
<td>2.2 < |eta| < 2.4</td>
</tr>
<tr>
<td></td>
<td>NoM 5-8</td>
<td>NoM 9-10</td>
<td>NoM 11-12</td>
<td>NoM 13-14</td>
<td>NoM 15-16</td>
<td>NoM 17-18</td>
<td>NoM 19-20</td>
<td>NoM 21-22</td>
</tr>
</table>

</FONT>
</BODY>
</HTML>

EOF




exit

