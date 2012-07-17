#!/bin/bash

usage='Usage: makeDeDxPlots.sh -f fileName -d directoryName'
  
args=`getopt fd: -- "$@"`
if test $? != 0
     then
         echo $usage
         exit 1
fi

eval set -- "$args"
for i
  do
  case "$i" in
      -f) shift; file_name=$2;shift;;
      -d) shift; dir_name=$2;shift;;
  esac
done

if [ "X"${file_name} == "X" ]
    then
    echo "invalid root file name! Please give a valid one."
    echo $usage
    exit
fi
if [ "X"${dir_name} == "X" ]
    then
    echo "invalid directory name! Please give a valid one."
    echo $usage
    exit
fi

myPwd=`pwd`
mkdir -p ${dir_name}
cd ${dir_name}

# run root in batch
root -b $myPwd/$file_name <<!
.x $myPwd/makeDeDxPlots.C
.q
!


cat > index.html <<EOF

<HTML>

<HEAD><TITLE>dE/dx Plots</TITLE></HEAD>
 
<BODY link="Red">
<FONT color="Black">
 
<Center>
<h2><FONT color="Blue"> TK E1 dE/dx distributions </FONT></h2>
</Center>

<table>
<tr>
<td></td>
<td>NoM 1-5</td>
<td>NoM 6-10</td>
<td>NoM 11-15</td>
<td>NoM 16-20</td>
<td>NoM 21+</td>
</tr>
<tr>
<td>Pt 15-50 GeV</td>
<td><a href="histForPDFNoMPtBinsIndex1.png"><img height="200" src="histForPDFNoMPtBinsIndex1.png"></a></td>
<td><a href="histForPDFNoMPtBinsIndex2.png"><img height="200" src="histForPDFNoMPtBinsIndex2.png"></a></td>
<td><a href="histForPDFNoMPtBinsIndex3.png"><img height="200" src="histForPDFNoMPtBinsIndex3.png"></a></td>
<td><a href="histForPDFNoMPtBinsIndex4.png"><img height="200" src="histForPDFNoMPtBinsIndex4.png"></a></td>
<td><a href="histForPDFNoMPtBinsIndex5.png"><img height="200" src="histForPDFNoMPtBinsIndex5.png"></a></td>
</tr>
<tr>
<td>Pt 50-150 GeV</td>
<td><a href="histForPDFNoMPtBinsIndex6.png"><img height="200" src="histForPDFNoMPtBinsIndex6.png"></a></td>
<td><a href="histForPDFNoMPtBinsIndex7.png"><img height="200" src="histForPDFNoMPtBinsIndex7.png"></a></td>
<td><a href="histForPDFNoMPtBinsIndex8.png"><img height="200" src="histForPDFNoMPtBinsIndex8.png"></a></td>
<td><a href="histForPDFNoMPtBinsIndex9.png"><img height="200" src="histForPDFNoMPtBinsIndex9.png"></a></td>
<td><a href="histForPDFNoMPtBinsIndex10.png"><img height="200" src="histForPDFNoMPtBinsIndex10.png"></a></td>
</tr>
<tr>
<td>Pt 150-300 GeV</td>
<td><a href="histForPDFNoMPtBinsIndex11.png"><img height="200" src="histForPDFNoMPtBinsIndex11.png"></a></td>
<td><a href="histForPDFNoMPtBinsIndex12.png"><img height="200" src="histForPDFNoMPtBinsIndex12.png"></a></td>
<td><a href="histForPDFNoMPtBinsIndex13.png"><img height="200" src="histForPDFNoMPtBinsIndex13.png"></a></td>
<td><a href="histForPDFNoMPtBinsIndex14.png"><img height="200" src="histForPDFNoMPtBinsIndex14.png"></a></td>
<td><a href="histForPDFNoMPtBinsIndex15.png"><img height="200" src="histForPDFNoMPtBinsIndex15.png"></a></td>
</tr>
<tr>
<td>Pt 300-450 GeV</td>
<td><a href="histForPDFNoMPtBinsIndex16.png"><img height="200" src="histForPDFNoMPtBinsIndex16.png"></a></td>
<td><a href="histForPDFNoMPtBinsIndex17.png"><img height="200" src="histForPDFNoMPtBinsIndex17.png"></a></td>
<td><a href="histForPDFNoMPtBinsIndex18.png"><img height="200" src="histForPDFNoMPtBinsIndex18.png"></a></td>
<td><a href="histForPDFNoMPtBinsIndex19.png"><img height="200" src="histForPDFNoMPtBinsIndex19.png"></a></td>
<td><a href="histForPDFNoMPtBinsIndex20.png"><img height="200" src="histForPDFNoMPtBinsIndex20.png"></a></td>
</tr>
<tr>
<td>Pt > 450 GeV</td>
<td><a href="histForPDFNoMPtBinsIndex21.png"><img height="200" src="histForPDFNoMPtBinsIndex21.png"></a></td>
<td><a href="histForPDFNoMPtBinsIndex22.png"><img height="200" src="histForPDFNoMPtBinsIndex22.png"></a></td>
<td><a href="histForPDFNoMPtBinsIndex23.png"><img height="200" src="histForPDFNoMPtBinsIndex23.png"></a></td>
<td><a href="histForPDFNoMPtBinsIndex24.png"><img height="200" src="histForPDFNoMPtBinsIndex24.png"></a></td>
<td><a href="histForPDFNoMPtBinsIndex25.png"><img height="200" src="histForPDFNoMPtBinsIndex25.png"></a></td>
</tr>
</table>

</FONT>
</BODY>
</HTML>

EOF


exit

