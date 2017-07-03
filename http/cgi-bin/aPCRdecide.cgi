#!/usr/local/bin/python3.4
#
# (C) Tyson R. Shepherd, PhD
#
import math
from primerLib import *
import os, sys
import cgi
import cgitb; cgitb.enable()
#
#
#
def primer_tm(temp, length):
	primers=[]
	primers.append([])
	primers.append([])
	primers.append([])
	primers.append([])
#	primers.append([])
	testOne = '';
	testTwo = '';
	count = 0;
	for p in range(18,23):
		for q in range(0, int(len(temp)/2)):
			testOne=temp[q:p+q];
			Tm1 = get_melting_temp(get_enthalpy(testOne), get_entropy(testOne), check_self_comp(testOne), GCcontent(testOne), len(testOne));
			if ((Tm1 >= 55) and (Tm1 <= 57)):
				if (((testOne[p-1]=='A') or (testOne[p-1]=='T')) and ((testOne[p-2]=='A') or (testOne[p-2]=='T')) and (GCcontent(testOne[:9]) > GCcontent(testOne[len(testOne)-9:]))):
					if (((GCcontent(testOne) >= 0.4) and (GCcontent(testOne) <= 0.6)) and ((testOne[0]=='C') or (testOne[0]=='G'))):
						for t in range(18,23):
							testTwo=temp[q+length-t:q+length];
							Tm2 = get_melting_temp(get_enthalpy(testTwo), get_entropy(testTwo), check_self_comp(testTwo), GCcontent(testTwo), len(testTwo));
							if ((Tm2 >= Tm1+1) and (Tm2 <= Tm1+3)):
								if ((GCcontent(testTwo) >= 0.4) and (GCcontent(testTwo) <= 0.6)):
									if (((testTwo[t - 1] == 'C') or (testTwo[t - 1] == 'G')) and ((testTwo[0] == 'C') or (testTwo[0] == 'G'))):
#										if (int((seqCheck(temp, testOne)))==1 and int((seqCheck(temp, testTwo)))==1):
										if 1==1:
											if not (testTwo in primers[1]):
												primers[0].append(testOne)
												primers[1].append(testTwo)
												primers[2].append(str(Tm1)[:4])
#												primers[3].append(str(GCcontent(temp[q:q+length])*100)[:4])
												primers[3].append(str(temp[q:q+length]))
	return primers

print("Content-Type: text/html\n\n")
print("""
<html>
<head>
<style>
   body {
        font-family:courier;
   }
   td {
	padding: 10px;
	border: 1px solid black;
   }
   th {
	text-align: left;
	padding: 10px;
	border: 1px solid black;
	background-color: #ffffe0;
   }
</style>
<title>Output</title>
</head>
<body >
<h3>Output</h3>
""");
#
# Process form to python
#
form = cgi.FieldStorage();
scLen = int(form["sLen"].value);
templType=int(form["templT"].value);
if (templType==1):
	fInStr='./m13.txt'
	if (scLen>7248):
		print('Scaffold must be less than 7248 for M13')
		print('</body>');
		print('</html>');
		exit();
elif (templType==2):
	fInStr='./lambda.txt'
	if (scLen>48000):
		print('Scaffold must be less than 48000 for Lambda')
		print('</body>');
		print('</html>');
		exit();
else:
	templDNA=str(form["utempl"].value)+str(form["utempl"].value);
	if (scLen>int(len(templDNA)/2-1)):
		print('Desired scaffold length is greater than the provided user sequence');
		print('</body>');
		print('</html>');
		exit();
if (templType!=3):
	fp=open(fInStr, 'r');
	templDNA=''
	for line in fp:
		newline=line.rstrip('\r\n')
		templDNA=templDNA+newline;
	fp.close();
print('<table width=100%>');
print('<tr><th>Forward Primer</th><th>Reverse Primer</th><th>Tm</th><th>Scaffold Sequence</th></tr>');
for i in range(scLen, scLen+1):
	print('<tr>');
	targetLen=i;
	primers=primer_tm(templDNA, targetLen);
	print('<td>');
	for j in primers[0]:
		print(str(j)+'<br>');
	print('</td><td>');
	for j in primers[1]:
		print(str(j)+'<br>');
	print('</td><td>');
	for j in primers[2]:
		print(str(j)+'<br>');
	print('</td><td>');
	for j in primers[3]:
		print(str(j)+'<br>');
	print('</td>');
#	for j in primers[4]:
#		print(str(j)+'<br>');
#	print('</td><td>');
	print('</tr>');
print('</table>');
print('<br><br>Template: <br>');
print('<table width=100% style="table-layout: fixed; width: 100%"><tr><td style="word-wrap: break-word">');
print(templDNA[:int(len(templDNA)/2)]+'<br>');
print('</td></tr></table><br>');
#print('Desired scaffold length: '+str(scLen)+'<br>');
#print('File to read from: '+fInStr+'<br>');
print('</body>');
print('</html>');
