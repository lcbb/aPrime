#!/usr/local/bin/python3.4
from primerLib import *
import string
import math

def primer_tm(temp, length):
	primers=[]
	primers.append([])
	primers.append([])
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
										primers[0].append(testOne)
										primers[1].append(testTwo)
	return primers

primers = [];
fp=open('/home/trsheph/aPCRprog/in.txt', 'r');
fout=open('/home/trsheph/aPCRprog/out.txt', 'w');
templDNA=''
for line in fp:
	newline=line.rstrip('\r\n')
	templDNA=templDNA+newline;
print(templDNA[:int(len(templDNA)/2)]);
fp.close();
scLen=1000;
for i in range(scLen, scLen+20):
	targetLen=i;
	primers=primer_tm(templDNA, targetLen);
	print(primers)
fout.close();
