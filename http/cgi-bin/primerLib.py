import string
import math
import Levenshtein

def get_melting_temp(enthalpy, entropy, selfcomp, gccont, plen):
	R=1.9872
#	R=3.327
	molarConc=0.0004
	nacl=0.05
	cfac=((4.29*gccont-3.95)*math.log(nacl)+0.94*(math.log(nacl))**2)*10**-5
	if (selfcomp):
		x = 1.0;
	else:
		x = 4.0;
	Tm = enthalpy * 1000 / (entropy + R*math.log(molarConc / x)) - 273.15;
#	Tm=1/((1/TmB)+cfac)
	Tm=Tm*0.96
#	Tm=Tm*0.947
	return Tm;

def seqCheck(templ, prim):
	j=0;
	lp=len(prim);
	txt3=prim[::-1]
	temLen=int(len(templ)/2)
	txt4='';
	for i in range(0, len(txt3)):
		if txt3[i]=='G':
			txt4+='C';
		elif txt3[i]=='C':
			txt4+='G';
		elif txt3[i]=='A':
			txt4+='T';
		elif txt3[i]=='T':
			txt4+='A';
	txt3=''
	txt3=txt4;
	txt4=''
	for i in range(0, len(templ)-lp-2):
		if Levenshtein.ratio(templ[i:i+lp+2], prim) > 0.83:
			j=j+1
	for i in range(0, len(templ)-lp-2):
		if Levenshtein.ratio(templ[i:i+lp+2], txt3) > 0.83:
			j=j+1
	if j>14:
		return 0;
	else:
		return 1;

def get_enthalpy(seq):
	nnenthalpy=0.0;
	i=0
	enthalpy = 0.2;
	for m in seq:
		i=i+1
		if ((i==1 or i==len(seq)) and (m=='A' or m=='T')):
			enthalpy=enthalpy+2.2
		if (i<len(seq)):
			n=seq[i]
			if (m == 'A'):
				if (n == 'A'):
					nnenthalpy = -7.6;
				elif (n == 'T'):
					nnenthalpy = -7.2;
				elif (n == 'C'):
					nnenthalpy = -8.4;
				elif (n == 'G'):
					nnenthalpy = -7.8;
			elif (m == 'T'):
				if (n == 'A'):
					nnenthalpy = -7.2;
				elif (n == 'T'):
					nnenthalpy = -7.6;
				elif (n == 'C'):
					nnenthalpy = -8.2;
				elif (n == 'G'):
					nnenthalpy = -8.5;
			elif (m == 'C'):
				if (n == 'A'):
					nnenthalpy = -8.5;
				elif (n == 'T'):
					nnenthalpy = -7.8;
				elif (n == 'C'):
					nnenthalpy = -8.0;
				elif (n == 'G'):
					nnenthalpy = -10.6;
			elif (m == 'G'):
				if (n == 'A'):
					nnenthalpy = -8.2;
				elif (n == 'T'):
					nnenthalpy = -8.4;
				elif (n == 'C'):
					nnenthalpy = -9.8;
				elif (n == 'G'):
					nnenthalpy = -8.0;
		else:
			nnenthalpy=0.0
		enthalpy = enthalpy + nnenthalpy
	return enthalpy;

def get_entropy(seq):
	R=1.9872
	molarConc=0.0004
	nacl=0.05
	cfac=((4.29*GCcontent(seq)-3.95)*math.log(nacl)+0.94*(math.log(nacl))**2)*10**-5
	mg=0.002
	aa=0.0000392
	ba=-0.00000911
	ca=0.0000626
	da=0.0000142
	ea=-0.000482
	fa=0.000525
	ga=0.0000831
	mgfac=aa+ba*math.log(mg)+GCcontent(seq)*(ca+da*math.log(mg))+(ea+fa*math.log(mg)+ga*(math.log(mg))**2)/(2*(len(seq)-1))
	nnentropy=0.0;
	i=0
	entropy = -5.7;
	for m in seq:
		i=i+1
		if ((i==1 or i==len(seq)) and (m=='A' or m=='T')):
			entropy=entropy+6.9;
		if (i<len(seq)):
			n=seq[i]
			if (m == 'A'):
				if (n == 'A'):
					nnentropy = -21.3;
				elif (n == 'T'):
					nnentropy = -20.4;
				elif (n == 'C'):
					nnentropy = -22.4;
				elif (n == 'G'):
					nnentropy = -21.0;
			elif (m == 'T'):
				if (n == 'A'):
					nnentropy = -21.3;
				elif (n == 'T'):
					nnentropy = -21.3;
				elif (n == 'C'):
					nnentropy = -22.2;
				elif (n == 'G'):
					nnentropy = -22.7;
			elif (m == 'C'):
				if (n == 'A'):
					nnentropy = -22.7;
				elif (n == 'T'):
					nnentropy = -21.0;
				elif (n == 'C'):
					nnentropy = -19.9;
				elif (n == 'G'):
					nnentropy = -27.2;
			elif (m == 'G'):
				if (n == 'A'):
					nnentropy = -22.2;
				elif (n == 'T'):
					nnentropy = -22.4;
				elif (n == 'C'):
					nnentropy = -24.4;
				elif (n == 'G'):
					nnentropy = -19.9;
		else:
			nnentropy=0.0
		entropy = entropy + nnentropy
	entropy=entropy+get_enthalpy(seq)*mgfac*1000
	entropy=entropy+0.368*len(seq)*math.log(nacl)
	return entropy;

def check_self_comp(seq):
	checker = 0;
	seqRC=[]
	seqRev=list(reversed(seq));
	seqOri=list(reversed(seqRev));
	for i in seqRev:
		if i=='A':
			seqRC+='T'
		elif i=='T':
			seqRC+='A'
		elif i=='C':
			seqRC+='G'
		elif i=='G':
			seqRC+='C'
	if seqOri==seqRC:
		checker=1;
	return checker;

def GCcontent(seq):
	content = (seq.count('G')+seq.count('C'))/len(seq);
	return content;

