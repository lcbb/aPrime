#!/usr/local/bin/python3.4
import sqlite3
import cgi
import cgitb; cgitb.enable()
def getRows(scLen):
	conn = sqlite3.connect('aPrimerM13DB.db')
	c=conn.cursor()
	t=(scLen,)
	n=0
	rct=c.execute("SELECT COUNT(*) from primers where length=?", t)
	rowCount=rct.fetchone()[0]
	rowsA=c.execute('SELECT * FROM primers WHERE length=?', t);
	while rowCount==0:
		n=n+1
		t=(scLen+n,)
		rct=c.execute("SELECT COUNT(*) from primers where length=?", t)
		rowCount=rct.fetchone()[0]
		rowsA=c.execute('SELECT * FROM primers WHERE length=?', t);
	conn.close
	return(rowsA)
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
<title>aPrime Output</title>
</head>
<body>
""");
form = cgi.FieldStorage();
scLen = int(form["sLen"].value);
templType=int(form["templT"].value);
#
# Process form to python
#
if (templType==1):
	fInStr='./m13.txt'
	templName='M13'
	templAmount='30ng'
	if (scLen>7248):
		print('Scaffold must be less than 7248 for M13')
		print('</body>');
		print('</html>');
		exit();
elif (templType==2):
	fInStr='./lambda.txt'
	templName='Lambda'
	templAmount='0.5ng'
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
valid=''
rowCalls=getRows(scLen)
print("<h3>aPrime Output, Length "+str(scLen)+" nts</h3>");
print('<table>');
print('<tr><th>Forward Primer (1uM)</th><th>Reverse Primer (20nM)</th><th>Validated?</th></tr>');
for row in rowCalls:
	if row[0]==row[8]:
		valid='Yes'
	else:
		valid=''
	print('<tr><td>'+row[2]+'</td><td>'+row[3]+'</td><td>'+valid+'</td></tr>');
print('</table>');
print("""<br><hr><br>
<H3>Protocol</H3>
On ice, mix:<br>
1X Quantabio AccuStart HIFI buffer<br>
2mM MgSO<sub>4</sub><br>
1uM Forward primer<br>
15-20nM Reverse primer<br>
0.6ng/ul of reaction of M13 or gBlock template<br>
ddH<sub>2</sub>O to reaction volume<br>
0.004ul of enzyme/ul of reaction<br>
<br>
Thermocycle:<br>
1. 94C for 1 minute<br>
2. 94C for 30 seconds<br>
3. 55-57C for 30 seconds<br>
4. 68C for 1 minute 30 seconds per kb<br>
5. Goto step 2 30 times<br>
6. 4C<br>
<br>
Purify:<br>
Load reaction with 1X loading buffer to low-melt agarose gel<br>
Run until bands are separated. ssDNA band will run faster than the dsDNA, have a square look, and will be a pinker hue if stained with SybrSafe.<br>
Cut the band out with a clean razor<br>
Extract using ZymoClean gel cleanup kit.<br>""");
print('</body>');
print('</html>');
