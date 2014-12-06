from string import *
from math import *
import numpy as np
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

f = open ('result.txt', 'r')

lines = f.readlines()

num = len(lines) - 3

p=[]; dotp=[]; x=[]; y=[]; z=[]; lum=[]; R=[];


for i in range (24, num):
	line = split(lines[i])
	try:
		p.append(float(line[0]))
		dotp.append(log10(abs(float(line[1]))))
		x.append(float(line[2]))
		y.append(float(line[3]))
		z.append(float(line[4]))
		R.append(sqrt(float(line[2])**2.0 + float(line[3])**2.0))
		lum.append(log10(float(line[8])))
	except:
		print line, 'is excluded'

standard = open('p_dotp.txt', 'r')

st_p=[]; st_dotp=[]; st_x=[]; st_y=[]; st_z=[]; st_lum=[]; st_R=[];

for lines in standard.readlines():
	line = split(lines)
	try:
		st_p.append(float(line[1]))
		st_dotp.append(log10(abs(float(line[2]))))
		st_z.append(float(line[3]))
		st_x.append(float(line[4]))
		st_y.append(float(line[5]))
		st_R.append(sqrt(float(line[4])**2.0 + float(line[5])**2.0))
	except:
		print line, 'is excluded'

	try:
		st_lum.append(log10(float(line[6])))
	except:
		print line, ' is excluded from luminosity analysis'


st_p1    = np.asarray(st_p)
st_dotp1 = np.asarray(st_dotp)
st_z1    = np.asarray(st_z)
st_R1    = np.asarray(st_R)
st_lum1  = np.asarray(st_lum)

st_p    = st_p1    [ st_dotp1 > -17.3  ] 
st_dotp	= st_dotp1 [ st_dotp1 > -17.3  ]
st_z    = st_z1    [ (st_dotp1 > -17.3) & (st_z1 != 1.76) & (st_z1 != -1.76) ]
st_R    = st_R1    [ st_dotp1 > -17.3  ]
#st_lum2  = st_lum1  [ st_dotp1 > -17.3  ]

#pbins = np.linspace(0.01, 5.0, 50)

#plt.hist(p, pbins, alpha=0.5, normed=True)
#plt.hist(st_p, pbins, alpha=0.5, normed=True)
#plt.xlabel('Periods (sec)')
#plt.ylabel('Relative number')
#plt.show()

#dotpbins = np.linspace(-20, -10, 50)

#plt.hist(dotp, dotpbins, alpha=0.5, normed=True)
#plt.hist(st_dotp, dotpbins, alpha=0.5, normed = True)
#plt.xlabel('Period derivatives (sec/sec)')
#plt.ylabel('Relative number')
#plt.show()

#bins = np.linspace(-2.0, 2.0, 50)

#plt.hist(z, bins, alpha=0.5, normed=True)
#plt.hist(st_z, bins, alpha=0.5, normed = True)
#plt.xlabel('Z distance (kpc)')
#plt.ylabel('Relative number')
#plt.show()

#bins = np.linspace(0, 20, 50)

#plt.hist(R, bins, alpha=0.5, normed = True)
#plt.hist(st_R, bins, alpha=0.5, normed=True)
#plt.xlabel('R (kpc)')
#plt.ylabel('Relative number')
#plt.show()

print 'K-S test for periods: ', ks_2samp(st_p, p)
print 'K-S test for period derivatives: ', ks_2samp(st_dotp, dotp)
print 'K-S test for z coordinate: ', ks_2samp(st_z, z)
print 'K-S test for R coordinate: ', ks_2samp(st_R, R)
print 'K-S test for Lum: ', ks_2samp(st_lum, lum)

ks1 = ks_2samp(st_p, p)
ks2 = ks_2samp(st_dotp, dotp)
ks3 = ks_2samp(st_z, z)
ks4 = ks_2samp(st_R, R)
ks5 = ks_2samp(st_lum, lum)

s1='K-S test P   : '+str('%5.3f' % ks1[0]) + ' i.e P=' + str('%6.2e' % ks1[1]) 
s2='K-S test dotP: '+str('%5.3f' % ks2[0]) + ' i.e P=' + str('%6.2e' % ks2[1])
s3='K-S test z   : '+str('%5.3f' % ks3[0]) + ' i.e P=' + str('%6.2e' % ks3[1])
s4='K-S test R   : '+str('%5.3f' % ks4[0]) + ' i.e P=' + str('%6.2e' % ks4[1])
s5='K-S test Lum : '+str('%5.3f' % ks5[0]) + ' i.e P=' + str('%6.2e' % ks5[1])



pbins = np.linspace   (0.01, 5.0, 50)
dotpbins = np.linspace(-20, -10, 50)
zbins = np.linspace(-2.0, 2.0, 50)
rbins = np.linspace(0, 20, 50)
lbins = np.linspace(-3, 4, 50)

font0 = FontProperties()
family = 'monospace'
font0.set_family(family)

fig=plt.figure(1, figsize=(12, 9), dpi=80)
#fig.text(.1, .1, 'Comparison')
#plt.rc('font', family='monospace')
plt.suptitle('Comparison of the PMBS and SMBS with synthesis.')
#plt.annotate('Comparison')

plt.subplot(321)
plt.hist(p, pbins, alpha=0.5, normed=True)
plt.hist(st_p, pbins, alpha=0.5, normed=True)
plt.xlabel('Period (sec)')
plt.ylabel('Relative number')

plt.subplot(322)
plt.hist(dotp, dotpbins, alpha=0.5, normed=True)
plt.hist(st_dotp, dotpbins, alpha=0.5, normed=True)
plt.xlabel('Period derivative (sec/sec)')
plt.ylabel('Relative number')

plt.subplot(323)
plt.hist(z, zbins, alpha=0.5, normed=True)
plt.hist(st_z, zbins, alpha=0.5, normed = True)
plt.xlabel('Z distance (kpc)')
plt.ylabel('Relative number')

plt.subplot(324)
plt.hist(R, rbins, alpha=0.5, normed = True)
plt.hist(st_R, rbins, alpha=0.5, normed=True)
plt.xlabel('R (kpc)')
plt.ylabel('Relative number')

plt.subplot(325)
plt.hist(lum, lbins, alpha=0.5, normed=True)
plt.hist(st_lum, lbins, alpha=0.5, normed=True)
plt.xlabel('Pseudo-radio luminosity (mJy kpc^2)')
plt.ylabel('Relative number')

plt.subplot(326)
plt.text(0.05, 0.75, s1, fontproperties=font0)
plt.text(0.05, 0.6,  s2, fontproperties=font0)
plt.text(0.05, 0.45, s3, fontproperties=font0)
plt.text(0.05, 0.30, s4, fontproperties=font0)
plt.text(0.05, 0.15, s5, fontproperties=font0)

plt.savefig('compar.pdf', format='pdf')
plt.show()



