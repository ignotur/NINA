from string import *
from math import *
import numpy as np
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

f = open ('result.txt', 'r')

lines = f.readlines()

num = len(lines) - 3

p=[]; dotp=[]; x=[]; y=[]; z=[]; lum=[]; R=[]; l=[];


for i in range (24, num):
	line = lines[i].split()
	try:
		p.append(float(line[0]))
		dotp.append(log10(abs(float(line[1]))))
		x.append(float(line[2]))
		y.append(float(line[3]))
		z.append(float(line[4]))
		R.append(sqrt(float(line[2])**2.0 + (float(line[3])-8.5)**2.0))
		lum.append(log10(float(line[8])))

		c_x = float(line[2])
		c_y = float(line[3])
		a1 = c_x
		a2 = c_y - 8.5
		b1 = 0
		b2 = -8.5
	
#		print 'look here -- ', (a1*b2 - b1*a2)/(a1*b1 + a2*b2)
	
		value = degrees(atan2 ( (a1*b2 - b1*a2), (a1*b1 + a2*b2) ) )

#		print 'look here -- ', value

		if (value < 0):
			value = value + 360.0
		l.append( value )
	except:
		print(line, 'is excluded')

standard = open('p_dotp.txt', 'r')

st_p=[]; st_dotp=[]; st_x=[]; st_y=[]; st_z=[]; st_lum=[]; st_R=[]; st_l=[]

for lines in standard.readlines():
	line = lines.split()
	try:
		st_p.append(float(line[1]))
		st_dotp.append(log10(abs(float(line[2]))))
		st_z.append(float(line[3]))
		st_x.append(float(line[4]))
		st_y.append(float(line[5]))
		st_R.append(sqrt(float(line[4])**2.0 + float(line[5])**2.0))
		c_x = float(line[4])
		c_y = float(line[5])
		a1 = c_x 
		a2 = c_y - 8.5
		b1 = 0.0
		b2 = -8.5

		value = degrees(atan2 ( (a1*b2 - b1*a2), (a1*b1 + a2*b2) ) )
		if (value < 0):
			value = value + 360.0
		st_l.append( value )

	except:
		print(line, 'is excluded')

	try:
		st_lum.append(log10(float(line[6])))
	except:
		print(line, ' is excluded from luminosity analysis')


st_p1    = np.asarray(st_p)
st_dotp1 = np.asarray(st_dotp)
st_z1    = np.asarray(st_z)
st_R1    = np.asarray(st_R)
st_lum1  = np.asarray(st_lum)
st_l1    = np.asarray(st_l)

st_p    = st_p1    [ st_dotp1 > -17.3  ] 
st_dotp	= st_dotp1 [ st_dotp1 > -17.3  ]
st_z    = st_z1    [ (st_dotp1 > -17.3) & (st_z1 != 1.76) & (st_z1 != -1.76) ]
st_R    = st_R1    [ (st_dotp1 > -17.3) & (st_z1 != 1.76) & (st_z1 != -1.76) ]
#st_l    = st_l1    [ st_dotp1 > -17.3  ]

z_new = np.asarray(z)
R_new = np.asarray(R)

b_list = np.degrees (z_new / R_new)  

st_z_new = np.asarray(st_z)
st_R_new = np.asarray(st_R)

st_b_list = np.degrees (st_z_new / st_R_new)

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

print('K-S test for periods: ', ks_2samp(st_p, p))
print('K-S test for period derivatives: ', ks_2samp(st_dotp, dotp))
print('K-S test for z coordinate: ', ks_2samp(b_list, st_b_list))
print('K-S test for R coordinate: ', ks_2samp(st_R, R))
print('K-S test for Lum: ', ks_2samp(st_lum, lum))

ks1 = ks_2samp(st_p, p)
ks2 = ks_2samp(st_dotp, dotp)
ks3 = ks_2samp(st_z, z)
ks4 = ks_2samp(st_l, l)
ks5 = ks_2samp(st_lum, lum)

s1='K-S test P   : '+str('%5.3f' % ks1[0]) + ' i.e P=' + str('%6.2e' % ks1[1]) 
s2='K-S test dotP: '+str('%5.3f' % ks2[0]) + ' i.e P=' + str('%6.2e' % ks2[1])
s3='K-S test b   : '+str('%5.3f' % ks3[0]) + ' i.e P=' + str('%6.2e' % ks3[1])
s4='K-S test l   : '+str('%5.3f' % ks4[0]) + ' i.e P=' + str('%6.2e' % ks4[1])
s5='K-S test Lum : '+str('%5.3f' % ks5[0]) + ' i.e P=' + str('%6.2e' % ks5[1])



pbins = np.linspace   (0.01, 5.0, 50)
dotpbins = np.linspace(-20, -10, 50)
bbins = np.linspace(-17.0, 17.0, 50)
rbins = np.linspace(0, 360, 50)
lbins = np.linspace(-3, 4, 50)

font0 = FontProperties()
family = 'monospace'
font0.set_family(family)

fig=plt.figure(1, figsize=(12, 13), dpi=80)
#fig.text(.1, .1, 'Comparison')
#plt.rc('font', family='monospace')
plt.suptitle('Comparison of the PMBS and SMBS with synthesis.')
#plt.annotate('Comparison')

plt.subplot(421)
plt.hist(p, pbins, alpha=0.5, normed=True)
plt.hist(st_p, pbins, alpha=0.5, normed=True)
plt.xlabel('Period (sec)')
plt.ylabel('Relative number')

plt.subplot(422)
plt.hist(dotp, dotpbins, alpha=0.5, normed=True)
plt.hist(st_dotp, dotpbins, alpha=0.5, normed=True)
plt.xlabel('Period derivative (sec/sec)')
plt.ylabel('Relative number')

plt.subplot(423)
plt.hist(b_list, bbins, alpha=0.5, normed=True)
plt.hist(st_b_list, bbins, alpha=0.5, normed = True)
plt.xlabel('Galactic latitude (deg)')
plt.ylabel('Relative number')

plt.subplot(424)
plt.hist(l, rbins, alpha=0.5, normed = True)
plt.hist(st_l, rbins, alpha=0.5, normed=True)
plt.xlabel('Galactic longitude (deg)')
plt.ylabel('Relative number')

plt.subplot(425)
plt.hist(lum, lbins, alpha=0.5, normed=True)
plt.hist(st_lum, lbins, alpha=0.5, normed=True)
plt.xlabel('Pseudo-radio luminosity (mJy kpc^2)')
plt.ylabel('Relative number')

plt.subplot(426)
plt.text(0.05, 0.75, s1, fontproperties=font0)
plt.text(0.05, 0.6,  s2, fontproperties=font0)
plt.text(0.05, 0.45, s3, fontproperties=font0)
plt.text(0.05, 0.30, s4, fontproperties=font0)
plt.text(0.05, 0.15, s5, fontproperties=font0)

logp    = np.log10(p)

plt.subplot(427)
plt.scatter(logp, dotp, s=1)
plt.xlabel('Log P (s) (Synthesis)')
plt.ylabel('Log dotP (s/s)')
plt.xlim([-1.5, 1.0])
plt.ylim([-18, -11])

plt.subplot(428)
plt.scatter(np.log10(st_p), st_dotp, s=1)
plt.xlabel('Log P (s) (PMBS & SMBS)')
plt.ylabel('Log dotP (s/s)')
plt.xlim([-1.5, 1.0])
plt.ylim([-18, -11])

plt.savefig('compar.pdf', format='pdf')
plt.show()



