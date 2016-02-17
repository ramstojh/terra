#This program computes the convective mass of any solar twin making a double interpolation between their mass and metallicity.
#The mass cover a range from 0.5 to 1.3
#The metallicity cover a range from -0.127 to 0.474

#----------------------------------------------------------------------------------
import numpy

#----------------------------------------------------------------------------------
data = numpy.loadtxt('./TABLES/model-lionel.txt')

Mx = float(raw_input('Mass of the solar twin: '))
Zx = float(raw_input('Metallicity [X/H]     : '))

for i in range(len(data[0])):
	if Zx == data[0][i]:
		Z1 = Zx
		Z2 = 0.0
		xp = i
		break
	elif data[0][i-1] < Zx < data[0][i]:
		Z1 = data[0][i-1]
		Z2 = data[0][i]
		xp = i

for j in range(10):
	if Mx == data[j][0]:
		M1 = Mx
		M2 = 0.0
		yp = j
		break
	elif data[j-1][0] < Mx < data[j][0]:
		M1 = data[j-1][0]
		M2 = data[j][0]
		yp = j

m1 = (Z2 - Zx)/(Z2 - Z1)
m2 = (Zx - Z1)/(Z2 - Z1)
b1 = (M2 - Mx)/(M2 - M1)
b2 = (Mx - M1)/(M2 - M1)

if Z1 == Zx and M1 == Mx:
	CM = data[yp][xp]
elif Z1 == Zx and M1 == data[yp-1][0]:
	CM11 = data[yp-1][xp]
	CM12 = 0
	CM21 = data[yp][xp]
	CM22 = 0
	CM = (m1*CM11 + m2*CM12)*b1 + (m1*CM21 + m2*CM22)*b2
elif Z1 == data[0][xp-1] and M1 == Mx:
	CM11 = data[yp][xp-1]
	CM12 = data[yp][xp]
	CM21 = 0
	CM22 = 0
	CM = (m1*CM11 + m2*CM12)*b1 + (m1*CM21 + m2*CM22)*b2
elif Z1 == data[0][xp-1] and M1 == data[yp-1][0]:
	CM11 = data[yp-1][xp-1]
	CM12 = data[yp-1][xp]
	CM21 = data[yp][xp-1]
	CM22 = data[yp][xp]
	CM = (m1*CM11 + m2*CM12)*b1 + (m1*CM21 + m2*CM22)*b2
print '====================================='
print 'The convective mass is:', '%0.3f' %CM
print '====================================='
