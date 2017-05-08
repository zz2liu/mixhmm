""" try to model the lrr sds and means from given data

Ref: http://zunzun.com/
"""
import math
from py25 import set_trace
from numpy import array, arange, exp, log2
from matplotlib import pyplot as pp
from py_cna._sidcon import lrr_from_cn

cns = arange(5) + 0.01
sds = array(map(float, '1.33 .29 .22 .20 .19'.split()))
means = array(map(float, '-3.53 -0.66 0  0.40 0.68'.split()))

def plot_f(f, grid=arange(0.01, 8, 0.01)):
    pp.plot(grid, array([f(e) for e in grid]))

def plot_together(f, x=cns, y=sds):
    plot_f(f)
    pp.plot(x, y, 'r+')
    pp.show()

def Weibull2D_model(x_in):
    temp = 0.0

    # coefficients
    a = 1.8152603328611427E-01
    b = -1.1484743092476073E+00
    c = 2.3602773953818241E+00
    d = 5.1995417765774243E-01

    temp = a - b*math.exp(-1.0 * c * math.pow(x_in, d))
    return temp

def SigmoidB_WO_2D_model(x_in):
	temp = 0.0

	# coefficients
	a = 1.2303032185919699E+02
	b = -1.8598708649540483E+00
	c = -3.9731435906372736E-01
	d = 1.9999741514476532E-01

	temp = a / (1.0 + math.exp(-1.0 * (x_in - b) / c)) + d
	return temp

def offset_exp(x):
    a = 5.6841104810567913E+02
    b = -2.5077450119788356E+00
    c = -6.2206743262043007E+00
    d = 0.19 #1.9997206580639676E-01
    return a * exp(b*x + c) + d

def exp_decay_to_certain(x):
    a, b, c = 1.3299, 2.5078, 0.19, 0.1999
    return ( a - c ) * exp( - b * x ) + c

def QuadraticLogarithmic2D_model(x_in):
	temp = 0.0 

	# coefficients
	a = -6.5284961622742876E-01
	b = 8.8882428613115572E-01
	c = 5.7328236763747202E-02

	temp = a + b*math.log(x_in) + c*math.pow(math.log(x_in), 2.0)
	return temp

def main():
    f = lambda x: 1.507142857e-1*x**2 - 8.398571429e-1 * x + 1.221428571
    f = lambda x: exp(-1.149147781 * x**2)
    f = lambda x: (4.345574856e-1)**x
    f = offset_exp
    f = Weibull2D_model
    f = SigmoidB_WO_2D_model
    f = exp_decay_to_certain

    f = lambda x: log2(x/2)
    f = lrr_from_cn
    #f= QuadraticLogarithmic2D_model

    plot_together(exp_decay_to_certain, y=sds)
    plot_together(QuadraticLogarithmic2D_model, y=means)

    txt = '''\
    0 1.33
    1 .29
    2 .22
    3 .20
    4 .19
    '''
txt = '''
-3.53 , 0
-0.66 , 1
0     , 2
0.40  , 3
0.68  , 4
'''
txt = '''\
 0.01  ,-3.53
 1  ,-0.66
 2  ,0    
 3  ,0.40 
 4  ,0.68 
 '''

if __name__ == '__main__':
    main()

