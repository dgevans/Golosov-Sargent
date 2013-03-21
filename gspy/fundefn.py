#FUNDEFN Defines a function family structure
#This is the file that mimics fundefn from the compecon files in matlab
#Inputs
#   bastype  : the string referencing function family ('cheb','spli' or 'lin')
#   n        : order of approximation along each dimension
#   a        : left endpoints of interpolation intervals in each dimension
#   b        : right endpoints of interpolation intervals in each dimension
#   order    : for 'spli' bastype, the order of the spline (default: 3 for cubic)
#   s1,s2... : additional column vectors for appending discrete variables to
#                the basis
#
# OUTPUT
#   fspace  : a function family structure
#
# USES: fundef
# See also: FUNDEF, FUNEVAL, FUNBAS.
#Original code done by Paul L. Fackler and Mario J. Miranda

def fundefn(bastype, n, a, b, order, varargin):
	if bastype.any():
		d = 0.
		params = np.empty(1, varargin.size)
	else:
		d = n.size
		if a.size != d:
			raise ValueError('a must be same dimension as n')
		if b.size != d:
			raise ValueError('b must be the same dimension as n')
		if any(a>b):
			raise ValueError('left endpoint must be less than right endpoint')
		if any(n<2):
			raise ValueError('n(i) must be greater than 1')
		#if nargin<5 and len(order):
			#order = 3
		