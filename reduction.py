# -*- coding: UTF-8 -*-

import matplotlib.pyplot as plt
import scipy
import scipy.stats as st
import math
import numpy as np

import pdb

#
PI = 3.1415926
g = 9.81
dt = 0.08

# curshing and tensile intensity of concrete and steel bar, MPa
Ints = {'Concrete':[26.8, 2.39], 'Steel':[540.0, 540.0]}

# fit distribution
def dist_fit(data, dist_name):
	#dist_names = ["beta", "lognorm", "gamma", "expon"]
	dist = getattr(st, dist_name)
	param = dist.fit(data)
	return param

# upward load
def upward_load(p_RPV):
	Rv = 2.387       	# radius of RPV,m
	Mv = 418e3       	# RPV mass, kg
	Av = PI*Rv**2    	# equivalent area of PRV cross section, m2
	K = 2083.2e6/1.24	# F=kH

	Vv = []
	Fv = []
	Hv = []

	for i in range(len(p_RPV)):
		Iv = max(0.0, p_RPV[i]*1e6*Av-Mv*g)*dt
		Vv.append(Iv/Mv)
		Fv.append(Iv/dt/1e6)              # MN
		Hv.append((math.sqrt((Mv*g)**2+Mv*K*Vv[i]**2)-Mv*g)/K)

	return Fv,Hv

# plot histogram and fitting distribution curve
def plot_hist(data, label_id, dfit_info, num_bins=20):
	#
	xlabel_text = ("Maxium pressure at the side wall/MPa", 
		"Maximum impulsion per area at the side wall/MPa-s",
		"Coolant kinetic energy at the side wall/MJ",
		"Upward force/MN",
		"Upward displacement/m")

	fig = plt.figure()
	#
	freq, bins, patch = plt.hist(data, num_bins, density=True, rwidth=1.0, edgecolor='black')

	#
	#xmin = max(0.1, min(data))
	xmin = bins[1]/2.0
	xmax = max(data)
	xstep = (xmax-xmin)/1000.0 
	x = scipy.arange(xmin, xmax, xstep)

	dist_name = dfit_info[0]
	dist = getattr(st, dist_name)
	param = dfit_info[1]
	y_fit = dist.pdf(x, *param[:-2], loc=param[-2], scale=param[-1])
	plt.plot(x, y_fit, 'r--', linewidth=3, label=dist_name)

	plt.xlabel(xlabel_text[label_id])
	plt.ylabel("Probability")
	plt.legend(loc="upper right")

	#plt.show()

# upper bounding limit of load
def load_lim(id):
	#
	D_coldleg = (698.5 + 69*2)/1000.0
	R_pit = 2.77
	R_out = 5.0
	H_pit = 6.15

	Dia_bar = 0.040            # diameter of steel rebar, m
	N_rebar = 25               # Number of steel rebar per unit area ???
	A_unit = 1.0*1.8           # ???

	xi = 0.803                 # decay factor

	Ir = 93.0                  # impulsion of the reference case, MN-s
	Er = 1850                  # kinetic energy of the reference case, MJ

	d_hole = 1.33
	l = 0.8
	h = 0.5*d_hole

	#
	A = (R_out - R_pit)*H_pit  # area of the cross section of the pit wall (c.f. schematics) 



	CDIF,TDIF,DIF= get_DIF(p_wall, impmax)

	if id == 0:
		# maximum pressure, MPa
		Pc = Ints['Concrete'][0]
		Ps = Ints['Steel'][0]

		Sc = 2*PI*R_pit*H_pit
		Ss = N_rebar*Dia_bar*2*PI*R_pit*H_pit        # ???
		S = 2*PI*R_pit*H_pit

		Xlim = (CDIF*Pc*Sc + DIF*Ps*Ss)/S

	elif id == 1 or id == 2:
		#
		Tc = Ints['Concrete'][1]
		Ts = Ints['Steel'][1]

		sigma_s = (N_rebar*PI*Dia_bar**2/4.0)/A_unit
		sigma_c = 1.0 - sigma_s

		F = (TDIF*Tc*sigma_c + DIF*Ts*sigma_s)*A
		P = F/(xi*R_pit*H_pit)

		if id == 1:
			# maximum impulsion per area, MPa-s
			Xlim = P*dt
		else:
			# coolant kinetic energy
			Xlim = (P*dt*A)/Ir*Er

	elif id == 3:
		# upward force
		Pc = Ints['Concrete'][0]
		Ps = Ints['Steel'][0]

		Ss = (N_rebar*Dia_bar)*(l*h)/A_unit*d_hole
		Sc = l*h
		S = Sc

		PX = (Pc*Sc + Ps*Ss)/S

		Xlim = 6.0*PX*d_hole*l

	else:
		# upward displacement
		ratio = 0.2
		Xlim = ratio*D_coldleg

	return Xlim

# plot fragile curve together with the cumulative probability curve of loads
def plot_frag(X, id, dfit_info):
	#
	#pdb.set_trace()
	xlabel_text = ("Maxium pressure at the side wall/MPa", 
		"Maximum impulsion per area at the side wall/MPa-s",
		"Coolant kinetic energy at the side wall/MJ",
		"Upward force/MN",
		"Upward displacement/m")

	mu = 0.5*X
	sigma = X/(2.0*st.norm.ppf(0.95))

	xmin = 0.01
	#xmax = max(data)
	xmax = X*100.0
	xstep = (xmax-xmin)/10000.0 
	x = scipy.arange(xmin, xmax, xstep)
	dist_name = dfit_info[0]
	dist = getattr(st, dist_name)
	param = dfit_info[1]

	y_fragile = st.norm.cdf(x, loc=mu, scale=sigma)
	y_cdf = dist.cdf(x, *param[:-2], loc=param[-2], scale=param[-1])

	#
	fig = plt.figure()
	ax_fra = fig.add_subplot(111)
	ax_fra.semilogx(x, y_fragile, linewidth=3, label="Failure probability")
	ax_fra.set_ylim(0.0, 1.1)
	ax_fra.set_xlabel(xlabel_text[i])
	ax_fra.set_ylabel("Failure probability")

	ax_load = ax_fra.twinx()
	ax_load.semilogx(x, y_cdf, 'r--', linewidth=3, label="Cumulative load probability")
	ax_load.set_ylim(0.0, 1.1)
	ax_load.set_ylabel("Cumulative load probability")

	fig.legend(loc="upper left", bbox_to_anchor=(0,1), bbox_transform=ax_fra.transAxes)

	plt.grid()
	#plt.show()

# dynamic increase factor
def get_DIF(p_wall, impmax):
	# characteristic time scale
	epsilon = [x/y for x,y in zip(p_wall, impmax)]
	eps = np.mean(epsilon)
	
	fc = Ints['Concrete'][0]
	ft = Ints['Concrete'][1]
	fy = Ints['Steel'][0]

	alpha = 1.0/(5.0 + 0.9*fc)
	delta = 1.0/(10.0 + 0.6*ft)
	gamma = 10.0**(6.156*alpha - 2.0)
	beta = 10.0**(7.112*delta - 2.33)

	if eps <= 30.0:
		CDIF = (eps/3e-5)**(1.026*alpha)
		TDIF = (eps/3e-5)**(1.016*delta)
	else:
		CDIF = gamma*(eps/3e-5)**(1.0/3.0)
		TDIF = beta*(eps/3e-6)**(1.0/3.0)

	DIF = (eps/1e-4)**(0.019 - 0.009*fy/414.0)

	return CDIF,TDIF,DIF

# evaluate the failure probability for each load
def pfail_load(i, X, dfit_info):
	#
	# P(Fail) = sum{P(Fail|Load_xi)*P(Load_xi)}
	# P(Fail|Load_xi): fragile curve
	# P(Load_xi) = pdf*dx

	mu = 0.5*X
	sigma = X/(2.0*st.norm.ppf(0.95))

	xmin = 0.01
	#xmax = max(data)
	xmax = X + 4*mu
	xstep = (xmax-xmin)/10000.0 
	x = scipy.arange(xmin, xmax, xstep)
	dist_name = dfit_info[0]
	dist = getattr(st, dist_name)
	param = dfit_info[1]

	y_fragile = st.norm.cdf(x, loc=mu, scale=sigma)
	y_pdf = dist.pdf(x, *param[:-2], loc=param[-2], scale=param[-1])

	#if i == 2:
	#	pdb.set_trace()

	pf = np.dot(y_fragile, y_pdf)*xstep

	return pf

# synthesized/overall containment failure probability
def pf_tot(pf):
	#
	# pf_tot = 1.0 - product((1.0-pf_i))
	temp = 1.0
	for i in range(len(pf)):
		temp *= (1.0 - pf[i])

	pf_cont = 1.0 - temp

	return pf_cont

# statistic quantities
def stats_var(data):
	#
	xmin = min(data)
	xmax = max(data)
	xptp = np.ptp(data)              # value range (peak to peak)
	xmean = np.mean(data)
	xvar = np.var(data, ddof=1)      # sample variance = sum((x-xmean)^2)/(n-1)
	xstd = np.std(data, ddof=1)
	xmid = np.median(data)

	return xmin, xmax, xptp, xmean, xvar, xstd, xmid  




if __name__ == '__main__':
	#
	load_file = "load.dat"
	stat_file = "stat.dat"
	prob_file = "fail_prob.dat"
	load_dir = "D:\\Research\\CNPE\\KY1606_315\\SE\\"
	load_path = load_dir + load_file
	prob_path = load_dir + prob_file
	stat_path = load_dir + stat_file

	# lateral side wall load
	p_wall = []         # maximum pressure at the lateral side wall, MPa
	p_RPV = []          # maximum pressure at the lower head of RPV, MPa
	liq_ke = []         # liquid coolant kinetic energy, MJ
	impmax = []         # maximum impulsion at the lateral side wall per unit area, MPa-s

	# upward load
	Fv = []
	Hv = []

	# statistic quantities/variables results
	sv = []

	# Tpye and its associated parameters for the fitting distribution
	dist_names = ["gamma", "gamma", "gamma", "gamma", "gamma"]
	dfit_info = []

	#
	pf = []

	# read loads on the lateral side wall from the external file
	with open(load_path, 'r') as f:
		lines = f.readlines()
		for line in lines[1:]:
			data = [float(x) for x in line.strip().split()]
			p_wall.append(data[0])
			liq_ke.append(data[1])
			p_RPV.append(data[2])
			impmax.append(data[3])

	# evaluate upward loads
	Fv,Hv = upward_load(p_RPV)

	# loop for each type of loads
	for i, data in enumerate([p_wall, impmax, liq_ke, Fv, Hv]):
		#
		# calculate the statistic quantities for each load
		sv.append(stats_var(data))

		# fit the distribution for each load
		dist_name = dist_names[i]
		dfit_info.append((dist_name, dist_fit(data, dist_name)))

		# plot histogram along with its distribution fitting curve
		#plot_hist(data, i, dfit_info[i], 20)

		# plot fragile curve
		Xlim = load_lim(i)                  # upper bound of loads
		print Xlim
		#plot_frag(Xlim, i, dfit_info[i])

		# compute the overall containment failure probability
		pf.append(pfail_load(i, Xlim, dfit_info[i]))
		pf_cont = pf_tot(pf)


	# write statistic results to the external file
	with open(stat_file, 'w') as f:
		#
		header = ["xmin", "xmax", "xptp", "xmean", "xvar", "xstd", "xmid"]
		f.write("%20s"*7 % tuple(header))
		f.write("\n")
		for i in range(len(sv)):
			f.write("%20.5f"*7 % sv[i])
			f.write("\n")

	# write the failure probabilities to the external file
	with open(prob_path, 'w') as f:
		for i in range(len(pf)):
			f.write("%15.4f\n" % pf[i])

		f.write("pf_cont = %15.4f\n" % pf_cont)
