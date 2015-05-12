
from ROOT import *
import array
import pandas as pd
import re

gROOT.Reset()

# Setting of general plotting style:
gStyle.SetCanvasColor(0)
gStyle.SetFillColor(0)

# Settings for statistics box:
gStyle.SetOptStat(1111)
gStyle.SetOptFit(1111)


charged = ['E','D','K','R','H'] 
peptide_bond = ['Q','P','N','A','T','S','V','G','M','C','I','L']
aromatic = ['Y','F','W']
polar = ['Q','N','T','S','M','C','W','H']
hydrophobic = ['A','I','L','F','V','P','G']

hp_l = dict([('I',4.5), ('F',2.8),('V',4.2),('L',3.8),('W',-0.9),('M',1.9),
	("A",1.8),('G',-0.4),('C',2.5),('Y',-1.3),('P',-1.6),('T',-0.7),('S',-0.8),
	('H',-3.2),('E',-3.5),('N',-3.5),('Q',-3.5),('D',-3.5),('K',-3.9),('R',-4.5)])

def sqr(x):
	return x*x

def print_full(x):
	import pandas as pd
	pd.set_option('display.max_rows', len(x))
	print(x)
	pd.reset_option('display.max_rows')
def read_protherm(input_csv):
	import pandas as pd

	#For graphical purposes
	#pd.set_option('display.line_width', 5000) 
	#pd.set_option('display.max_columns', 60) 

	import numpy as np
	#import matplotlib.pyplot as plt

	pd.options.display.mpl_style = 'default'

	'Read the data'
	data = pd.read_csv(input_csv, sep=';', decimal = ',', converters={'ddG': lambda x: float(x.replace(',','.'))})

	
	'Return the dataframe'
	return data
def df_list_import(path):
	import os
	import pandas as pd

	'Here we get a list of files that has the word data_frame in our working directory'
	lsdir = os.listdir(path)

	word = 'data_frame.csv'
	df_index = []
	for i in range(len(lsdir)):
		if word in lsdir[i]:
			df_index.append(lsdir[i])

	df_list = []
	for i in range(len(df_index)):
		df = read_protherm(path+df_index[i])
		df_list.append(df)

	#df_list = df_list[pd.notnull(df_list['ddG'])] #Remove rows with null data in ddG

	return pd.concat(df_list).reset_index()
def only_one():
	import re
	
	#df = read_protherm("df_output/1STN_data_frame.csv")
	df = df_list_import("df_output/")

	#return print_full(df['Mut_AA']), print_full(df['WT_AA'])

	'Here we remove all double and multiple mutations'
	df = df[df['Mut_AA'].map(len) < 2]
	df.reset_index()

	'Not really necessary - just in case WILD is not removed from above'
	df = df[df.Mutation != 'WILD']

	return df
def choose_amino(letter):
	df = only_one()
	grouped = df.groupby('Mut_AA')
	ddG_group = array.array('f')
	dlogp_group = array.array('f')

	for j in range(len(letter)):
		c_letter = letter[j]
		aminoacid = grouped.groups[c_letter]
		for i in range(len(aminoacid)):
			index = aminoacid[i]
			ddG_group.append(float(df['ddG'][index]))
			dlogp_group.append(float(df['d_log(p)'][index]))

	return ddG_group, dlogp_group, len(dlogp_group)
def amino_change(original, target):
	'This function returns a data frame'
	df = only_one()
	org_reg = "(" + ")|(".join(original) + ")"
	org_tar = "(" + ")|(".join(target) + ")"
	
	list_df = []
	for index in df.index.tolist():

		if re.match(org_reg, df['WT_AA'][index]):

			list_df.append(df.ix[[index]])

	org_df = pd.concat(list_df)

	list_org_df = []
	for index2 in org_df.index.tolist():

		if re.match(org_tar, org_df['Mut_AA'][index2]):

			list_org_df.append(org_df.ix[[index2]])

	df_org_tar = pd.concat(list_org_df)


	return df_org_tar	
def residuals(p0,x,y):
	res_list = []
	for i in range(len(x)):
		res_list.append(p0*x[i]-y[i])

	return res_list
def hydro_change(hp_l):
	df = only_one()
	dHB = []
	for index in df.index.tolist():
		X = df['WT_AA'][index]
		Y = df['Mut_AA'][index]
		dHB.append(hp_l[X]-hp_l[Y])

	df['dHB'] = dHB
	
	return df

def runs_test():
	r = residuals(p0,x,y)
	N_negative = 0.0
	N_positive = 0.0
	runs_start = 0
	runs = 1

	for i in range(len(r)):
		if r[i] < 0:
			if runs_start == 1:
				pass
			else:
				runs += 1
			N_positive += 1.0
			runs_start = 1

		else:
			if runs_start == 2:
				pass
			else:
				runs += 1

			N_negative += 1.0
			runs_start = 2

	N_total = N_positive+N_negative
	mean_Wald_Wolf = ((2.0*N_positive*N_negative)/(N_total))+1.0
	Variance_Wald = 2.0*N_positive*N_negative*(2.0*N_positive*N_negative-N_total)/(sqr(N_total)*(N_total-1.0))
	sigma = sqrt(Variance_Wald)

	test = (mean_Wald_Wolf-runs)/sigma

	return mean_Wald_Wolf,sigma,test







df = only_one()
y = array.array('f',list(df['d_log(p)'].astype(float)))
x = array.array('f',list(df['ddG'].astype(float)))
Npoints = len(x)



fitting = '[0]*x'
g = TGraph(Npoints,x,y)
fit = TF1("fit",fitting,-8,3)
g.Fit("fit","C")


p0 = fit.GetParameter(0)


####### residuals vs delta Hydrobhobic #######

'here we lend some parameters from the first fit'
res = residuals(p0,x,y)

df_h = hydro_change(hp_l)

r = array.array('f',res)
x_h = array.array('f',list(df_h['dHB'].astype(float)))
Npoints_h = len(x_h)

'Ok, let us plot it'

c0 = TCanvas("c0","",50,50,800,500)
#c0.Divide(3,4)

##############
g_h = TGraph(Npoints_h,x_h,r)
fit2 = TF1("fit2",'[0]*x+[1]',-8,10)
g_h.Fit("fit2","C")

g_h.SetTitle('Data. Size, N: %s.'%Npoints_h)
g_h.SetLineWidth(1)
g_h.SetMarkerStyle(20)
g_h.SetMarkerSize(0.2)
g_h.GetYaxis().SetTitle("residuals")
g_h.GetXaxis().SetTitle("change in hydrophobicity")

#c0.cd(1)

g_h.Draw("AP")
fit2.Draw('same')


c0.Update()
c0.SaveAs("residual-hydro-fit.ps")

print runs_test()

	
raw_input( ' ... Press enter to exit ... ' )
