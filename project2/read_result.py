

from ROOT import *
import array
import pandas as pd
import re
import numpy as np

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
secondary = ['H','B','E','G','I','T','S','-']

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
		try:
			aminoacid = grouped.groups[c_letter]
			for i in range(len(aminoacid)):
				index = aminoacid[i]
				ddG_group.append(float(df['ddG'][index]))
				dlogp_group.append(float(df['d_log(p)'][index]))
		except KeyError:
			pass

	return ddG_group, dlogp_group, len(dlogp_group)

def secondary_structure(letter):
	df = only_one()
	grouped = df.groupby('dssp')
	ddG_group = array.array('f')
	dlogp_group = array.array('f')

	for j in range(len(letter)):
		c_letter = letter[j]
		try:
			aminoacid = grouped.groups[c_letter]
			for i in range(len(aminoacid)):
				index = aminoacid[i]
				ddG_group.append(float(df['ddG'][index]))
				dlogp_group.append(float(df['d_log(p)'][index]))
		except KeyError:
			pass

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

def shuffle_data():
	import array
	import random
	df = df_list_import("df_output/")
	all_data = array.array('f',list(df['ddG'].astype(float)))
	random.shuffle(all_data)
	
	k = 0
	list1 = array.array('f')
	list2 = array.array('f')
	for i in range(len(all_data)):
		if k == 0:
			list1.append(all_data[i])
			k = 1
		else:
			list2.append(all_data[i])
			k = 0

	return list1, list2






fitting = "[0]*x"

df = only_one()
y = array.array('f',list(df['d_log(p)'].astype(float)))
x = array.array('f',list(df['ddG'].astype(float)))
Npoints = len(x)



print "The data consist of %s points."%Npoints



save_name = 'Giant_Plot'

####### HERE STARTS PLOTTING ########

c0 = TCanvas("c0","",50,50,1400,900)
c0.Divide(2,2)

##############
g = TGraph(Npoints,x,y)
fit = TF1("fit",fitting,-8,3)
g.Fit("fit","C")


'Correlation'
corr = g.GetCorrelationFactor()
z = 0.5*np.log((1+corr)/(1-corr))
z_std = 1/sqrt(Npoints-3)
cor_text = " Correlation, z = %s +- %s"%(z,z_std)

this_plot = 'Figure S1) All_data'

g.SetTitle('%s. Size, N: %s.'%(this_plot,Npoints)+cor_text)
g.SetLineWidth(1)
g.SetMarkerStyle(20)
g.SetMarkerSize(0.2)
g.GetYaxis().SetTitle("dlog(p)")
g.GetXaxis().SetTitle("ddG")

c0.cd(1)

g.Draw("AP")
fit.Draw('same')



#############

ddG_ch = choose_amino(charged)[0]
dlogp_ch = choose_amino(charged)[1]
N_ch = choose_amino(charged)[2]




g_ch = TGraph(N_ch,ddG_ch,dlogp_ch)
fit3 = TF1("fit3",fitting,-8,3)
g_ch.Fit("fit3","")

corr = g_ch.GetCorrelationFactor()
z = 0.5*np.log((1+corr)/(1-corr))
z_std = 1/sqrt(N_ch-3)
cor_text = " Correlation, z = %s +- %s"%(z,z_std)

g_ch.SetTitle('Figure S2) charged. Size, N: %s. %s'%(N_ch,cor_text))
g_ch.SetLineWidth(1)
g_ch.SetMarkerStyle(20)
g_ch.SetMarkerSize(0.2)
g_ch.GetYaxis().SetTitle("dlog(p)")
g_ch.GetXaxis().SetTitle("ddG")

c0.cd(2)

g_ch.Draw("AP")
fit3.Draw('same')

#############

ddG_pb = choose_amino(peptide_bond)[0]
dlogp_pb = choose_amino(peptide_bond)[1]
N_pb = choose_amino(peptide_bond)[2]

g_pb = TGraph(N_pb,ddG_pb,dlogp_pb)
fit4 = TF1("fit4",fitting,-8,3)
g_pb.Fit("fit4","")

corr = g_pb.GetCorrelationFactor()
z = 0.5*np.log((1+corr)/(1-corr))
z_std = 1/sqrt(N_pb-3)
cor_text = " Correlation, z = %s +- %s"%(z,z_std)

g_pb.SetTitle('Figure S3) peptide bond mutation. Size, N: %s. %s'%(N_pb,cor_text))
g_pb.SetLineWidth(1)
g_pb.SetMarkerStyle(20)
g_pb.SetMarkerSize(0.2)
g_pb.GetYaxis().SetTitle("dlog(p)")
g_pb.GetXaxis().SetTitle("ddG")

c0.cd(3)

g_pb.Draw("AP")
fit4.Draw('same')

#############

ddG_am = choose_amino(aromatic)[0]
dlogp_am = choose_amino(aromatic)[1]
N_am = choose_amino(aromatic)[2]

g_am = TGraph(N_am,ddG_am,dlogp_am)
fit5 = TF1("fit5",fitting,-8,3)
g_am.Fit("fit5","")

corr = g_am.GetCorrelationFactor()
z = 0.5*np.log((1+corr)/(1-corr))
z_std = 1/sqrt(N_am-3)
cor_text = " Correlation, z = %s +- %s"%(z,z_std)

g_am.SetTitle('Figure S4 )Aromatic mutation. Size, N: %s. %s'%(N_am,cor_text))
g_am.SetLineWidth(1)
g_am.SetMarkerStyle(20)
g_am.SetMarkerSize(0.2)
g_am.GetYaxis().SetTitle("dlog(p)")
g_am.GetXaxis().SetTitle("ddG")

c0.cd(4)

g_am.Draw("AP")
fit5.Draw('same')

#############

c3 = TCanvas("c3","",50,50,1400,900)
c3.Divide(2,2)

cp_df = amino_change(charged,peptide_bond)

dlogp_cp = array.array('f',list(cp_df['d_log(p)'].astype(float)))
ddG_cp = array.array('f',list(cp_df['ddG'].astype(float)))
N_cp = len(ddG_cp)

g_cp = TGraph(N_cp,ddG_cp,dlogp_cp)
fit7 = TF1("fit7",fitting,-8,3)
g_cp.Fit("fit7","")



corr = g_cp.GetCorrelationFactor()
z = 0.5*np.log((1+corr)/(1-corr))
z_std = 1/sqrt(N_cp-3)
cor_text = " Correlation, z = %s +- %s"%(z,z_std)

g_cp.SetTitle('Figure S5) Mutation: "peptide to charged amino acid"'+'. Size, N: %s. %s'%(N_cp,cor_text))
g_cp.SetLineWidth(1)
g_cp.SetMarkerStyle(20)
g_cp.SetMarkerSize(0.2)
g_cp.GetYaxis().SetTitle("dlog(p)")
g_cp.GetXaxis().SetTitle("ddG")

c3.cd(1)

g_cp.Draw("AP")
fit7.Draw('same')

#############

spec_df = amino_change(peptide_bond,charged)

dlogp_sp = array.array('f',list(spec_df['d_log(p)'].astype(float)))
ddG_sp = array.array('f',list(spec_df['ddG'].astype(float)))
N_sp = len(ddG_sp)

g_sp = TGraph(N_sp,ddG_sp,dlogp_sp)
fit6 = TF1("fit6",fitting,-8,3)
g_sp.Fit("fit6","")


corr = g_sp.GetCorrelationFactor()
z = 0.5*np.log((1+corr)/(1-corr))
z_std = 1/sqrt(N_sp-3)
cor_text = " Correlation, z = %s +- %s"%(z,z_std)

g_sp.SetTitle('Figure S6) Mutation: "charged to peptide amino acid"'+'. Size, N: %s. %s'%(N_sp,cor_text))
g_sp.SetLineWidth(1)
g_sp.SetMarkerStyle(20)
g_sp.SetMarkerSize(0.2)
g_sp.GetYaxis().SetTitle("dlog(p)")
g_sp.GetXaxis().SetTitle("ddG")

c3.cd(2)

g_sp.Draw("AP")
fit6.Draw('same')


#############

ddG_hp = choose_amino(hydrophobic)[0]
dlogp_hp = choose_amino(hydrophobic)[1]
N_hp = choose_amino(hydrophobic)[2]

g_hp = TGraph(N_hp,ddG_hp,dlogp_hp)
fit9 = TF1("fit9",fitting,-8,3)
g_hp.Fit("fit9","")


corr = g_hp.GetCorrelationFactor()
z = 0.5*np.log((1+corr)/(1-corr))
z_std = 1/sqrt(N_hp-3)
cor_text = " Correlation, z = %s +- %s"%(z,z_std)

g_hp.SetTitle('Figure S7) hydrophobic. Size, N: %s. %s'%(N_hp,cor_text))
g_hp.SetLineWidth(1)
g_hp.SetMarkerStyle(20)
g_hp.SetMarkerSize(0.2)
g_hp.GetYaxis().SetTitle("dlog(p)")
g_hp.GetXaxis().SetTitle("ddG")

c3.cd(3)

g_hp.Draw("AP")
fit9.Draw('same')

#############

ddG_p = choose_amino(polar)[0]
dlogp_p = choose_amino(polar)[1]
N_p = choose_amino(polar)[2]

g_p = TGraph(N_p,ddG_p,dlogp_p)
fit10 = TF1("fit10",fitting,-8,3)
g_p.Fit("fit10","")


corr = g_p.GetCorrelationFactor()
z = 0.5*np.log((1+corr)/(1-corr))
z_std = 1/sqrt(N_p-3)
cor_text = " Correlation, z = %s +- %s"%(z,z_std)

g_p.SetTitle('Figure S8) polar. Size, N: %s. %s'%(N_p,cor_text))
g_p.SetLineWidth(1)
g_p.SetMarkerStyle(20)
g_p.SetMarkerSize(0.2)
g_p.GetYaxis().SetTitle("dlog(p)")
g_p.GetXaxis().SetTitle("ddG")

c3.cd(4)

g_p.Draw("AP")
fit10.Draw('same')

#############
c4 = TCanvas("c4","",50,50,1400,900)
c4.Divide(2,2)

h_p_df = amino_change(hydrophobic,polar)

dlogp_h_p = array.array('f',list(h_p_df['d_log(p)'].astype(float)))
ddG_h_p = array.array('f',list(h_p_df['ddG'].astype(float)))
N_h_p = len(ddG_h_p)

g_h_p = TGraph(N_h_p,ddG_h_p,dlogp_h_p)
fit11 = TF1("fit11",fitting,-8,3)
g_h_p.Fit("fit11","")


corr = g_h_p.GetCorrelationFactor()
z = 0.5*np.log((1+corr)/(1-corr))
z_std = 1/sqrt(N_h_p-3)
cor_text = " Correlation, z = %s +- %s"%(z,z_std)

g_h_p.SetTitle('Figure S9) Mutation: "hydrophobic to polar amino acid"'+'. Size, N: %s. %s'%(N_h_p,cor_text))
g_h_p.SetLineWidth(1)
g_h_p.SetMarkerStyle(20)
g_h_p.SetMarkerSize(0.2)
g_h_p.GetYaxis().SetTitle("dlog(p)")
g_h_p.GetXaxis().SetTitle("ddG")

c4.cd(1)

g_h_p.Draw("AP")
fit11.Draw('same')

#############

p_h_df = amino_change(polar,hydrophobic)

dlogp_p_h = array.array('f',list(p_h_df['d_log(p)'].astype(float)))
ddG_p_h = array.array('f',list(p_h_df['ddG'].astype(float)))
N_p_h = len(ddG_p_h)

g_p_h = TGraph(N_p_h,ddG_p_h,dlogp_p_h)
fit12 = TF1("fit12",fitting,-8,3)
g_p_h.Fit("fit12","")


corr = g_p_h.GetCorrelationFactor()
z = 0.5*np.log((1+corr)/(1-corr))
z_std = 1/sqrt(N_p_h-3)
cor_text = " Correlation, z = %s +- %s"%(z,z_std)

g_p_h.SetTitle('Figure S10) Mutation: "polar to hydrophobic amino acid"'+'. Size, N: %s. %s'%(N_p_h,cor_text))

g_p_h.SetLineWidth(1)
g_p_h.SetMarkerStyle(20)
g_p_h.SetMarkerSize(0.2)
g_p_h.GetYaxis().SetTitle("dlog(p)")
g_p_h.GetXaxis().SetTitle("ddG")

c4.cd(2)

g_p_h.Draw("AP")
fit12.Draw('same')



#############
c1 = TCanvas("c1","",50,50,1400,900)
c1.Divide(2,3)



#################
ddG_hp = secondary_structure(['H'])[0]
dlogp_hp = secondary_structure(['H'])[1]
N_hp = secondary_structure(['H'])[2]

g_dsspH = TGraph(N_hp,ddG_hp,dlogp_hp)
fit13 = TF1("fit13",fitting,-8,3)
g_dsspH.Fit("fit13","")


corr = g_dsspH.GetCorrelationFactor()
z = 0.5*np.log((1+corr)/(1-corr))
z_std = 1/sqrt(N_hp-3)
cor_text = " Correlation, z = %s +- %s"%(z,z_std)

g_dsspH.SetTitle('Figure S11) Helix. Size, N: %s. %s'%(N_hp,cor_text))
g_dsspH.SetLineWidth(1)
g_dsspH.SetMarkerStyle(20)
g_dsspH.SetMarkerSize(0.2)
g_dsspH.GetYaxis().SetTitle("dlog(p)")
g_dsspH.GetXaxis().SetTitle("ddG")


c4.cd(3)

g_dsspH.Draw("AP")
fit13.Draw('same')

#################


#################
ddG_hp = secondary_structure(['T'])[0]
dlogp_hp = secondary_structure(['T'])[1]
N_hp = secondary_structure(['T'])[2]

g_dsspB = TGraph(N_hp,ddG_hp,dlogp_hp)
fit14 = TF1("fit14",fitting,-8,3)
g_dsspB.Fit("fit14","")


corr = g_dsspB.GetCorrelationFactor()
z = 0.5*np.log((1+corr)/(1-corr))
z_std = 1/sqrt(N_hp-3)
cor_text = " Correlation, z = %s +- %s"%(z,z_std)

g_dsspB.SetTitle('Figure S13) Turn. Size, N: %s. %s'%(N_hp,cor_text))
g_dsspB.SetLineWidth(1)
g_dsspB.SetMarkerStyle(20)
g_dsspB.SetMarkerSize(0.2)
g_dsspB.GetYaxis().SetTitle("dlog(p)")
g_dsspB.GetXaxis().SetTitle("ddG")


c4.cd(4)

g_dsspB.Draw("AP")
fit14.Draw('same')

#################
ddG_hp = secondary_structure(['E'])[0]
dlogp_hp = secondary_structure(['E'])[1]
N_hp = secondary_structure(['E'])[2]

g_dsspE = TGraph(N_hp,ddG_hp,dlogp_hp)
fit15 = TF1("fit15",fitting,-8,3)
g_dsspE.Fit("fit15","")


corr = g_dsspE.GetCorrelationFactor()
z = 0.5*np.log((1+corr)/(1-corr))
z_std = 1/sqrt(N_hp-3)
cor_text = " Correlation, z = %s +- %s"%(z,z_std)

g_dsspE.SetTitle('Figure S14) Strand. Size, N: %s. %s'%(N_hp,cor_text))
g_dsspE.SetLineWidth(1)
g_dsspE.SetMarkerStyle(20)
g_dsspE.SetMarkerSize(0.2)
g_dsspE.GetYaxis().SetTitle("dlog(p)")
g_dsspE.GetXaxis().SetTitle("ddG")


c1.cd(1)

g_dsspE.Draw("AP")
fit15.Draw('same')

#################
ddG_hp = secondary_structure(['-'])[0]
dlogp_hp = secondary_structure(['-'])[1]
N_hp = secondary_structure(['-'])[2]

g_dsspT = TGraph(N_hp,ddG_hp,dlogp_hp)
fit16 = TF1("fit16",fitting,-8,3)
g_dsspT.Fit("fit16","")


corr = g_dsspT.GetCorrelationFactor()
z = 0.5*np.log((1+corr)/(1-corr))
z_std = 1/sqrt(N_hp-3)
cor_text = " Correlation, z = %s +- %s"%(z,z_std)

g_dsspT.SetTitle('Figure S15) Non. Size, N: %s. %s'%(N_hp,cor_text))
g_dsspT.SetLineWidth(1)
g_dsspT.SetMarkerStyle(20)
g_dsspT.SetMarkerSize(0.2)
g_dsspT.GetYaxis().SetTitle("dlog(p)")
g_dsspT.GetXaxis().SetTitle("ddG")


c1.cd(2)

g_dsspT.Draw("AP")
fit16.Draw('same')

#################

Hist1 = TH1F("Hist1", "Figure S16) logP, MuMu",50,-6,6)
Hist2 = TH1F("Hist2", "Figure S17) ddG",50,-6,6)

for i in range(len(y)):

	Hist1.Fill(y[i])
	Hist2.Fill(x[i])


Hist1.SetLineWidth(2)
Hist1.SetLineColor(kBlue)
Hist1.SetAxisRange(0., 60,"Y")
c1.cd(3)
Hist1.Draw("")


Hist2.SetLineWidth(2)
Hist2.SetLineColor(kBlue)
Hist2.SetAxisRange(0., 100,"Y")
c1.cd(4)
Hist2.Draw("")

save0 = 'Images/'+'save0'+".ps"
save1 = 'Images/'+'save1'+".ps"
save2 = 'Images/'+'save2'+".ps"
save3 = 'Images/'+'save3'+".ps"

c0.Update()
c0.SaveAs(save0)

c1.Update()
c1.SaveAs(save1)

c3.Update()
c3.SaveAs(save2)

c4.Update()
c4.SaveAs(save3)




#print print_full(df)

	
raw_input( ' ... Press enter to exit ... ' )

