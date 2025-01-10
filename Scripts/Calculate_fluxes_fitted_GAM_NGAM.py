#!/usr/bin/env python

"""
Author Sara Benito-Vaquerizo. 
Script to calculate the metabolic fluxes and biomass yields obtained with the fitted GAM/NGAM values and 
experimentally retrieved growth rates for the Calvin cycle and reductive glycine pathway scenarios.
Python 3.9 is used as the programming language. Cobrapy and functions are used from https://github.com/opencobra/m_model_collection/
"""

# Import statements
import os
import warnings
import re
from itertools import chain
from sys import argv
import sympy
import scipy
import scipy.io
import cobra
from cobra import Model, Reaction, Metabolite
import pandas
from cobra.util.solver import linear_reaction_coefficients
import numpy as np
import pandas as pd
from pandas import DataFrame
from contextlib import suppress
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)
from cobra.sampling import OptGPSampler, ACHRSampler, sample
from cobra.medium import minimal_medium
import math
from math import sqrt
import statistics as stats
from statistics import stdev
import optlang
import cobra.util.solver as sutil


def add_new_biomass(model,x):
	"""
	Addition of new biomass reaction changing GAM


	Input:
	model: model, cobrapy model structure

	Output:
	model, cobrapy model structure with the added pathway
	"""

# Define the reactions

	r1 = Reaction("new_biomass")
	r1.name = "Biomass synthesis reaction"
	r1.subsystem = "Biomass"
	r1.lower_bound = 0.  # This is the default
	r1.upper_bound = 1000. # This is the default
	r1.add_metabolites({ model.metabolites.get_by_id("PHOSPHOLIPID_c"):- 0.0495,
	model.metabolites.get_by_id("DNA_c"):-0.031,
	model.metabolites.get_by_id("PEPTIDO_c"):-0.06,
	model.metabolites.get_by_id("CAV_c"): -0.03,
	model.metabolites.get_by_id("RNA_c"): -0.06,
	model.metabolites.get_by_id("LPS_c"): -0.034,
	model.metabolites.get_by_id("atp_c"):-x,
	model.metabolites.get_by_id("h2o_c"): -x,
	model.metabolites.get_by_id("PROTEIN_c"): -0.68,
	model.metabolites.get_by_id("CARBO_c"): -0.055,
	model.metabolites.get_by_id("pi_c"): x,
	model.metabolites.get_by_id("BIOMASS_c"):1.0,
	model.metabolites.get_by_id("h_c"): x,
	model.metabolites.get_by_id("adp_c"): x})

	model.add_reactions({r1})
	

	return model
	
def add_ftl(Model):
	"""
	Addition of new pathways to the model


	Input:
	model: model, cobrapy model structure

	Outpu:
	model, cobrapy model structure with the added pathway
	"""

# Define the reactions

	r1 = Reaction("Ftl")
	r1.name = "formate thf ligase"
	r1.subsystem = "growth methylo"
	r1.lower_bound = 0.  # This is the default
	r1.upper_bound = 1000. # This is the default
	r1.add_metabolites({ model.metabolites.get_by_id("for_c"): -1.0,
	model.metabolites.get_by_id("thf_c"):-1.0,
	model.metabolites.get_by_id("atp_c"):-1.0,
	model.metabolites.get_by_id("adp_c"): +1.0,
	model.metabolites.get_by_id("pi_c"): +1.0,
	model.metabolites.get_by_id("10fthf_c"): +1})
	
	r2 = Reaction("Fch")
	r2.name = "methenyl-thf cyclohydrolase"
	r2.subsystem = "growth methylo"
	r2.lower_bound = 0.  # This is the default
	r2.upper_bound = 1000. # This is the default
	r2.add_metabolites({ model.metabolites.get_by_id("10fthf_c"): -1.0,
	model.metabolites.get_by_id("h_c"):-1.0,
	model.metabolites.get_by_id("methf_c"): +1.0,
	model.metabolites.get_by_id("h2o_c"): +1})	


	model.add_reactions({r1,r2})
	
	

	return model



def add_all_reactions(model):
	"""
	Add the new reactions
	Input
	model, cobrapy model structure
	Output
	model, cobrapy model structure with the new reactions added
	"""
	add_ftl(model)

	return model

	
if __name__ == "__main__":
	model=cobra.io.read_sbml_model('../Data/RehMBEL1391_sbml_L3V1.xml')
	model.reactions.EX_formate_e.lower_bound=-100
	model.reactions.EX_formate_e.upper_bound=-5
	model.reactions.EX_fru_e.lower_bounds=-0.00
	model.reactions.EX_fru_e.upper_bounds=-0.00
	model.reactions.Maintenance.bounds=[3,3] #Use the fitted NGAM value
	model.reactions.EX_pbhb_e.bounds=[0,0]
	add_all_reactions(model) #Uncomment for rGlyP scenario
	model.reactions.RBPC.bounds=[0,0] #Uncomment for rGlyP scenario
	model.reactions.Biomass.knock_out()
	x=135 #Launch it at this fitted GAM
	add_new_biomass(model,x)
	#model.reactions.new_biomass.bounds=[0.0495,0.0495] #Uncomment for Calvin cycle scenario. Edit for a different growth rate
	model.reactions.new_biomass.bounds=[0.0495,0.0495] #Uncomment for rGlyP scenario. Edit for a different growth rate. Now is the growth rate of 14 h dt
	#model.reactions.GLYAMT.bounds=[0,1000] #Non-rev for Calvin cycle scenario
	
	#Model constraints to avoid futile cycles and applied known phenotypes
	model.reactions.THMDt2.bounds=[0,0]
	model.reactions.THRtr.bounds=[0,0]
	model.reactions.FOMETRi.bounds=[0,0]
	model.reactions.ASO3t2.bounds=[0,0]
	model.reactions.SUCCtr.bounds=[0,0]
	model.reactions.PTAr.bounds=[0,1000]
	model.reactions.PPAKr.bounds=[0,1000]
	model.reactions.PROt4.bounds=[0,0]
	model.reactions.SERt4.bounds=[0,0]
	model.reactions.GLUt4.bounds=[0,0]
	model.reactions.get_by_id('3HBCDH').bounds=[0,0]
	model.reactions.NAt3_1g.bounds=[0,0]
	model.reactions.CITt7.bounds=[0,0]
	model.reactions.URAt2.bounds=[0,0]
	model.reactions.BENZOTt.bounds=[0,0]
	model.reactions.GLUABUTt7.bounds=[0,0]
	model.reactions.INSt2.bounds=[0,0]
	model.reactions.ADNt2.bounds=[0,0]
	model.reactions.ADK3.bounds=[0,0]
	model.reactions.ADK4.bounds=[0,0]
	model.reactions.ASPALAt.bounds=[0,0]
	model.reactions.THRA.bounds=[0,0]
	model.reactions.MGSA.bounds = [0.0, 0.0]
	model.reactions.MDH2.bounds = [0.0, 0.0]
	model.reactions.POX.bounds = [0.0, 0.0]
	model.reactions.ICL.bounds=[0,0]
	model.reactions.PPAKr.bounds=[0,0]
	model.reactions.HACD1.bounds=[0,0]
	model.reactions.ALRTg.bounds=[0,0]
	model.reactions.URIt2.bounds=[0,0]
	model.reactions.PTA2.bounds=[0,0]
	model.reactions.CYTDtr.bounds=[0,0]
	model.reactions.ADPT.bounds=[0,0]
	model.reactions.ALRTgp.bounds=[0,0]
	model.reactions.ADPRT3.bounds=[0,0]
	model.reactions.ADPRT4.bounds=[0,0]	
	model.reactions.CYTDt2.bounds=[0,0]
	model.reactions.PCT1.bounds=[0,0]
	model.reactions.URIt2.bounds=[0,0]
	model.reactions.get_by_id('P5CD4').bounds=[0,0]
	model.reactions.P5CD5.bounds=[0,0]
	model.reactions.CYTTS3.bounds=[0,0]
	model.reactions.ISOVC.bounds=[0,0]
	model.reactions.HPYRI.bounds=[0,0]
	model.reactions.CYTTS5.bounds=[0,0]
	model.reactions.CYTTS1.bounds=[0,0]
	model.reactions.G3PD2.bounds=[0,0]
	model.reactions.NADTRHD.bounds=[0,0]
	model.reactions.ALCDgl.bounds=[0,100]
	model.reactions.ALCD19.bounds=[-100,0]
	model.reactions.EX_acal_e.bounds = [0.0, 0.0]
	model.reactions.EX_acac_e.bounds = [0.0, 0.0]
	model.reactions.EX_pyr_e.bounds = [0.0, 0.0]
	model.reactions.EX_cit_e.bounds = [0.0, 0.0]
	model.reactions.EX_icit_e.bounds = [0.0, 0.0]
	model.reactions.EX_fum_e.bounds = [0.0, 0.0]
	model.reactions.EX_mlt_e.bounds = [0.0, 0.0]
	model.objective="EX_formate_e"
	#print(model.optimize())
	#print(model.summary())
	solution=model.optimize()
	print('The biomass yield at 14 dt is:', round((-model.reactions.new_biomass.flux/model.reactions.EX_formate_e.flux)*1000, 2))
	with open('../Results/fluxes_WT_rGlyP.txt','w') as f:#Change name for CBB
		for i in model.reactions:
			f.write('{}\t{}\n'.format(i.id,solution.fluxes[i.id]))	

