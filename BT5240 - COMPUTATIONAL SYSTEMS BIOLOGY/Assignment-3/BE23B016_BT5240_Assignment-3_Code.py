!pip install cameo --quiet
!pip install cobra --quiet
!pip install fastsl --quiet
!pip install swiglpk --quiet

# Importing cobra and cameo
import cobra
import cameo
from cobra.io import read_sbml_model # Import read_sbml_model from cobra.io
from cobra.io import load_model # Import load_model from cobra.io for using in-built models
print(cobra.__version__)
# ====================================================================================================
from cobra.io import read_sbml_model
from cobra.io import write_sbml_model
from cobra.medium import minimal_medium
from cobra.flux_analysis import pfba
from cobra.flux_analysis import flux_variability_analysis

import pandas as pd
import numpy as np
import urllib.request
import matplotlib.pyplot as plt
# ====================================================================================================
# Question 1
# Part 1
urllib.request.urlretrieve("http://bigg.ucsd.edu/static/models/iJO1366.xml", "iJO1366.xml")

model_aerobic = read_sbml_model("iJO1366.xml")
#======================================
model_anaerobic = read_sbml_model("iJO1366.xml")
o2_Exchange = model_anaerobic.reactions.get_by_id("EX_o2_e")
o2_Exchange.lower_bound = 0
# ====================================================================================================
# For the aerobic case
fba_sol_aerobic = model_aerobic.optimize()

# This model now will use the calculated minimal media as it's medium
min_medium_aerobic = minimal_medium(model_aerobic, fba_sol_aerobic.objective_value)
model_aerobic.medium = min_medium_aerobic

# Peforming FBA (Aerobic Case)
fba_sol_aerobic = model_aerobic.optimize()

print(f"Biomass reaction flux: {fba_sol_aerobic.fluxes['BIOMASS_Ec_iJO1366_core_53p95M']}")
print(fba_sol_aerobic.status)
print("="*100)

# Peforming pFBA (Aerobic Case)
pfba_aerobic_sol = pfba(model_aerobic)

print(f"Biomass reaction flux: {pfba_aerobic_sol.fluxes['BIOMASS_Ec_iJO1366_core_53p95M']}")
print(pfba_aerobic_sol.status)
print("="*100)
# ---------------------------------------------------------------------------------------------------------------------------
# Plotting the fluxes
top_fluxes_aerobic_fba = fba_sol_aerobic.fluxes[abs(fba_sol_aerobic.fluxes) > 1e-6].abs().sort_values().tail(20)

fba_sol_aerobic.fluxes[top_fluxes_aerobic_fba.index].plot.barh()
plt.xlabel("Flux Value")
plt.title("Top 20 Fluxes (Aerobic, FBA)")
plt.show()
# ---------------------------------------------------------------------------------------------------------------------------
top_fluxes_aerobic_pfba = pfba_aerobic_sol.fluxes[abs(pfba_aerobic_sol.fluxes) > 1e-6].abs().sort_values().tail(20)

pfba_aerobic_sol.fluxes[top_fluxes_aerobic_pfba.index].plot.barh()
plt.xlabel("Flux Value")
plt.title("Top 20 Fluxes (Aerobic, pFBA)")
plt.show()
# ---------------------------------------------------------------------------------------------------------------------------
# Differences in fluxes between FBA and pFBA
Aerobic_flux_diff = {}
for rxn_id in fba_sol_aerobic.fluxes.index:
    flux_dif = fba_sol_aerobic.fluxes[rxn_id] - pfba_aerobic_sol.fluxes[rxn_id]
    Aerobic_flux_diff[rxn_id] = (flux_dif)
    
Aerobic_flux_difference_df = pd.DataFrame({
    "Reaction_ID": Aerobic_flux_diff.keys(),
    "Flux differences (FBA - pFBA)": Aerobic_flux_diff.values()
})

Aerobic_flux_difference_df = Aerobic_flux_difference_df.sort_values(by="Flux differences (FBA - pFBA)",ascending=False)

Aerobic_flux_difference_df.head(10)
# ====================================================================================================
# Peforming FVA (Aerobic Case)
fva_aerobic_results = flux_variability_analysis(model = model_aerobic, reaction_list = model_aerobic.reactions, fraction_of_optimum=1.0)
fva_aerobic_results

# Creating a column range and sorting to find reactions with highly variable flux allowability
fva_aerobic_results["range"] = abs(fva_aerobic_results["maximum"] - fva_aerobic_results["minimum"])
fva_sorted = fva_aerobic_results.sort_values(by="range", ascending=False)
fva_sorted.head(50)
# ====================================================================================================
# For the anaerobic case
fba_sol_anaerobic = model_anaerobic.optimize()
print(f"Objetive value (biomass reaction flux: {fba_sol_anaerobic.objective_value})")

min_medium_anaerobic = minimal_medium(model_anaerobic, fba_sol_anaerobic.objective_value)
model_anaerobic.medium = min_medium_anaerobic

# Peforming FBA (Anaerobic Case)
fba_sol_anaerobic = model_anaerobic.optimize()

print(f"Biomass reaction flux: {fba_sol_anaerobic.fluxes['BIOMASS_Ec_iJO1366_core_53p95M']}")
print(fba_sol_anaerobic.status)
print("="*100)

# Peforming pFBA (Anaerobic Case)
pfba_anaerobic_sol = pfba(model_anaerobic)

print(f"Biomass reaction flux: {pfba_anaerobic_sol.fluxes['BIOMASS_Ec_iJO1366_core_53p95M']}")
print(pfba_anaerobic_sol.status)
print("="*100)

#---------------------------------------------------------------------------------------------------------------------------
# Plotting the fluxes
top_fluxes_anaerobic_fba = fba_sol_anaerobic.fluxes[abs(fba_sol_anaerobic.fluxes) > 1e-6].abs().sort_values().tail(20)

fba_sol_anaerobic[top_fluxes_anaerobic_fba.index].plot.barh()
plt.xlabel("Flux Value")
plt.title("Top 20 Fluxes (Anaerobic, FBA)")
plt.show()
#---------------------------------------------------------------------------------------------------------------------------
top_fluxes_anaerobic_pfba = pfba_anaerobic_sol.fluxes[abs(pfba_anaerobic_sol.fluxes) > 1e-6].abs().sort_values().tail(20)

pfba_anaerobic_sol[top_fluxes_anaerobic_pfba.index].plot.barh()
plt.xlabel("Flux Value")
plt.title("Top 20 Fluxes (Anaerobic, pFBA)")
plt.show()
#---------------------------------------------------------------------------------------------------------------------------
# Differences in fluxes between FBA and pFBA
Anaerobic_flux_diff = {}
for rxn_id in fba_sol_anaerobic.fluxes.index:
    flux_dif = fba_sol_anaerobic.fluxes[rxn_id] - pfba_anaerobic_sol.fluxes[rxn_id]
    Anaerobic_flux_diff[rxn_id] = (flux_dif)
    
Anaerobic_flux_difference_df = pd.DataFrame({
    "Reaction_ID": Anaerobic_flux_diff.keys(),
    "Flux differences (FBA - pFBA)": Anaerobic_flux_diff.values()
})

Anaerobic_flux_difference_df = Anaerobic_flux_difference_df.sort_values(by="Flux differences (FBA - pFBA)",ascending=False)

Anaerobic_flux_difference_df.head(10)
# ====================================================================================================
# Peforming FVA (Anaerobic Case)
fva_anaerobic_results = flux_variability_analysis(model = model_anaerobic, reaction_list = model_anaerobic.reactions, fraction_of_optimum=1.0)
fva_anaerobic_results

# Creating a column range and sorting to find reactions with highly variable flux allowability
fva_anaerobic_results["range"] = abs(fva_anaerobic_results["maximum"] - fva_anaerobic_results["minimum"])
fva_sorted = fva_anaerobic_results.sort_values(by="range", ascending=False)
fva_sorted.head(50)
# ====================================================================================================
# Part 2
# Intitalizations
model_aerobic = read_sbml_model("iJO1366.xml")
#======================================
model_anaerobic = read_sbml_model("iJO1366.xml")
o2_Exchange = model_anaerobic.reactions.get_by_id("EX_o2_e")
o2_Exchange.lower_bound = 0
# ====================================================================================================
Weights = []
w1 = (1-0.03, 0.03)
w2 = (1-0.05, 0.05)
w3 = (1-0.065, 0.065)

Weights.append(w1)
Weights.append(w2)
Weights.append(w3)

biomass = model_aerobic.reactions.get_by_id("BIOMASS_Ec_iJO1366_core_53p95M")
product = model_aerobic.reactions.get_by_id("EX_etoh_e")
# ====================================================================================================
# In the aerobic case
print("IN THE AEROBIC CASE")
for wb, wp in Weights:
    print(f"The weights are {wb} for biomass and {wp} for the product")
    model_aerobic.objective = {biomass: wb,product: wp}
    sol = model_aerobic.optimize()
    
    # Printing the reactions with high fluxes to gauge the differences in flux distribution:
    top_fluxes = sol.fluxes[abs(sol.fluxes) > 1e-6].abs().sort_values().tail(20)

    sol.fluxes[top_fluxes.index].plot.barh()
    plt.xlabel("Flux Value")
    plt.title(f"Top 20 Fluxes (Aerobic, [biomass: {wb},product: {wp}])")
    plt.show()

    print(f"The value of the new objective function is: {sol.objective_value}")
    print(f"Biomass reaction flux: {sol.fluxes['BIOMASS_Ec_iJO1366_core_53p95M']}")
    print(f"Ethanol reaction flux: {sol.fluxes['EX_etoh_e']}")
    print(sol.status)
    print("-"*100)

print("="*100)
print("IN THE ANAEROBIC CASE")
for wb, wp in Weights:
    print(f"The weights are {wb} for biomass and {wp} for the product")
    model_anaerobic.objective = {biomass: wb,product: wp}
    sol = model_anaerobic.optimize()
    
    # Printing the reactions with high fluxes to gauge the differences in flux distribution:
    top_fluxes = sol.fluxes[abs(sol.fluxes) > 1e-6].abs().sort_values().tail(20)

    sol.fluxes[top_fluxes.index].plot.barh()
    plt.xlabel("Flux Value")
    plt.title(f"Top 20 Fluxes (Anaerobic, [biomass: {wb},product: {wp}])")
    plt.show()
    
    print(f"The value of the new objective function is: {sol.objective_value}")
    print(f"Biomass reaction flux: {sol.fluxes['BIOMASS_Ec_iJO1366_core_53p95M']}")
    print(f"Ethanol reaction flux: {sol.fluxes['EX_etoh_e']}")
    print(sol.status)
    print("-"*100)
# ====================================================================================================   
# Question 2
# Setting up the models
urllib.request.urlretrieve("http://bigg.ucsd.edu/static/models/iNJ661.xml", "iNJ661.xml")
urllib.request.urlretrieve("http://bigg.ucsd.edu/static/models/iCN718.xml", "iCN718.xml")
    
Tuber_model = read_sbml_model("iNJ661.xml")
Actino_model = read_sbml_model("iCN718.xml")
# ====================================================================================================  
# Part 1
print(f"The number of genes in the Mycobacterium tuberculosis: {len(list(Tuber_model.genes))}")
print(f"The number of reactions in the Mycobacterium tuberculosis: {len(list(Tuber_model.reactions))}")
print(f"The number of metabolites in the Mycobacterium tuberculosis: {len(list(Tuber_model.metabolites))}")
print(f"The reaction-to-gene ratio for Mycobacterium tuberculosis is: {len(list(Tuber_model.reactions))/len(list(Tuber_model.genes))}")
print("="*100)
print(f"The number of genes in the Acinetobacter baumannii: {len(list(Actino_model.genes))}")
print(f"The number of reactions in the Acinetobacter baumannii: {len(list(Actino_model.reactions))}")
print(f"The number of metabolites in the Acinetobacter baumannii: {len(list(Actino_model.metabolites))}")
print(f"The reaction-to-gene ratio for Acinetobacter baumannii is: {len(list(Actino_model.reactions))/len(list(Actino_model.genes))}")
# ====================================================================================================  
# Part 2
# Part a, identifying essential genes
sol = Tuber_model.optimize()
Threshold = (sol.objective_value * 0.05)
MTB_essential_gene_set = cobra.flux_analysis.find_essential_genes(Tuber_model, threshold = Threshold)
print(MTB_essential_gene_set)

print("="*100)

sol = Actino_model.optimize()
Threshold = (sol.objective_value * 0.05)
ACB_essential_gene_set = cobra.flux_analysis.find_essential_genes(Actino_model, threshold = Threshold)
print(ACB_essential_gene_set)
# ==================================================================================================== 
MTB_temp = set()
ACB_temp = set()
for gene in MTB_essential_gene_set:
    MTB_temp.add(gene.id)
for gene in ACB_essential_gene_set:
    ACB_temp.add(gene.id)

print(MTB_temp)    
print("*"*125)
print(ACB_temp)
# ==================================================================================================== 
# Part B
# Determining Essential Reactions for both models
MTB_essential_reactions = set()

for gene in MTB_essential_gene_set:
    for rxn in gene.reactions:
        MTB_essential_reactions.add(rxn.id)

ACB_essential_reactions = set()

for gene in ACB_essential_gene_set:
    for rxn in gene.reactions:
        ACB_essential_reactions.add(rxn.id)

# Obtaining the common reactions
common_reactions = MTB_essential_reactions.intersection(ACB_essential_reactions)
print(common_reactions)
# ==================================================================================================== 
# Trying to determine potentially common genes
common_genes = pd.DataFrame(columns=['Reaction Name', 'Mycobacterium tuberculosis', 'Acinetobacter baumannii', 'Pathway'])
n = 0

for rxn in common_reactions:
    rxn_MTB = Tuber_model.reactions.get_by_id(rxn)
    temp_1 = []
    reaction = Tuber_model.reactions.get_by_id(rxn)
    pathway_1 = reaction.subsystem
    for gene in rxn_MTB.genes:
        if gene.id in MTB_temp:
            temp_1.append(gene.id)
    
    rxn_ACB = Actino_model.reactions.get_by_id(rxn)
    temp_2 = []
    reaction = Actino_model.reactions.get_by_id(rxn)
    pathway_2 = reaction.subsystem
    for gene in rxn_ACB.genes:
        if gene.id in ACB_temp:
            temp_2.append(gene.id)
            
    common_genes.loc[n] = [rxn, temp_1, temp_2, [pathway_1, pathway_2]]
    n+=1
            
common_genes.head(50)
# ==================================================================================================== 
# Question 3
# Part 1
# Part a, b
# Setting up the model
urllib.request.urlretrieve("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

# Initial Checking 
model_Ecoli = read_sbml_model("e_coli_core.xml")
print(f"The number of genes in the E.coli: {len(list(model_Ecoli.genes))}")
print(f"The number of reactions in the E.coli: {len(list(model_Ecoli.reactions))}")
print(f"The number of metabolites in the E.coli: {len(list(model_Ecoli.metabolites))}")
print("="*100)


# Validating the model
print(cobra.io.sbml.validate_sbml_model("e_coli_core.xml"))
print("="*100)

# Post validation
print(f"The number of genes in the E.coli: {len(list(model_Ecoli.genes))}")
print(f"The number of reactions in the E.coli: {len(list(model_Ecoli.reactions))}")
print(f"The number of metabolites in the E.coli: {len(list(model_Ecoli.metabolites))}")
print("="*100)
# ==================================================================================================== 
# Part c
import collections
import collections.abc

collections.Mapping = collections.abc.Mapping
collections.MutableMapping = collections.abc.MutableMapping
collections.Sequence = collections.abc.Sequence

from cameo import load_model
from cameo.strain_design.deterministic.linear_programming import OptKnock

model_Ecoli = read_sbml_model("e_coli_core.xml")
model_Ecoli.exchanges
# ==================================================================================================== 
# Setting up the requirements
product = model_Ecoli.reactions.get_by_id("EX_etoh_e")
biomass = model_Ecoli.reactions.get_by_id("BIOMASS_Ecoli_core_w_GAM")

# Executing OptKnock
optknock = OptKnock(model_Ecoli)

single_ko_results = optknock.run(
        target = product,
        biomass = biomass,
        max_knockouts=1,
        knock_out="genes"
    )

print("\nOptKnock result single knockout:\n")
print(single_ko_results.data_frame)
print("="*100)

optknock = OptKnock(model_Ecoli)

double_ko_results = optknock.run(
        target = product,
        biomass = biomass,
        max_knockouts=2,
        knock_out="genes"
    )

print("\nOptKnock result double knockout:\n")
print(double_ko_results.data_frame)
print("="*100)

optknock = OptKnock(model_Ecoli)

triple_ko_results = optknock.run(
        target = product,
        biomass = biomass,
        max_knockouts=3,
        knock_out="genes"
    )

print("\nOptKnock result triple knockout:\n")
print(triple_ko_results.data_frame)
print("="*100)
# ====================================================================================================
# Part 2
# patch for fastsl (!sed doesnt work locally)
file_path = r"d:\anaconda\lib\site-packages\fastsl\genes.py"

with open(file_path, "r") as f:
    file_data = f.read()

file_data = file_data.replace(
    "solver_tol = model.solver.configuration.tolerances.optimality",
    "solver_tol = 1e-9"
)

with open(file_path, "w") as f:
    f.write(file_data)

import fastsl

print("Exchange Reactions of E.coli core model \n")

for rxn in model_Ecoli.exchanges:
    met_name = list(rxn.metabolites.keys())[0].name
    lb, ub = rxn.lower_bound, rxn.upper_bound
    
    print(f"{rxn.id:<15} {met_name:<30} {lb:>8.1f} {ub:>8.1f}")
# ====================================================================================================
# Setting up the models for the question
Glu_10_model = model_Ecoli.copy()
# Setting Glucose as the sole carbon source, assumed it is consumed since it is a carbon source
Glu_10_model.reactions.get_by_id("EX_glc__D_e").lower_bound = -10
Glu_10_model.reactions.get_by_id("EX_glc__D_e").upper_bound = -10
#-------------------------------------------------------------------------------------------------------------------
Ace_10_model = model_Ecoli.copy()
# Setting Acetate as the sole carbon source, assumed it is consumed since it is a carbon source
Ace_10_model.reactions.get_by_id("EX_ac_e").lower_bound = -10
Ace_10_model.reactions.get_by_id("EX_ac_e").upper_bound = -10
Ace_10_model.reactions.get_by_id("EX_glc__D_e").lower_bound = 0
#-------------------------------------------------------------------------------------------------------------------
Min_model = model_Ecoli.copy()
# Setting the model to have minimal media
fba_sol = Min_model.optimize()

min_medium = minimal_medium(Min_model, fba_sol.objective_value)
Min_model.medium = min_medium

Glu_10_model.id = "glu_10"
Ace_10_model.id = "ace_10"
Min_model.id = "min_media"

write_sbml_model(Glu_10_model, "glu_10.xml")
write_sbml_model(Ace_10_model, "ace_10.xml")
write_sbml_model(Min_model, "min_media.xml")
# ====================================================================================================
def get_gene_reactions(model, gid):
    try:
        g = model.genes.get_by_id(gid)
        return [r.id for r in g.reactions]
    except:
        return[]
# ====================================================================================================
# Checking Fast SL normally
!fast-sl e_coli_core.xml --genes

path = r".\fast-sl-results\e_coli_core\genes\e_coli_core_double_lethal_genes.csv"

df = pd.read_csv(path, header=None)
df.head()
# ====================================================================================================
# Case 1
!fast-sl glu_10.xml --genes

path = r".\fast-sl-results\glu_10\genes\glu_10_double_lethal_genes.csv"

df_1 = pd.read_csv(path, header=None)
df_1.head()
sl_1 = set([row for row in df_1.itertuples(index=False, name=None)])
print(sl_1)
print(len(sl_1))
# ====================================================================================================
# Case 2
!fast-sl ace_10.xml --genes

path = r".\fast-sl-results\ace_10\genes\ace_10_double_lethal_genes.csv"

df_2 = pd.read_csv(path, header=None)
df_2.head()
sl_2 = set([row for row in df_2.itertuples(index=False, name=None)])
print(sl_2)
print(len(sl_2))
# ====================================================================================================
# Case 3
!fast-sl min_media.xml --genes

path = r".\fast-sl-results\min_media\genes\min_media_double_lethal_genes.csv"

df_3 = pd.read_csv(path, header=None)
df_3.head()

sl_3 = set([row for row in df_3.itertuples(index=False, name=None)])
print(sl_3)
print(len(sl_3))
# ====================================================================================================
# To find common synthetic lethals
common_syn_lethals = ((sl_1.intersection(sl_2)).intersection(sl_3))
print(common_syn_lethals)
print(len(common_syn_lethals))
# ====================================================================================================