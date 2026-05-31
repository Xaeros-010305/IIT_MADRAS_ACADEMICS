# Getting all the imports
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from scipy import stats
import scipy.io as sio
import random
import pprint
import math
# =================================================================================================================================================================
# Question 1
# Part (a) Using mmread to read the mtx file
# mmread, loads a zero based indexing eventhough the data is 1 based indexing.
adj_mat = sio.mmread("bio-yeast.mtx")
G = nx.from_scipy_sparse_array(adj_mat)
# Since the graph is in 0 based indexing, I am relabeling it back to 1 based indexing to match the dataset
G = nx.relabel_nodes(G, lambda x: x + 1)

print("Number of nodes:", G.number_of_nodes())
print("Number of edges:", G.number_of_edges())
print("The density of the above graph is:",nx.density(G))
print("The diameter of the above graph is:",nx.diameter(G))
print("="*125)
# =================================================================================================================================================================
# Part (b)
dc = (nx.degree_centrality(G))
bc = (nx.betweenness_centrality(G))
cc = (nx.closeness_centrality(G))

dc_dist = [v for v in dc.values()]
bc_dist = [v for v in bc.values()]
cc_dist = [v for v in cc.values()]

# Plotting the distribution of the different centrality measures
plt.hist(dc_dist, bins = 15)
plt.xlabel("Degree centrality")
plt.ylabel("Number of nodes")
plt.title("Centrality Measure (Degree centrality)")
plt.show()

plt.hist(bc_dist, bins = 15)
plt.xlabel("Betweeness centrality")
plt.ylabel("Number of nodes")
plt.title("Centrality Measure (Betweeness centrality)")
plt.show()

plt.hist(cc_dist, bins = 15)
plt.xlabel("Closeness centrality")
plt.ylabel("Number of nodes")
plt.title("Centrality Measure (Closeness centrality)")
plt.show()
# =================================================================================================================================================================
# Part (c)
# Floored to consider <= 5%
Threshold = int((5*G.number_of_nodes())/100)

print("Top 5% of nodes in the network with the highest degree centrality shown as a list of (degree_centrality, node_ID)")
dc_list = [(dc[key],key) for key in dc.keys()]
dc_list.sort(reverse=True)
dc_putative={k:v for v,k in dc_list[:Threshold]}
print(dc_list[:Threshold])
print("="*125)

print("Top 5% of nodes in the network with the highest betweenness centrality shown as a list of (betweenness_centrality, node_ID)")
bc_list = [(bc[key],key) for key in bc.keys()]
bc_list.sort(reverse=True)
bc_putative={k:v for v,k in bc_list[:Threshold]}
print(bc_list[:Threshold])
print("="*125)

print("Top 5% of nodes in the network with the highest closeness centrality shown as a list of (closeness_centrality, node_ID)")
cc_list = [(cc[key],key) for key in cc.keys()]
cc_list.sort(reverse=True)
cc_putative={k:v for v,k in cc_list[:Threshold]}
print(cc_list[:Threshold])
print("="*125)
# =================================================================================================================================================================
# Part (d)

True_largest = len(max(nx.connected_components(G), key=len))
# Removing nodes based on degree centrality (node removal without recomputation of centrality measure)
adj_mat = sio.mmread("bio-yeast.mtx")
G1 = nx.from_scipy_sparse_array(adj_mat)
G1 = nx.relabel_nodes(G1, lambda x: x + 1)

for n in dc_putative.keys():
    G1.remove_node(n)

largest_cc_dcr = max(nx.connected_components(G1), key=len)
print(f"The size of the largest connected component before removal of any nodes is: {True_largest}")
print(f"The size of the largest connected component after removal of putative proteins (degree centrality based) is: {len(largest_cc_dcr)}")

# Removing nodes based on between centrality (node removal without recomputation of centrality measure)
adj_mat = sio.mmread("bio-yeast.mtx")
G1 = nx.from_scipy_sparse_array(adj_mat)
G1 = nx.relabel_nodes(G1, lambda x: x + 1)

for n in bc_putative.keys():
    G1.remove_node(n)
    
largest_cc_bcr = max(nx.connected_components(G1), key=len)
print(f"The size of the largest connected component before removal of any nodes is: {True_largest}")
print(f"The size of the largest connected component after removal of putative proteins (betweeness centrality based) is: {len(largest_cc_bcr)}")

# Removing nodes based on between centrality (node removal without recomputation of centrality measure)
adj_mat = sio.mmread("bio-yeast.mtx")
G1 = nx.from_scipy_sparse_array(adj_mat)
G1 = nx.relabel_nodes(G1, lambda x: x + 1)

for n in cc_putative.keys():
    G1.remove_node(n)
    
largest_cc_ccr = max(nx.connected_components(G1), key=len)
print(f"The size of the largest connected component before removal of any nodes is: {True_largest}")
print(f"The size of the largest connected component after removal of putative proteins (closeness centrality based) is: {len(largest_cc_ccr)}")
# =================================================================================================================================================================
# Part (e)
Threshold = int((5*G.number_of_nodes())/100)
largest_cc_distribution = []

for n in range (100):
    
    adj_mat = sio.mmread("bio-yeast.mtx")
    G1 = nx.from_scipy_sparse_array(adj_mat)
    G1 = nx.relabel_nodes(G1, lambda x: x + 1)  # Nodes follow 1 based indexing
    occured = {}
    
    while len(occured) < Threshold:
        x = random.randint(1, 1458)
        if occured.get(x) is None:
            G1.remove_node(x)
            occured[x] = True
            
    largest_cc = max(nx.connected_components(G1), key=len)
    largest_cc_distribution.append(len(largest_cc))


plt.hist(largest_cc_distribution, bins = 15)
plt.xlabel("Largest Connected Component Size")
plt.ylabel("Frequency")
plt.title("Random Deletion of Nodes, Lagest Connected Component Size Frequency Distribution")
plt.show()

# We identify the distribution above as the control to obtain the Z-score of the 3 deletion cases
# Testing Degree Centrality, alpha = 0.05, two-tailed test
cntrl = pd.Series(largest_cc_distribution)

cntrl_mean, cntrl_std = cntrl.mean(), cntrl.std()
z_score = ((len(largest_cc_dcr))-cntrl_mean)/cntrl_std

alpha = 0.05
z_crit = stats.norm.ppf(1 - (alpha/2))
p_value = stats.norm.sf(abs(z_score))*2

print(f'Calculated Z score = {round(z_score,3)}\nCritical Z value  = {round(z_crit,3)}\n')

if z_crit> abs(z_score):
    print(f"Fail to reject the null hypothesis (p = {p_value})")
else:
    print(f"Reject the null hypothesis (p = {p_value})")

# Testing Betweeness Centrality, alpha = 0.05, two-tailed test
cntrl = pd.Series(largest_cc_distribution)

cntrl_mean, cntrl_std = cntrl.mean(), cntrl.std()
z_score = ((len(largest_cc_bcr))-cntrl_mean)/cntrl_std

alpha = 0.05
z_crit = stats.norm.ppf(1 - (alpha/2))
p_value = stats.norm.sf(abs(z_score))*2

print(f'Calculated Z score = {round(z_score,3)}\nCritical Z value  = {round(z_crit,3)}\n')

if z_crit> abs(z_score):
    print(f"Fail to reject the null hypothesis (p = {p_value})")
else:
    print(f"Reject the null hypothesis (p = {p_value})")

# Testing Closeness Centrality, alpha = 0.05, two-tailed test
cntrl = pd.Series(largest_cc_distribution)

cntrl_mean, cntrl_std = cntrl.mean(), cntrl.std()
z_score = ((len(largest_cc_ccr))-cntrl_mean)/cntrl_std

alpha = 0.05
z_crit = stats.norm.ppf(1 - (alpha/2))
p_value = stats.norm.sf(abs(z_score))*2

print(f'Calculated Z score = {round(z_score,3)}\nCritical Z value  = {round(z_crit,3)}\n')

if z_crit> abs(z_score):
    print(f"Fail to reject the null hypothesis (p = {p_value})")
else:
    print(f"Reject the null hypothesis (p = {p_value})")
# =================================================================================================================================================================
# Part (f), determining nodes with low degree, and high betweeness centrality
adj_mat = sio.mmread("bio-yeast.mtx")
G1 = nx.from_scipy_sparse_array(adj_mat)
G1 = nx.relabel_nodes(G1, lambda x: x + 1)

degrees = [d for n, d in G1.degree()]

plt.hist(degrees, bins = 20)
plt.xlabel("Degree (k)")
plt.ylabel("Number of nodes N(k)")
plt.title("Degree Distribution")
plt.show()

# From the above we decide that any node with 5 or lower edges is considered to have a relatively low degree
degree_dist = [n for n, d in G1.degree() if d<=5]
bc = (nx.betweenness_centrality(G1))
bw_low_deg = []
for node in degree_dist:
    bw_low_deg.append((bc[node],node))

bw_low_deg.sort(reverse = True)
print(bw_low_deg[:5])

# A more constrained test just to check
# From the above we decide that any node with 3 or lower edges is considered to have a relatively low degree
degree_dist = [n for n, d in G1.degree() if d<=3]
bc = (nx.betweenness_centrality(G1))
bw_low_deg = []
for node in degree_dist:
    bw_low_deg.append((bc[node],node))

bw_low_deg.sort(reverse = True)
print(bw_low_deg[:5])
# =================================================================================================================================================================
# Part (g)
# Reading the edgelist
G_yeast = nx.read_edgelist("yeast.edgelist")
print("Number of nodes:", G_yeast.number_of_nodes())
print("Number of edges:", G_yeast.number_of_edges())

# Finding nodes with the highest closeness centrality
cc_yeast = (nx.closeness_centrality(G_yeast))
cc_y_list = [(cc_yeast[key],key) for key in cc_yeast.keys()]
cc_y_list.sort(reverse=True)
print(f"The protein {cc_y_list[0][1]} has the highest closeness centrality = {cc_y_list[0][0]}")
print(f"The protein {cc_y_list[1][1]} has the second highest closeness centrality = {cc_y_list[1][0]}")
print(f"The protein {cc_y_list[2][1]} has the third highest closeness centrality = {cc_y_list[2][0]}")
print(f"The protein {cc_y_list[3][1]} has the forth highest closeness centrality = {cc_y_list[3][0]}")
print(f"The protein {cc_y_list[4][1]} has the fifth highest closeness centrality = {cc_y_list[4][0]}")
# ******************************************************************************************************************************************************************
# Question 2
# Part (a)
# Loading Zachary's Karate Club in graph form
G = nx.karate_club_graph()

# Setting only upper bound on communities (not enforcing any lower bounds)
c1 = nx.community.greedy_modularity_communities(G, best_n= 2)
modularity_1 = nx.community.modularity(G, c1)

print("For part (a), (i), maximum number of communities = 2")
for k in range (len(c1)):
    print(f"Community {k+1}: {set(c1[k])}, it's size is: {len(c1[k])}")
print(f"The modularity value is: {modularity_1}")
print("="*125)

c2 = nx.community.greedy_modularity_communities(G, best_n= 3)
modularity_2 = nx.community.modularity(G, c2)

print("For part (a), (ii), maximum number of communities = 3")
for k in range (len(c2)):
    print(f"Community {k+1}: {set(c2[k])}, it's size is: {len(c2[k])}")
print(f"The modularity value is: {modularity_2}")
print("="*125)

c3 = nx.community.greedy_modularity_communities(G, best_n= 4)
modularity_3 = nx.community.modularity(G, c3)

print("For part (a), (iii), maximum number of communities = 4")
for k in range (len(c3)):
    print(f"Community {k+1}: {set(c3[k])}, it's size is: {len(c3[k])}")
print(f"The modularity value is: {modularity_3}")
print("="*125)

# To analyse the algorithm for Question 2, part (e) we also set a lower cutoff bound for number of communities
G = nx.karate_club_graph()

# Setting only upper bound on communities (not enforcing any lower bounds)
c1 = nx.community.greedy_modularity_communities(G, best_n= 2, cutoff=2)
modularity_1 = nx.community.modularity(G, c1)

print("For part (a), (i), number of communities = 2")
for k in range (len(c1)):
    print(f"Community {k+1}: {set(c1[k])}, it's size is: {len(c1[k])}")
print(f"The modularity value is: {modularity_1}")
print("="*125)

c2 = nx.community.greedy_modularity_communities(G, best_n= 3, cutoff=3)
modularity_2 = nx.community.modularity(G, c2)

print("For part (a), (ii), number of communities = 3")
for k in range (len(c2)):
    print(f"Community {k+1}: {set(c2[k])}, it's size is: {len(c2[k])}")
print(f"The modularity value is: {modularity_2}")
print("="*125)

c3 = nx.community.greedy_modularity_communities(G, best_n= 4, cutoff=4)
modularity_3 = nx.community.modularity(G, c3)

print("For part (a), (iii), number of communities = 4")
for k in range (len(c3)):
    print(f"Community {k+1}: {set(c3[k])}, it's size is: {len(c3[k])}")
print(f"The modularity value is: {modularity_3}")
print("="*125)
# =================================================================================================================================================================
# Part (b)
# Loading Zachary's Karate Club in graph form
G = nx.karate_club_graph()
iter_coms = {}

gn_coms = nx.community.girvan_newman(G)
iter = 0
for c in gn_coms:
    iter += 1
    sizes = [len(com) for com in c]
    sizes.sort()
    modularity = nx.community.modularity(G, c)
    
    print(f"Partition number {iter}")
    print(f"Number of communities: {len(c)}")
    for k in range (len(c)):
        print(f"Community {k+1}: {set(c[k])}")
    print(f"The sizes of the community are: {sizes}")
    print(f"The modularity is: {modularity}")
    print("="*125)
    
    iter_coms[iter] = c

    if iter >= 5:
        break
# =================================================================================================================================================================
# Part (c), looking at the acutal classifications, visualizing the graph
# Loading Zachary's Karate Club in graph form
G = nx.karate_club_graph()

for n in G.nodes(): 
    print(n,"=",G.nodes[n]['club'])

nx.draw_circular(G, with_labels=True)
plt.show()

G = nx.karate_club_graph()

Hi = set()
Off = set()

for n in G.nodes(): 
    if G.nodes[n]['club'] == "Mr. Hi":
        Hi.add(n)
    else:
        Off.add(n)
   
print("The true communities as per labels are:")
print(f"Mr.Hi community: {Hi}")
print(f"Officer community: {Off}")

Greed_mod_hi = {0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 16, 17, 19, 21}
Greed_mod_off = {8, 14, 15, 18, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33}

GN_hi = {0, 1, 3, 4, 5, 6, 7, 10, 11, 12, 13, 16, 17, 19, 21}
GN_off  ={2, 8, 9, 14, 15, 18, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33}

accuracy_greedy = 100 - ((len(Greed_mod_hi-Hi) + len(Greed_mod_off-Off))*(100/34))
accuracy_GN = 100 - ((len(GN_hi-Hi) + len(GN_off-Off))*(100/34))

print(f"The accuracy of the Greedy modularity algorithm in detecting the communities is: {accuracy_greedy}")
print(f"The accuracy of the Girvan-Newman algorithm in detecting the communities is: {accuracy_GN}")

# True and maximum modularity of the algorithms
true_coms = [Hi,Off]
true_modularity = nx.community.modularity(G, true_coms)
print(f"The true modularity is: {true_modularity}")
print("="*125)

# Maximum modularity was obtained in the 4th iteration of Girvan-Newman
max_mod_coms = iter_coms[4]
max_mod = nx.community.modularity(G, max_mod_coms)
print(f"The max modularity in Girvan-Newman is: {max_mod}")
print(f"The max modularity using the Greedy modularity method is: {modularity_3}")
print("="*125)
# ******************************************************************************************************************************************************************
# Question 3
# Part (a)
H_Data = pd.read_excel('q3-edges.xlsx', sheet_name="Hospital")
M_Data = pd.read_excel('q3-edges.xlsx', sheet_name="Metro")

G_H = nx.from_pandas_edgelist(H_Data, 'v1', 'v2', edge_attr='Weight')
G_M = nx.from_pandas_edgelist(M_Data, 'v1', 'v2', edge_attr='Weight')


degree_H = [d for n, d in G_H.degree()]
degree_M = [d for n, d in G_M.degree()]

# Plotting the histogram
plt.hist(degree_H, bins=10)
plt.xlabel("Degree (k)")
plt.ylabel("Number of nodes N(k)")
plt.title("Degree Distribution (Hospital)")
plt.show()

# Plotting the histogram
plt.hist(degree_M, bins=10)
plt.xlabel("Degree (k)")
plt.ylabel("Number of nodes N(k)")
plt.title("Degree Distribution (Metro) ")
plt.show()
# =================================================================================================================================================================
# Part (b)
# To plot average clustering coefficient vs degree

data_H, data_M = {},{}

for n,d in G_H.degree():
    if data_H.get(d) is None:
        data_H[d] = (nx.clustering(G_H, n),)
    else:
        data_H[d] = data_H[d]+(nx.clustering(G_H, n),)

for n,d in G_M.degree():
    if data_M.get(d) is None:
        data_M[d] = (nx.clustering(G_M, n),)
    else:
        data_M[d] = data_M[d]+(nx.clustering(G_M, n),)

H_x, H_y = [],[]

for deg in data_H.keys():
    H_x.append(deg)
    H_y.append(sum(data_H[deg])/len(data_H[deg]))

plt.bar(H_x, H_y)

plt.title("Average Clustering Coefficient vs Degree (Hospital)")
plt.xlabel("Degree")
plt.xticks(range(1,30,2))
plt.xlim(0,30)
plt.ylabel("Clustering Coefficient")

plt.show()
print("="*125)

M_x, M_y = [],[]

for deg in data_M.keys():
    M_x.append(deg)
    M_y.append(sum(data_M[deg])/len(data_M[deg]))

plt.bar(M_x, M_y)

plt.title("Average Clustering Coefficient vs Degree (Metro)")
plt.xlabel("Degree")
plt.xticks(range(1,30,2))
plt.xlim(0,30)
plt.ylabel("Clustering Coefficient")

plt.show()
print("="*125)
# =================================================================================================================================================================
# Part (c)

avg_deg_H = (2*G_H.number_of_edges())/G_H.number_of_nodes()
avd_clu_H = nx.average_clustering(G_H)
cpl_H = nx.average_shortest_path_length(G_H)

avg_deg_M = (2*G_M.number_of_edges())/G_M.number_of_nodes()
avd_clu_M = nx.average_clustering(G_M)
cpl_M = nx.average_shortest_path_length(G_M)

print(f"The average node degree in the microbial association networks from the Hospital is: {avg_deg_H}")
print(f"The average clustering coefficient in the microbial association networks from the Hospital is: {avd_clu_H}")
print(f"The average characteristic path length  in the microbial association networks from the Hospital is: {cpl_H}")
print("="*125)
print(f"The average node degree in the microbial association networks from the Metro is: {avg_deg_M}")
print(f"The average clustering coefficient in the microbial association networks from the Metro is: {avd_clu_M}")
print(f"The average characteristic path length  in the microbial association networks from the Metro is: {cpl_M}")
# =================================================================================================================================================================
# Part (d)
# The intersection file is obtained using cytoscape and is loaded here to visualize the how the normalized correlation values differ between environments.
Intersection = pd.read_csv('Intersection_Q3.csv')
nodes = np.array(Intersection.iloc[:, 1])
xH = np.array(Intersection.iloc[:, 5])
xM = np.array(Intersection.iloc[:, 6])

xDiff = xH-xM

# Horizontal bar plots for xH 
plt.figure(figsize=(10,8))
bars = plt.barh(nodes, xH, color = "yellow", edgecolor = "black")
plt.xlim(-0.35, 0.35);

# Horizontal bar plots vs xM
plt.figure(figsize=(10,8))
bars = plt.barh(nodes, xM, color = "red", edgecolor = "black")
plt.xlim(-0.35, 0.35);

# Horizontal bar plots vs xDiff
plt.figure(figsize=(10,8))
bars = plt.barh(nodes, xDiff, color = "orange", edgecolor = "black")
plt.xlim(-0.35, 0.35);
# ******************************************************************************************************************************************************************