#! /usr/bin/env Rscript

#Compute topological parameters of Experimentally verified, hierarchical random graphs (Al-Anzi et. al, JTB, 2017). 
#Sherif Gerges <sherif_gerges@g.harvard.edu>






#Setup, install packages.
library(igraph)

#Import the Network
input=read.csv('import_adjacency_matrix.csv',header=TRUE,row.names=1,check.names=FALSE)

#Convert adjacency matrix to an edgelist, needed for iGraph
input2= which(input==1,arr.ind=TRUE)
final=data.frame(r=rownames(input)[input2[,"row"]],
                      i=1,
                      c=colnames(input)[input2[,"col"]])
write.csv(final, file = "edgelist.csv")

#Read the edgelist.csv file, convert to data frame. 
grph =read.csv('worm_edgelist.csv',header=FALSE,check.names=FALSE)
g=graph.data.frame(grph, directed=FALSE, vertices=NULL)

#calculate parameters of network
transitivity(g)
average.path.length(g)
wtc=walktrap.community(g)
modularity(g,membership(wtc))

#Fit HRG model to network. 
hrg_model = hrg.fit (g, hrg = NULL, start = FALSE, steps = 0)

#Loop through 20000 simulation of HRG.
L_counter = 0
Clustering_counter = 0
M_counter = 0
degrees = 0

degrees_file <- NULL
clustering_file <- NULL
modularity_file <- NULL
path_length_file <- NULL

for(i in 1:20000){
  hrg_model_final = sample_hrg(hrg_model)
  
  model_wtc=walktrap.community(hrg_model_final)
  M_counter = M_counter + modularity(hrg_model_final, membership(model_wtc))
  Clustering_counter = Clustering_counter + transitivity(hrg_model_final)
  L_counter = L_counter + average.path.length(hrg_model_final)
  degrees = degree(hrg_model_final)
  
  degrees_file <- rbind(degrees_file, degrees)
  clustering_file <- rbind(clustering_file, Clustering_counter)
  modularity_file <- rbind(modularity_file, M_counter)
  path_length_file <- rbind(path_length_file, L_counter)
  
  
  degrees = 0
  L_counter = 0
  Clustering_counter = 0
  M_counter = 0
}

#Write files containing all the data. 
write.csv(degrees_file,file="_degree_distribution_file.csv")
write.csv(clustering_file,file="_clustering.csv")
write.csv(modularity_file,file="_modularity.csv")
write.csv(path_length_file,file="_path_length.csv")

#Convert the file into data frame to print out format needed. This needs to be done for each metric. 
frequency_table1 =read.csv('_degree_distribution_file.csv',header=FALSE,check.names=FALSE)
frequency_table2 =read.csv('_clustering.csv',header=FALSE,check.names=FALSE)
frequency_table3 =read.csv('_modularity.csv',header=FALSE,check.names=FALSE)
frequency_table4 =read.csv('_path_length.csv',header=FALSE,check.names=FALSE)

#Convert to table. 
final_file1 = as.data.frame(table(unlist(frequency_table1)))
final_file2 = as.data.frame(table(unlist(frequency_table2)))
final_file3 = as.data.frame(table(unlist(frequency_table3)))
final_file4 = as.data.frame(table(unlist(frequency_table4)))

#Write table, this is the final file for the cloud. 
write.csv(final_file1,file="_degree_frequency.csv")
write.csv(final_file2,file="_clustering.csv")
write.csv(final_file3,file="_modularity.csv")
write.csv(final_file4,file="_path_length.csv")

#end
