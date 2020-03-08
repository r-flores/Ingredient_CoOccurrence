######Data Generation#####
#
## Set variables
#

#Preferred working directory
workingDir = "~/ingredient_network/"

#Where data and code are stored
dataDir = "~/ingredient_network/data/"
codeDir = "~/ingredient_network/code/"

# Names you would like to give these files, including
# The ingredients file output from the R script (after pre-processing)
hash_infile = "ingredients.tab"

# The ingredient hash file (counting occurrence of ingredients)
hash_outfile = "ingredients.hash.txt"

#The name of the co-occurrence network 
network_rawfile = "ingredients.network.sif"

#The name of the co-occurrence network that has been filtered for only highly co-occurring ingredients
parsed_network_rawfile = "ingredients.network.parsed.sif"
#The threshold to be used to determine when an ingredient is "highly co-occurring"
edge_weight_threshold = 2000

setwd(paste0(workingDir,"data"))


##### Select one of the three datasets to use #####
#
## 1 The below commands are for the latest Open Foods Database Dataset
#
ret = download.file("https://static.openfoodfacts.org/data/en.openfoodfacts.org.products.csv","OFD_raw_data.csv")
system('wc -l OFD_raw_data.csv')
system('cut -f 1,32,35,36,37 OFD_raw_data.csv | grep \'United States\' | uniq | wc -l')
system('cut -f 1,32,35 OFD_raw_data.csv | grep \'\\d*\tUnited States\t.*\' | cut -f 1,3 > raw_ingredients.txt ')
system('cut -f 1 raw_ingredients.txt | uniq | wc -l')
raw_ingredients = read.csv("raw_ingredients.txt",sep="\t",header = TRUE)
names(raw_ingredients) <- c("barcode","ingredients_text")
ingredients <- as.data.frame(raw_ingredients[-which(raw_ingredients$ingredients_text == ""), ])
# You can now remove the raw_ingredients variable 
# if you are feeling confident in the reproducibility of your project
# To do this, uncomment the command below
#remove(raw_ingredients)

#
## 2 The below commands will download FoodData Central Dataset USDA new data source (December 2019)
#
ret = download.file("https://fdc.nal.usda.gov/fdc-datasets/FoodData_Central_branded_food_csv_2019-12-17.zip","FCD_raw_data.zip")
system('unzip FCD_raw_data.zip')
library(readr)
Products <- read_csv("branded_food.csv")
system('wc -l branded_food.csv')
ingredients <- Products[c(1,4)]
colnames(ingredients) <- c("barcode", "ingredients_text")


#####Removal of special characters (Extension of Data Pre-processing)##############
# If you are using the second or third dataset start here
ingredients$ingredients_text <- tolower(ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('\\d+[.,]*\\d*\\s{0,1}%','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('organic ','',ingredients$ingredients_text, ignore.case = TRUE)
ingredients$ingredients_text <-gsub('org ','',ingredients$ingredients_text, ignore.case = TRUE)
ingredients$ingredients_text <-gsub('certified organic','',ingredients$ingredients_text, ignore.case = TRUE)
ingredients$ingredients_text <-gsub('[éèë]','e',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('ï','i',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('[âà]','a',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('ô','o',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('_','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('?','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('\\[\\]','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('ingredient[s]*\\s{0,1}\\:*','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('ingrédient[s]*\\s{0,1}\\:*','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('ingrédiente[ns]*\\s{0,1}\\:*','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('amount','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('serving','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('nourishment','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('nourishment','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('made from ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('only ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('pure ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('contains less than','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('vital ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('cultured','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('pasteurized','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('distilled','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('california','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('grade a+','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('extra','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('virgin','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('free range','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('french','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('high fiber','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('low fat','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('whole ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('rolled ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('expeller ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('pressed ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('raw ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('evaporated ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('cow\'s milk','milk',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('non-gmo','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('refined ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('milled ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('soft ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('hard ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('toasted ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('sliced ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('mediterranean ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('peeled ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('cored ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('dry ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('roasted ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('mountain ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('spring ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('diced ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('sparkling','carbonated',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('cooked ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('freshly made from: ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('concentrated ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('stone-ground ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('partially hydrogenated','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('steel-cut ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('steel cut ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('thick cut ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('thick-cut ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('gourmet ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('fresh ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('frozen ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('natural ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('non-genetically modified ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('genetically modified','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('and nothing else','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('love ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('nothing ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('italian ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('imported ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('hand-picked ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('dried ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('crushed ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('added to preserve freshness','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('for freshness','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('to preserve freshness','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('mechanically','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('mini','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('naturally','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('h2o','water',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('\\*','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('as a preservative','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('juice from','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('exclusively','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('ground','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('shaved','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('kosher','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('ripe','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('\\.','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('nitrogen-flushed to maintain freshness','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('farm-raised','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('boneless','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('skinless','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('wild caught','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('parboiled','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('preparboiled','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('individually wrapped','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('all nat','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('wq all nat','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('young','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('aged','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('premium','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('added to retain color','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('blend of','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('in shell','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('with their juices','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('selected ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('with calcium ascorbate to promote whiteness','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('to promote freshness','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('not from concentrate','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('for color','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('real','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('propellant-free','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('gluten-free','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('grass-fed','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('grass fed','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('modified','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('pitted','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('reconstituted','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('cage\\-free','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('cage free','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('free\\-range','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('freeze-','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('wild-harvested','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('contains:','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('flavors','flavor',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('less than','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('of:','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('filtered','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('contents:','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('pre-washed','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('in-shell','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('in-the-shell','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('carragenan','carrageenan',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('carrageena','carrageenan',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('carrageenen','carrageenan',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('carrageen','carrageenan',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('enriched','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('unbleached','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('bleached','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('mono and diglycerides','monoglycerides and diglycerides',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('or less of the following','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('or less of each of the following','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('or less of','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('or less','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('distributed by','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('the','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('contains','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('with','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('textured','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('more of','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('following','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('product','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('strach','starch',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('emulsidiers','emulsifiers',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('occurring','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('or of following','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('\\(\\)','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('defatted','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('fractionated','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('hydrogenated','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('preserved','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('dry','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('modified','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('partially','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('cultured','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('or less of each of the following','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('bragg\'s ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('bart \\& judy\'s','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('nagai\'s','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('best annie\'s friends','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('annie\'s','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('what\'s inside?','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('newman\'s own','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('dr\\.\\s*kefir\'s','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('rockin\' poppin\'','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('hellmann\'s','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('us grown ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('new zealand ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('orville redenbacher\'s ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('opal ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('of modena','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('(sourced from the united kingdom)','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub(' may contain a blend of united states, brazilian, mexican, belize and south african concentrates.','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('from emilia romagna.','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('spanish','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('himalayan','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('irish\\-style','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('scott\'s','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('pedro ximenez','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('best','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('friends','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('teacher\'s ','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('barley\'s','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('baker\'s','bakers',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('brewer\'s','brewers',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('sheep\'s','sheeps',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('kroger co','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('confectioner\'s','confectioners',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('trader joe\'s','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('m\\&m\'s','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('m\\,m\'s','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('&','and',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('\\]',')',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('\\[',')',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('\\}',')',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('\\{',')',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub(';',',',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub(':',',',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('-','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('-','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('@','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('#','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('!','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('\\/',' or ',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('\\|','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('\\%','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('†','',ingredients$ingredients_text)
ingredients$ingredients_text <-gsub('`','',ingredients$ingredients_text)

######Creaate a simple Co-Occurence Network#####
write.table(ingredients, file = hash_infile,sep='\t',eol='\n',quote=FALSE,append=FALSE,row.names = FALSE)

#Create an ingredient hash and an ingredient network with the scripts provided on this repository
make_hash_command = paste0("python ",codeDir,"ingredient_hash.py ","-i ",dataDir,hash_infile," ","-o ",dataDir,hash_outfile)
system(make_hash_command)
make_network_command = paste0("python ",codeDir,"ingredient_network.py ","-i ",dataDir,hash_infile," ","-o ",dataDir,network_rawfile)
system(make_network_command)

#install.packages("igraph")
library(igraph)
g_data <- read.csv(network_rawfile,sep = " ",header = F,as.is = T,col.names = c("v1","v2","weight"))

#make sure weight is treated as numeric
g_data$weight <- as.numeric(as.character(g_data$weight))

g <- graph_from_data_frame(g_data)

##############Basic Analysis with R###############
# Check that the graph is weighted
# Should say TRUE
is_weighted(g)

# Check that the graph is not directed
# Should say FALSE
is.directed(g)

#Is the graph totally connected?
# Should say FALSE
is.connected(g)

#
# Remove loops and multiple edges
#
g <-simplify(g,remove.multiple = T,remove.loops = T)

#Get number of nodes, edges, and edge density
length(V(g))
length(E(g))
graph.density(g)

# Get the transitivity or clustering coefficient
cc <- transitivity(g,type="global")

# Is the graph totally connected? We expect not.
is.connected(g)

# Centrality Measures
# Assortativity
#ass <- assortativity_degree(g, directed=F)

# Closeness
#close <- closeness(g, mode="all", weights=NA)

# Eigenvector centrality
#eigen <- eigen_centrality(g, directed=FALSE, weights=NA)

#
#  Note: The commands below take time since the network is large, 
# so only uncomment if you can step away
#

# Get the diameter, using edge weights
#diam <- diameter(g,directed=FALSE,weights=E(g)$weight)

# Betweenness centrality
#bet <- betweenness(g, directed=FALSE, weights=NA)

deg.dist <- degree.distribution(g,cumulative = TRUE,mode="all")

#Plot the degree distribution as calculated
#This will plot the graph shown in Figure 1A
plot(deg.dist, main = "Degree Distribution, only top 1000 nodes by degree shown",
     col="darkblue",pch=10,cex=0.75,xlab="Degree", ylab="Cumulative Frequency",xlim=c(0,1000))

# Plot the log-log of the degree distribution
# This is just a re-statement of the figure in 1A but is not presented in the manuscript
plot(log(deg.dist,base=10), main = "Degree Distribution of the Ingredient Co-Occurence Network",
     col="darkblue",pch=10,cex=0.75,xlab="Log10(Degree)", ylab="Log10(Cumulative Frequency)")

#
# Compare the degree distribution of the ingredient co-occurrence network against
# a random Erdos Reyni (random) graph and a Barabasi preferential attachment simulated graph
#

# Create the graphs to compare against, an "ab" (Albert-Barabasi) network and an "er"(Erdös-Reyni) network
random.er.graph = erdos.renyi.game(n=as.numeric(length(V(g))), p.or.m = 1/1175,directed = F, loops = F)
random.ab.graph = barabasi.game(n=length(V(g)),m=37,directed=F)

# Compare vertex count
length(V(g))
length(V(random.ab.graph))
length(V(random.er.graph))

# Compare the edge count
length(E(g))
length(E(random.ab.graph))
length(E(random.er.graph))

# Compare the histograms
hist(degree(g),xlim=c(0,6000))
hist(degree(random.ab.graph))
hist(degree(random.er.graph))

# Generate degree distributions for the random graphs
random.er.degdist = degree.distribution(random.er.graph,cumulative=T,mode='all')
random.ab.degdist = degree.distribution(random.ab.graph,cumulative=T,mode='all')

# Perform the Kolmogorov Smirnov tests of the real graph against the random graphs
ks.test(deg.dist,random.er.degdist)
ks.test(deg.dist,random.ab.degdist)

#--------------------------------------------------
#Get the top X nodes by degree and their edge density
#

#Obtain the degrees for all nodes in the graph
deg <- as.data.frame(degree(g, mode="all"))
names(deg)= c("degree")
deg$vertex_id = V(g)

#Sort the dataframe by decreasing order (highest degree nodes will be first)
deg <- deg[order(deg$degree,decreasing = T),]

#We begin by examining the top 3 nodes and then increase that size by intervals, see below
top_x = 3

#Create new empty lists to hold incoming data
hub_list = list()   #A list for the hub nodes
k_list = list()     #A list for keeping track of what values of k we have tried
den_list = list()   #A list for the densities of each induced subgraph we examine
stop = FALSE        #A stop condition to tell us when we have sampled enough of the network

# This loop will run until we have sampled the top 5000 nodes in intervals
while(stop == FALSE){
  #Sanity check print statement to know your progress (its not fast, but not too slow either)
  print(c("at ",top_x))
  
  #Create a list of the top x nodes by degree
  for(i in c(1:top_x)){
    hub_list = c(hub_list,deg[i,]$vertex_id)
  }
  
  #Get the induced subgraph of the top x hub nodes in the network
  subg <- induced_subgraph(g,vids=c(hub_list))
  subg <- as.undirected(subg, mode= "collapse")
  
  #Reset the hub_list if you havent already
  hub_list=list()
  
  #Add this value of k to the list for use in plotting
  k_list = c(k_list,top_x)
  
  #Get the density of the induced subgraph used above and add it to our list of subgraph densities for k  
  den_list = c(den_list,as.numeric(graph.density(subg,loops = F)))
  
  #Start with k = 3, then go to use the top 50 nodes, increase by 50 nodes at each step until we hit 5000 nodes  
  if(top_x ==3){
    top_x = 50
  }else{
    top_x = top_x+50
  }
  if(top_x >= 5000){
    stop=TRUE
  }
}

#Round the density list for better plotting
round(den_list, digits=2)

#Plot the sampled density of the top k nodes
#This will result in the plot shown in Figure 1B
plot(k_list,den_list,type = "o",
     xlab="Top k Nodes by Degree",
     ylab="Density of the Induced Subgraph of the Top k Nodes",
     col="darkgreen",cex=0.90,pch=19,
     main = "Sampled Edge Density of the Top k Nodes"
)

#Top 400 nodes are 90%+ density; we see the difference at the split between 400 and 450 nodes.
#90% density threshold is chosen arbitrarily
den_list[9]
k_list[9]

den_list[10]
k_list[10]

top_x = 5000 #<-for edge density
hub_list = list()
for(i in c(1:top_x)){
  hub_list = c(hub_list,deg[i,]$vertex_id)
}
#Get the induced subgraph of the top i hub nodes in the network
subg <- induced_subgraph(g,vids=c(hub_list))
subg <- as.undirected(subg, mode= "collapse")

#This number should be 400
length(V(subg))

#How many edges in the induced subgraph by top 400 nodes? Whats the graph density?
graph.density(subg)

#-------------------------------------------------------------------------------
#Get the induced subgraph for the (top_x) nodes by degree and their edge density
#
deg <- as.data.frame(degree(g, mode="all"))
names(deg)= c("degree")
deg$vertex_id = V(g)
deg <- deg[order(deg$degree,decreasing = T),]
top_x = 200 #<-Change for Network Nodes
hub_list = list()
for(i in c(1:top_x)){
  hub_list = c(hub_list,deg[i,]$vertex_id)
}

#Get the induced subgraph of the top i hub nodes in the network
subg <- induced_subgraph(g,vids=c(hub_list))
subg <- as.undirected(subg, mode= "collapse")
write.graph(subg,file = paste0("top_",top_x,"_network.txt"),format="ncol")
write.csv(degree(subg),file = paste0("top_",top_x,"_node_degree.txt"))

#simple plot of the top 20 nodes
plot(subg)

# ---------------------------------------------------------------------
# Parse the network so only edges with weights higher than 2000 are used
# This step may take a few minutes
#
parse_network_command = paste0("python ",codeDir,"parse_network_by_weight.py ","-i ",dataDir,network_rawfile," ","-o ",dataDir,parsed_network_rawfile," ","-t ",edge_weight_threshold)
system(parse_network_command)

#
# Read in the parsed network
#
g_data <- read.csv(parsed_network_rawfile,sep = " ",header = F,as.is = T,col.names = c("v1","v2","weight"))
g <- graph_from_data_frame(g_data,directed = F)

# Check that the graph is weighted
# Should say TRUE
is_weighted(g)

# Check that the graph is not directed
# Should say FALSE
is.directed(g)

#Is the graph totally connected?
is.connected(g)

#Remove loops and multiple edges
g <-simplify(g,remove.multiple = T,remove.loops = T)

#Get number of nodes
length(V(g))

#Get number of edges
length(E(g))

#Edge Density
graph.density(g,loops = FALSE)

#Transitivity
cc <- transitivity(g,type="global")

#Diameter, using weights
#diam <- diameter(g,directed=FALSE)

#Diameter, not using weights
#diam <- diameter(g,directed=FALSE,weights=NA)

#assortativity
#ass <- assortativity_degree(g, directed=F)
#close <- closeness(g, mode="all", weights=NA)
#eigen <- eigen_centrality(g, directed=FALSE, weights=NA)
#bet <- betweenness(g, directed=FALSE, weights=NA)


#
# Examine and plot the degree distribution
#
deg.dist <- degree.distribution(g,cumulative = TRUE,mode="all")
plot(deg.dist, main = "Degree Distribution, Parsed Co-Occurence Network",
     col="darkblue",pch=10,cex=0.75,xlab="Degree", ylab="Cumulative Frequency",xlim=c(0,1000))



cliq<- max_cliques(g)
high_clique<- clique_num(g)

# Examining the cliq variable allows us to see that the maximal clique size found is 16, 
# which is used below to investigate their contents


cliq2<-cliques(g,min=high_clique,max=high_clique)
length(cliq2)
count<- 1
for(i in cliq2){
  cat("cliq2:", count, "\n")
  print(i)
  count<- count + 1}
