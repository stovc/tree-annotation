######################
#        LOAD
######################


source('/home/sovchinn/tree-annotation/treeannotation.R')


setwd('/home/sovchinn/Documents/tree-analysis/ClpC1_NTD')


# load protein tree
tree <- read.iqtree(file="tree.contree")

# load the annotation dataframe
csv <- read.csv("filtered_clustered.csv")

# csv <- csv[!duplicated(csv$ID),] - remove duplicates by ID. don't use this step. might be useful for debugging later

# load the annotation data of filtered bu t not clustered hits (for the bacterial tree with phyletic pattern)
filtered_data <- read.csv("filtered_hits.csv")[, c('ID', 'assembly', 'taxid')]

# load the cluster data. which protein inherits properties from which representative
cluster_data <- read.csv("cluster_dict.csv")

# load the domain data.list of domains: which molecule they belong, domain name, start & end coordinates
domain_data <- read.csv("domains.csv")

# load genomic context data
context_data <- read_context_data()

# assign taxonomy
tree <- traverse_properties(tree, csv)

# load species tree
org_tree_full <- read.tree(file="org_trees/org_tree_full.nwk")
org_data_full <- read_delim("org_trees/org_tree_full_data.csv", delim=';')

org_tree_genus <- read.tree(file="org_trees/org_tree_genus.nwk")
org_data_genus <- read_delim("org_trees/org_tree_genus_data.csv", delim=';')

org_tree_family <- read.tree(file="org_trees/org_tree_family.nwk")
org_data_family <- read_delim("org_trees/org_tree_family_data.csv", delim=';')

org_tree_order <- read.tree(file="org_trees/org_tree_order.nwk")
org_data_order <- read_delim("org_trees/org_tree_order_data.csv", delim=';')

org_tree_class <- read.tree(file="org_trees/org_tree_class.nwk")
org_data_class <- read_delim("org_trees/org_tree_class_data.csv", delim=';')

org_tree_phylum <- read.tree(file="org_trees/org_tree_phylum.nwk")
org_data_phylum <- read_delim("org_trees/org_tree_phylum_data.csv", delim=';')

######################
#      ANNOTATE
######################
# Tree All
tree_all <- annotate_tree(
  tree, 
  start='all', 
  expand=c('Proteobacteria'), 
  threshold=1)

# Tree Actinobacteria
tree_act <- annotate_tree(
  tree, 
  start = 'Actinobacteria', 
  expand = c('Actinobacteria', 'Actinomycetia', 'Corynebacteriales'), 
  threshold = 1)

# Tree Cyanobacteria
tree_cyan <- annotate_tree(
  tree, 
  start='Cyanobacteria', 
  expand=c('Cyanobacteria', ''), 
  threshold=1)


#####################################################################

# assign paralogs

# which node of the protein tree correspond to a branch of which paralog
paralog_df = data.frame(
  node    = c(511,    1002,     1012,     689,    792), 
  paralog = c("ClpC", "ClpC2", "ClpC3", "ClpCB", "ClpB")
)

# add the paralog data on the annotated tree. 2nd argument stands for the basic paralog name
# annotated tree is needed to make a list of paralogs
tree_par <- assign_paralogs(tree_all, "ClpCB", paralog_df)
tab_tree_par = as_tibble(tree_par) # and get it in the tabular form

# form the list of paralog values of the clustered proteins and their taxid numbers to annotate the species tree
paralogs = tab_tree_par[!is.na(tab_tree_par$taxid), c('label', 'taxid', 'paralog')]
paralogs$representative = paralogs$label

# paralogs + cluster data -> all_paralogs
cluster_data2 <- full_join(cluster_data, paralogs, by='representative')
cluster_data2 <- mutate(cluster_data2, id = ifelse(is.na(id), representative, id)) # why can id be NA???
cluster_data2 <- cluster_data2[!duplicated(cluster_data2$id),]  # where do duplicates come from?
cluster_data2 <- cluster_data2[, c('id', 'paralog')]
names(cluster_data2)[1] <- 'ID'
all_paralogs <- full_join(filtered_data, cluster_data2, by='ID')
all_paralogs <- all_paralogs[, c('assembly', 'taxid', 'paralog')]

counted <- all_paralogs %>% group_by(assembly, paralog) %>% mutate(count=n())
counted <- distinct(counted)
counted <- pivot_wider(counted, names_from = 'paralog', values_from = 'count')

counted[is.na(counted)] <- 0

counted2 <- subset(counted, select = -c(assembly))
# apparently, I summurize for the case when I have several asssemblies for the same genome (should I better use median?)
counted2 <- counted2 %>% group_by(taxid) %>% summarize_all(mean)
colnames(counted2)[1] <- 'label'
counted2$label = as.character(counted2$label)

# 

org_tree_full_a = annotate_org_tree(org_tree_full, org_data_full, counted2)

org_tree_family_a = copy_annotation(org_tree_family, org_tree_full_a)
org_tree_phylum_a = copy_annotation(org_tree_phylum, org_tree_full_a)

######################
#        PLOT
######################

# BASIC
# basic with bootstraps
plot_tree(tree_all, "circular", "branch.length", tips="product", circle="bootstrap",
          filename="tree-base.pdf", width=150, height=100)
# length
plot_tree(tree_all, "circular", "branch.length", tips="product", circle="length",
          filename="tree-base-len.pdf", width=150, height=100)

plot_tree(tree_all, "equal_angle", "branch.length", context=NA,
          filename="equal_angle.pdf", width=100, height=67)

plot_tree(tree_par, "circular", "branch.length", aes_color="paralog", 
          width=150, height=100, filename="tree-paralogs.pdf")

# TAXALINK
plot_tree(tree_all, "inward_circular", "none", taxalink=T, label_nodes=T, collapse=NA, circle="bootstrap",
          width=50, height=33, filename="tree-taxalink-all.svg")
# taxalink len
plot_tree(tree_all, "inward_circular", "none", taxalink=T, label_nodes=T, collapse=NA, circle="length",
          width=50, height=33, filename="tree-taxalink-all-len.svg")

plot_tree(tree_act, "inward_circular", "none", taxalink=T, label_nodes=T, collapse="Actinobacteria",
          width=100, height=67, filename="tree-taxalink-act.svg")
plot_tree(tree_cyan, "inward_circular", "none", taxalink=T, label_nodes=T, collapse="Cyanobacteria",
          width=100, height=67, filename="tree-taxalink-cyan.svg")

plot_tree(tree, "circular", "branch.length", collapse="Actinobacteria", filename="tree-C.pdf")
plot_tree(tree, "rectangular", "branch.length", tips="product", filename="tree-tips.pdf")
plot_tree(tree, "rectangular", "none", label_ancectors="product", filename="tree-product_anc.pdf")

# domain architecture
plot_tree(tree_all, "rectangular", "none", domains=".", tips="product", filename="tree-domains.svg",
          height=200)
plot_tree(tree_act, "rectangular", "none", domains=".", tips="product", filename="tree-domains.svg")

# context
plot_tree(tree_all, "rectangular", "none", context = '.', legend='none', tips="product",
          filename="tree-cont.svg", height=200)

plot_tree(tree_act, "rectangular", "none", collapse="Actinobacteria", 
          context = '.', legend='none', 
          filename="tree-act-cont2.svg", width=250, height=900)

plot_tree(tree_cyan, "rectangular", "none", collapse="Cyanobacteria", 
          context = '.', legend='none', 
          filename="tree-cont-cyan.svg", width=250, height=160)

# orgtree
plot_org_tree(org_tree_family_a,
              collapse=NA, width=100, height=150,
              filename='org_tree_fam.svg')

plot_org_tree(org_tree_phylum_a, width=50, height=30,
              filename='org_tree_phy.svg')

plot_org_tree(org_tree_full_a, width=150, height=450,
              filename='org_tree_full.svg')

plot_org_tree(org_tree, org_tree2,
              width=250, height=500,
              filename='org_tree.svg')
