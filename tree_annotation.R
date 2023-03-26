library(tidyr)
library(dplyr)
library(readr)
library(data.table)
library("ggplot2")
library("treeio")
library("ggtree")
library("gggenes")
library("ape")
library(ggtreeExtra)


# IMPORT

read_context_data <- function(path) {
  context_data <- read_delim(path, delim=';')
  #context_data <- read.csv("cont.csv")
  
  context_data$start = as.numeric(context_data$start)
  context_data$end = as.numeric(context_data$end)
  context_data$strand = as.numeric(context_data$strand)
  grouped = context_data %>% group_by(molecule)
  grouped = grouped %>% mutate(med = median(start))
  print(nrow(grouped))
  grouped <- grouped %>% filter(abs(start - med) < 4000)
  print(nrow(grouped))
  grouped <- grouped %>% filter(abs(end - med) < 4000)
  print(nrow(grouped))
  context_data = ungroup(grouped)
  print(nrow('ungrouped'))
  context_data$start = context_data$start - context_data$med
  print(nrow('st - med'))
  context_data$end = context_data$end - context_data$med
  print(nrow('end - med'))
  names(context_data)[2] = 'motif'
  context_data
}


# TREE ANNOTATION

annotate_tree <- function(tree, csv) {
  tree <- traverse_property(tree, csv, "phylum")
  tree <- traverse_property(tree, csv, "class")
  tree <- traverse_property(tree, csv, "order")
  tree <- traverse_property(tree, csv, "family")
  tree <- traverse_property(tree, csv, "genus")
  tree <- traverse_property(tree, csv, "species")
  
  tree <- traverse_property(tree, csv, "product")
  tree <- traverse_property(tree, csv, "gene")
  
  tree <- traverse_property(tree, csv, "assembly", recursive = F)
  tree <- traverse_property(tree, csv, "length", recursive = F)
  tree <- traverse_property(tree, csv, "taxid", recursive = F)
  
  tree
}


traverse_property <- function(tree, df, taxon, recursive = T) {
  
  initial_tree <- tree
  
  tree <- as_tibble(tree)  # convert the tree to a tabular form
  
  # extract taxonomy info from dataframe
  dat <- df[, c("ID", taxon)]
  
  names(dat)[2] <- "tax"
  
  # assign taxons to tips
  tree <- tree %>% 
    merge(dat, by.x = "label", by.y = "ID", all.x = T)
  
  tree["tax_anc"] <- NA
  
  # select tips
  assignment_tree <- tree
  if (recursive) {
    assignment_tree <- traverse(tree)
    assignment_tree <- assignment_tree %>%
      select(node, tax, tax_anc)
    names(assignment_tree)[2] <- taxon
    names(assignment_tree)[3] <- paste(taxon, "anc", sep="_")
  } else {
    assignment_tree <- assignment_tree %>%
      select(node, tax)
    names(assignment_tree)[2] <- taxon
  }
  
  tree <- full_join(initial_tree, assignment_tree, by = "node")
  tree
}


traverse <- function(tree) {
  # traverse the tree from leafs to root to assign parents the taxonomic values of children if they coinside
  
  # select all parents with their children as a list
  parents <- split(tree$node, tree$parent)
  
  # convert the parent list to data frame (parent - child1 - child2)
  parent <- names(parents)
  child1 <- sapply(parent, function(x, parents) {parents[[x]][[1]]}, parents=parents)
  child2 <- sapply(parent, function(x, parents) {parents[[x]][[2]]}, parents=parents)
  parent <- as.numeric(parent)
  parents_df <- data.frame(parent, child1, child2)
  
  # assign children taxa to the parent dataframe
  tax_dat <- tree[, c("node", "tax")]
  names(tax_dat)[2] <- "tax1"
  parents_df <- merge(parents_df, tax_dat, by.x = "child1", by.y = "node", all.x = T)
  names(tax_dat)[2] <- "tax2"
  parents_df <- merge(parents_df, tax_dat, by.x = "child2", by.y = "node", all.x = T)
  
  # assign taxonomy to parents
  parents_df <- parents_df %>% 
    mutate(tax = ifelse(is.na(tax1) | is.na(tax2), NA,
                        ifelse(tax1 == tax2, tax1, "None"))) %>% 
    mutate(child_anc1 = ifelse(tax1 != tax2 & tax1 != "None" & !is.na(tax2), child1, NA)) %>%
    mutate(tax_anc1 = ifelse(tax1 != tax2 & tax1 != "None" & !is.na(tax2), tax1, NA)) %>%
    mutate(child_anc2 = ifelse(tax2 != tax1 & tax2 != "None" & !is.na(tax1), child2, NA)) %>%
    mutate(tax_anc2 = ifelse(tax2 != tax1 & tax2 != "None" & !is.na(tax1), tax2, NA))
  
  # MERGE TAX
  assignment_df <- parents_df %>% 
    select(parent, tax)
  
  names(assignment_df)[1] <- "node"
  
  tree <- tree %>% 
    merge(assignment_df, by = "node", all.x = T) %>% 
    mutate(tax = ifelse(is.na(tax.y), tax.x, tax.y)) %>% 
    select(!c(tax.x, tax.y))
  
  # MERGE TAX_ANC1
  assignment_df <- parents_df %>% 
    select(child_anc1, tax_anc1)
  
  names(assignment_df)[1] <- "node"
  names(assignment_df)[2] <- "tax_anc"
  
  tree <- tree %>% 
    merge(assignment_df, by = "node", all.x = T) %>% 
    mutate(tax_anc = ifelse(is.na(tax_anc.y), tax_anc.x, tax_anc.y)) %>% 
    select(!c(tax_anc.x, tax_anc.y))
  
  # MERGE TAX_ANC2
  assignment_df <- parents_df %>% 
    select(child_anc2, tax_anc2)
  
  names(assignment_df)[1] <- "node"
  names(assignment_df)[2] <- "tax_anc"
  
  tree <- tree %>% 
    merge(assignment_df, by = "node", all.x = T) %>% 
    mutate(tax_anc = ifelse(is.na(tax_anc.y), tax_anc.x, tax_anc.y)) %>% 
    select(!c(tax_anc.x, tax_anc.y))
  
  print(sum(is.na(tree["tax"])))
  
  if (sum(is.na(tree["tax"])) > 1) {
    tree <- traverse(tree)
  } else {
    tree
  }
  
  tree
}


# Annotate taxonomy

calibrate_taxonomy <- function(tree, start, expand, threshold=1) {
  
  tree <- assign_taxonomy(tree, start, expand)
  
  tree <- filter_taxonomy(tree, threshold)
  tree <- filter_genomes(tree)
  
  tree
}


assign_taxonomy <- function(tree, start, tax_vect) {
  
  tab_tree <- as_tibble(tree)
  
  if (start == "all") {
    tab_tree <- tab_tree %>% 
      mutate(taxon = ifelse(!(phylum %in% tax_vect), phylum, 
                            ifelse(!(class %in% tax_vect), class,
                                   ifelse(!(order %in% tax_vect), order, 
                                          ifelse(!(family %in% tax_vect), family, genus)))))
  } else {
    tab_tree <- tab_tree %>% 
      mutate(taxon = ifelse(phylum != start, "None", 
                            ifelse(!(phylum %in% tax_vect), phylum,
                                   ifelse(!(class %in% tax_vect), class, 
                                          ifelse(!(order %in% tax_vect), order, 
                                                 ifelse(!(family %in% tax_vect), family, 
                                                        ifelse(!(genus %in% tax_vect), genus, species)))))))
  }
  
  
  
  
  
  assignment_tree <- tab_tree
  
  assignment_tree <- assignment_tree %>%
    select(node, taxon)
  
  tree <- full_join(tree, assignment_tree, by = "node")
  tree
}


filter_genomes <- function(tree) {
  
  tab_tree <- as_tibble(tree)
  
  tab_tree <- tab_tree %>%
    mutate(filt_genome = ifelse(taxon == 'None', NA, assembly))
  
  assignment_tree <- tab_tree
  
  assignment_tree <- assignment_tree %>%
    select(node, filt_genome)
  
  tree <- full_join(tree, assignment_tree, by = "node")
  tree
}


filter_taxonomy <- function(tree, threshold) {
  
  tab_tree <- as_tibble(tree)
  
  tab_tree <- tab_tree %>% 
    group_by(taxon) %>%
    mutate(taxonomy = ifelse(n() > threshold, taxon, "Other")) %>% 
    ungroup()
  
  assignment_tree <- tab_tree
  
  assignment_tree <- assignment_tree %>%
    select(node, taxonomy)
  
  tree <- full_join(tree, assignment_tree, by = "node")
  tree
}




traverse_mean <- function(tree) {
  
  assignment_df <- tree[, c(-2:-6)] %>% 
    group_by(parent) %>% 
    summarize(across(everything(), ~ mean(.x))) %>% 
    ungroup()
  
  names(assignment_df)[1] = 'node'
  
  tree_outer <- tree %>% 
    anti_join(assignment_df, by = "node")
  
  tree_inner <- tree[, 1:6] %>% 
    right_join(assignment_df, by = "node")
  
  tree2 <- bind_rows(tree_outer, tree_inner)
  
  nas = sum(is.na(tree2[, -1:-6]))
  
  print(nas)
  
  if (nas > 13) {
    tree2 <- traverse_mean(tree2)
  } else {
    tree2
  }
  
  tree2
}








#' assign paralog annotation to all descendants of nodes from dataframe
#' 
#' nonspecified nodes are assigned from from the 'protein' argument
assign_paralogs <- function(tree, protein, paralod_df) {
  # get tabular tree for further manipulations
  tab_tree <- as_tibble(tree)
  
  # assign the main paralog (i.e. paralog = "ClpP")
  nodes <- tab_tree$node
  d <- data.frame(node = nodes, paralog = protein)
  tab_tree <- full_join(tab_tree, d, by = "node")
  
  # itereate over node:paralog_name pairs and assign them to the subtries
  for (i in 1:nrow(paralog_df)) {
    node = paralod_df[i, 'node']
    paralog = paralod_df[i, 'paralog']
    
    offs <- offspring(tab_tree, node)
    offs$paralog <- paralog
    offs <- offs[, c('node', 'paralog')]
    tab_tree <- full_join(tab_tree, offs, by='node')
    tab_tree <- tab_tree %>% 
      mutate(paralog = ifelse(is.na(paralog.y), paralog.x, paralog.y)) %>% 
      select(-paralog.x, -paralog.y)
  }
  tab_tree <- tab_tree[, c('node', 'paralog')]
  tree <- full_join(tree, tab_tree, by='node')
  
  tree
}


annotate_org_tree <- function(org_tree, org_data, prot_tree) {
  
  # PART I - MAKE PARALOG ASSIGNMENT DF
  
  # add the paralog data on the annotated tree. 2nd argument stands for the basic paralog name
  # annotated tree is needed to make a list of paralogs
  
  tab_tree_par = as_tibble(prot_tree) # and get it in the tabular form
  
  # form the list of paralog values of the clustered proteins and their taxid numbers to annotate the species tree
  paralogs = tab_tree_par[!is.na(tab_tree_par$taxid), c('label', 'taxid', 'paralog')]
  paralogs$representative = paralogs$label
  
  # paralogs + cluster data -> paralogs_all_nodes
  paralogs_representative <- full_join(cluster_data, paralogs, by='representative')
  paralogs_representative <- mutate(paralogs_representative, id = ifelse(is.na(id), representative, id)) # why can id be NA???
  paralogs_representative <- paralogs_representative[!duplicated(paralogs_representative$id),]  # where do duplicates come from?
  paralogs_representative <- paralogs_representative[, c('id', 'paralog')]
  names(paralogs_representative)[1] <- 'ID'
  
  paralogs_all_orgs <- full_join(filtered_data, paralogs_representative, by='ID')
  paralogs_all_orgs <- paralogs_all_orgs[, c('assembly', 'taxid', 'paralog')]
  
  paralogs_counted <- paralogs_all_orgs %>% group_by(assembly, paralog) %>% mutate(count=n())
  paralogs_counted <- distinct(paralogs_counted)
  paralogs_counted <- pivot_wider(paralogs_counted, names_from = 'paralog', values_from = 'count')
  
  paralogs_counted[is.na(paralogs_counted)] <- 0
  paralogs_counted <- subset(paralogs_counted, select = -c(assembly))
  
  # apparently, I summurize for the case when I have several asssemblies for the same genome (should I better use median?)
  paralogs_summed <- paralogs_counted %>% group_by(taxid) %>% summarize_all(mean)
  colnames(paralogs_summed)[1] <- 'label'
  paralogs_summed$label = as.character(paralogs_summed$label)
  
  # PART II - ANNOTATE ORGANISM TREE
  
  # reshape org_tree_data df
  names(org_data)[1] <- 'label'  # could I have left it as 'id'???
  org_data$label <- as.character(org_data$label)
  org_data[org_data$label == '1', 'label'] <- ''  # WHY???
  
  # annotate org_tree with org_data
  org_data <- org_data[org_data$label %in% as_tibble(org_tree)$label, ]
  org_tree_a <- full_join(org_tree, org_data, by='label')
  
  tab_org_tree <- as_tibble(org_tree_a)
  
  # tt2 -- annotation assignment dataframe: label - paralog1 ... paralogN; counts
  
  paralogs_summed <- paralogs_summed[paralogs_summed$label %in% tab_org_tree$label, ]
  
  tt2 <- tab_org_tree
  tt2 <- full_join(tt2, paralogs_summed, by='label')
  tt2 <- tt2[!is.na(tt2$node), ]
  
  tt2 <- tt2 %>% mutate(across(-1:-6, ~ ifelse(rank == 'species' & is.na(.x), 0, .x)))
  tt2 <- tt2[c(-1, -3:-6)]
  
  org_tree_a <- full_join(org_tree_a, tt2, by='node')
  tab_org_tree <- as_tibble(org_tree_a)
  
  tab_org_tree2 <- traverse_mean(tab_org_tree)
  
  org_tree2 <- full_join(org_tree, tab_org_tree2)
  
  org_tree2
}


copy_annotation <- function(org_tree_prunned, org_tree_full_a) {
  tab_org_tree_full <- as_tibble(org_tree_full_a)
  
  assigniment_tabtree <- tab_org_tree_full[tab_org_tree_full$label %in% as_tibble(org_tree_prunned)$label, -c(1,2,3)]
  
  org_tree_prunned2 <- full_join(org_tree_prunned, assigniment_tabtree, by='label')
  
  org_tree_prunned2
}
