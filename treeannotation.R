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

annotate_tree <- function(tree, start, expand, threshold) {
  
  tree <- assign_tax(tree, start, expand)
  
  tree <- filter_taxonomy(tree, threshold)
  tree <- filter_genomes(tree)
  
  tree
}


annotate_org_tree <- function(org_tree, org_data, paralogs) {
  names(org_data)[1] <- 'label'  # could I have left it as 'id'???
  org_data$label <- as.character(org_data$label)
  org_data[org_data$label == '1', 'label'] <- ''  # WHY???
  
  # org_tree annotated with org_data
  org_data <- org_data[org_data$label %in% as_tibble(org_tree)$label, ]
  org_tree_a <- full_join(org_tree, org_data, by='label')
  
  tab_org_tree <- as_tibble(org_tree_a)
  
  # tt2 -- annotation assignment dataframe: label - paralog1 ... paralogN; counts
  View(org_data)
  View(tab_org_tree)
  View(paralogs)
  print(class(tab_org_tree$label))
  print(class(paralogs$label))
  
  paralogs <- paralogs[paralogs$label %in% tab_org_tree$label, ]
  
  tt2 <- tab_org_tree
  tt2 <- full_join(tt2, paralogs, by='label')
  tt2 <- tt2[!is.na(tt2$node), ]
  tt2[is.na(tt2)] <- 0
  
  tt2 <- tt2 %>% mutate(across(-1:-6, ~ ifelse(rank == 'species' & is.na(.x), 0, .x)))
  tt2 <- tt2[c(-1, -3:-6)]
  
  View(tt2)
  
  org_tree_a <- full_join(org_tree_a, tt2, by='node')
  tab_org_tree <- as_tibble(org_tree_a)
  
  tab_org_tree2 <- traverse_mean(tab_org_tree)
  
  org_tree2 <- full_join(org_tree, tab_org_tree2)
  
  org_tree2
}


copy_annotation <- function(org_tree_prunned, org_tree_full_a) {
  tab_org_tree_full <- as_tibble(org_tree_full_a)
  View(tab_org_tree_full)
  
  View(as_tibble(org_tree_prunned))
  
  assigniment_tabtree <- tab_org_tree_full[tab_org_tree_full$label %in% as_tibble(org_tree_prunned)$label, -c(1,2,3)]
  
  View(assigniment_tabtree)
  
  org_tree_prunned2 <- full_join(org_tree_prunned, assigniment_tabtree, by='label')
  
  View(as_tibble(org_tree_prunned2))
  
  org_tree_prunned2
}


assign_full_taxonomy <- function(tree, csv) {
  tree <- assign_taxonomy(tree, csv, "phylum")
  tree <- assign_taxonomy(tree, csv, "class")
  tree <- assign_taxonomy(tree, csv, "order")
  tree <- assign_taxonomy(tree, csv, "family")
  tree <- assign_taxonomy(tree, csv, "genus")
  tree <- assign_taxonomy(tree, csv, "species")
  
  tree <- assign_taxonomy(tree, csv, "product")
  tree <- assign_taxonomy(tree, csv, "gene")
  
  tree <- assign_taxonomy(tree, csv, "assembly", recursive = F)
  tree <- assign_taxonomy(tree, csv, "length", recursive = F)
  tree <- assign_taxonomy(tree, csv, "taxid", recursive = F)
  
  tree
}

assign_paralogs <- function(tree, protein, paralod_df) {
  #' assign paralog annotation to all descendants of nodes from dataframe
  #' nonspecified nodes are assigned from from the 'protein' argument
  
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


assign_tax <- function(tree, start, tax_vect) {
  
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

assign_taxonomy <- function(tree, df, taxon, recursive = T) {
  
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
    mutate(filt_tax = ifelse(n() > threshold, taxon, "Other")) %>% 
    ungroup()
  
  assignment_tree <- tab_tree
  
  assignment_tree <- assignment_tree %>%
    select(node, filt_tax)
  
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
    mutate(filt_tax = ifelse(n() > threshold, taxon, "Other")) %>% 
    ungroup()
  
  assignment_tree <- tab_tree
  
  assignment_tree <- assignment_tree %>%
    select(node, filt_tax)
  
  tree <- full_join(tree, assignment_tree, by = "node")
  tree
}


read_context_data <- function() {
  context_data <- read_delim("genome_context.csv", delim=';')
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


traverse_properties <- function(tree, csv) {
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



##########################################

plot_tree <- function(tree, layout, branch_length, aes_color="filt_tax",
                      taxalink=F, tips=NA, label_nodes=F, circle=NA, label_ancectors=NA, 
                      collapse=NA, clade_labels=NA, 
                      domains=NA, context=NA,
                      width=100, height=100, legend="left",
                      filename="tree.pdf") {
  
  COLOR_VECT = c("red", "#3cb44b", "blue",              # red green blue
                 "orange", "#808000", "lawngreen",      # orange olive lime
                 "#42d4f4", "#911eb4", "#f032e6",       # cyan purple magenta
                 "#800000", "#469990", "navy",          # maroon teal navy
                 "#9A6324", "#eedd00", "#fabed4",       # brown yellow pink
                 "#ffd8b1", "#fffac8", "#aaffc3",       # apricot, beige, mint
                 "#dcbeff", "cornflowerblue", "#bfef45", # lawander --- lime
                 "maroon", "steelblue")                 # lavender, cfb
  
  tab_tree <- as_tibble(tree)  # tabular form of the tree
  
  # construct color list
  tax <- split(tab_tree$node, tab_tree$filt_tax)
  tax <- tax[order(sapply(tax,length),decreasing=T)]
  
  u_tax <- names(tax)
  u_tax <- u_tax[u_tax != "None"]
  u_tax <- u_tax[u_tax != "Other"]
  u_tax <- c(u_tax, "Other", "None")
  
  n_tax <- length(u_tax)
  
  #colors = c(COLOR_VECT, rep("gray", n_tax - 1 - length(COLOR_VECT)), "black")
  colors = c(COLOR_VECT[1:(n_tax - 1)], 'black')
  colors = c(COLOR_VECT[1:(n_tax - 2)], 'gray', 'black')
  names(colors) = u_tax
  
  # plot tree
  if (taxalink == T) {
    xlim = 150
  }
  else {
    xlim = NULL
  }
  
  p <- ggtree(tree, aes_string(color=aes_color), layout=layout, branch.length=branch_length, xlim=xlim)
  if (aes_color == "filt_tax") {
    p <- p + scale_color_manual(values=colors)
    p <- p + guides(colour = guide_legend(ncol = 1)) + 
      guides(colour = guide_legend(ncol = 1))  # guide_legend(override.aes = list(size = 20), ncol = 1)) - like that to make legend elements bigger
    print('MET')
  } else {
    print('NOT MET')
  }
  
  # Bootstraps
  if (circle == "bootstrap") {
    p <- p + geom_point2(aes(subset=!isTip, 
                             fill=cut(UFboot, c(0, 70, 90, 100))), 
                         shape=21, size=2, colour="black") +
      scale_fill_manual(values=c("white", "grey", "black"), guide='legend', 
                        name='Bootstrap', 
                        breaks=c('(90,100]', '(70,90]', '(0,70]'), 
                        labels=expression(BS>=90,70 <= BS * " < 90", BS < 70))
  }
  
  # Length
  if (circle == "length") {
    p <- p + geom_point2(aes(subset=isTip, fill=length), size=4, colour="black", shape=21) +
      scale_fill_viridis_c()
  }
  
  #p <- p + geom_tippoint(
  #  aes(
  #    color = filt_tax, size=factor(length)
  #    )
  #  )
  #geom_label(aes(label = UFboot, fill = UFboot), size=3, color="black")
  
  # Label tips
  if (!is.na(tips)) {
    #p <- p + geom_tiplab(aes_string(label=tips), color='black', angle=0, align=TRUE, linesize=.1)
    p <- p + geom_tiplab(aes_string(label=tips), hjust=-.1, size=3, align=T)
  }
  
  # Label nodes
  if (label_nodes == T) {
    p <- p + geom_text2(aes(subset=!isTip, label=node), color='black', hjust=-.1, size=1)
    #p <- p + geom_tiplab(hjust=-.1, size=1)
  }
  
  # Collapse
  if (!is.na(collapse)) {
    tax_ancestors <- tab_tree[!is.na(tab_tree$phylum_anc) & tab_tree$phylum_anc != collapse, ]
    
    tax_ancestors <- tax_ancestors %>%
      select(node, phylum)
    
    for(i in tax_ancestors[["node"]]) {
      #tax <- tax_ancestors %>% 
      #  filter(node == i) %>% 
      #  select(phylum) %>% 
      #  pull(phylum)
      print(i)
      #color = colors[tax]
      #print(color)
      p <- collapse(p, i)
    }
  }
  
  # Clade label tax ancestors
  if (!is.na(label_ancectors)) {
    column = paste(label_ancectors, 'anc', sep='_')
    # make a dataframe of taxon ancestors
    View(tab_tree)
    tax_ancestors <- tab_tree[!is.na(tab_tree[[column]]),]
    
    tax_ancestors <- tax_ancestors %>%
      select(node, !!label_ancectors)
    # & phylum != 'Actinobacteria'
    View(tax_ancestors)
    
    for(i in tax_ancestors[["node"]]) {
      tax <- tax_ancestors %>% 
        filter(node == i) %>% 
        select(!!label_ancectors) %>% 
        pull(!!label_ancectors)
      print(tax)
      #color = colors[tax]
      #print(color)
      p <- p + geom_cladelabel(node=i, label=tax, align=TRUE, fontsize=4)
    }
    # geom_tiplab(size=1)
  }
  
  print('geom cladelabel done')
  
  # Clade label from clade_labels
  if (!is.na(clade_labels)) {
    for(i in node) {
      print(i)
      print(clade_labels[clade_labels$node == i, 'label'])
      print(clade_labels[clade_labels$node == i, 'angle'])
      p <- p + geom_cladelabel(i, 
                               label=clade_labels[clade_labels$node == i, 'label'], 
                               angle=clade_labels[clade_labels$node == i, 'angle'],
                               fontsize=16)
    }
  }
  
  # Taxalink
  if (taxalink == T) {
    # prepare taxa link
    genomes <- split(tab_tree$label, tab_tree$filt_genome)
    genome_names <- names(genomes)
    
    from = c()
    to = c()
    distance = c()
    
    for (i in genomes) {
      if (length(i) > 1) {
        print(i)
        comb <- combn(i, 2)
        from <- c(from, comb[1,])
        to <- c(to, comb[2,])
        #distance <- tab_tree[]
      }
    }
    
    taxalink_map <- data.frame(from, to)
    label_taxon <- tab_tree %>% select(label, filt_tax)
    names(label_taxon)[1] <- "from"
    taxalink_map <- merge(taxalink_map, label_taxon, by = "from")
    
    # add taxalink layer
    p <- p + geom_taxalink(data = taxalink_map, 
                           mapping=aes(taxa1=from, taxa2=to, color=filt_tax), 
                           alpha=0.2, size=0.5, curvature=1, ncp=10) + 
      coord_cartesian(clip = 'off') + 
      theme_tree2(plot.margin=unit(c(1,70,1,1), "cm"))
  }
  
  # Label tips
  if (!is.na(tips)) {
    p <- p + geom_tiplab(aes_string(label=tips), size=4, color='black', angle=0, align=TRUE, linesize=.1)
  }
  
  # Domain architecture
  if (!is.na(domains)) {
    p <- p + geom_facet(mapping = aes(xmin = start, xmax = end, fill = domain),
                        data = domain_data,
                        geom = geom_motif,
                        panel = 'Domains',
                        on = domains, label = 'domain', align = 'left',
                        arrowhead_width = grid::unit(0, "mm"),
                        arrow_body_height = grid::unit(3, "mm"),
                        arrowhead_height = grid::unit(3, "mm")) + 
      coord_cartesian(clip = 'off') + 
      theme_tree2(plot.margin=unit(c(1,70,1,1), "cm"))
  }
  
  # genomic context
  if (!is.na(context)) {
    p <- p + geom_facet(data = context_data,
                        mapping = aes(xmin = start, xmax = end, 
                                      fill = motif, forward = strand),
                        geom = geom_motif, 
                        panel = 'Genomic context',
                        on = context, label = 'motif', align = 'left',
                        arrowhead_width = grid::unit(3, "mm"),
                        arrow_body_height = grid::unit(3, "mm"),
                        arrowhead_height = grid::unit(3, "mm")) + 
      coord_cartesian(clip = 'off') + 
      theme_tree2(plot.margin=unit(c(1,70,1,1), "cm"))
  }
  
  # Apply theme
  p <- p + theme(
    legend.key.size = unit(2, "cm"),
    legend.key.width = unit(1.5,"cm"),
    legend.text = element_text(size = 48), 
    legend.position=legend)
  
  if (is.na(domains) & is.na(context)) {
    print('printing plot')
    print(p)
  }
  else {
    print('printing facets')
    facet_widths(p, widths=c(1,3))
  }
  
  print('printed')
  ggsave(filename, width = width, height = height, units = "cm", limitsize = F)
  
  print('saved')
}




plot_org_tree <- function(tree, levels,
                          width=100, height=100, legend='right',
                          filename='org_tree.svg') {
  # prepare heatmap data
  tab_tree <- as_tibble(tree)
  
  data_heatmap <- tibble(
    tab_tree[, c(-3, -5:-6)]
  )  # parent, node, label, c(paralogs)
  
  data_heatmap <- data_heatmap[!(data_heatmap$node %in% tab_tree$parent), ]
  
  data_heatmap <- data_heatmap[, c(-1, -2)]
  
  data_heatmap <- pivot_longer(data_heatmap, -1, names_to='paralog')
  data_heatmap <- data_heatmap[is.na(data_heatmap$value) == F,]
  
  # plot tree
  p <- ggtree(tree, aes_string(color='rank'), layout='rectangular', branch.length='none')
  
  print('constructed')
  
  # Apply theme

  data_heatmap$x = unclass(factor(data_heatmap$paralog, levels=levels))
  mapping <- data_heatmap[duplicated(data_heatmap$paralog) == F, ]
  
  y <- max(p$data$y)
  
  p <- p + geom_text2(aes(label=paste(name)), size=4) +
    geom_facet(data = data_heatmap, geom = geom_tile,
               mapping = aes(x=x, fill=value),
               colour='black',
               panel = 'Paralogs') +
    geom_facet(data = data_heatmap, geom = geom_text,
               mapping = aes(x=x, label=round(value,1)),
               colour='white', size=4,
               panel = 'Paralogs') +
    geom_facet(data = mapping, geom = geom_text,
               mapping = aes(x=x, y=0, label=paralog),
               colour='black', size=4,
               panel = 'Paralogs')
  #p <- p + geom_tiplab(label='name')

  p <- p + theme(
    legend.key.size = unit(1, "cm"),
    legend.key.width = unit(1, "cm"),
    legend.text = element_text(size = 12), 
    legend.position=legend) + scale_x_ggtree()
    
  print('themed')
  
  print('heatmapped')
  
  print(p)
  
  print('printed')
  
  ggsave(filename, width = width, height = height, units = "cm", limitsize = F)
  
  print('saved')
}
