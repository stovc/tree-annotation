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
                           alpha=0.4, size=0.3, curvature=1, ncp=10, outward=F)
  }
  
  # Domain architecture
  if (!is.na(domains)) {
    p <- p + geom_facet(data = domain_data,
                        mapping = aes(xmin = start, xmax = end, 
                                      fill = domain),
                        geom = geom_motif,
                        panel = 'Domains',
                        on = domains, label = 'domain', align = 'left',
                        arrowhead_width = grid::unit(0, "mm"),
                        arrowhead_height = grid::unit(3, "mm")) +
      xlim_tree(100)
  }
  
  # genome context
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
      xlim_tree(100)
  }
  
  # Legend
  #p <- p + theme(                This works when you have ~5000 tips
  #  legend.key.size = unit(2, "cm"),
  #  legend.key.width = unit(1.5,"cm"),
  #  legend.text = element_text(size = 16), 
  #  legend.position=legend)
  
  p <- p + theme(legend.position=legend)
  
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
  
  # geom_text2(aes(label=paste(name)), size=4)
  
  size_scale = c(18, 14, 10, 7, 5, 3, 2)
  names(size_scale) <- c('superkingdom', 'phylum', "class", "order", "family", "genus", "species")
  
  color_scale = c("darkred", "darkorange", "yellow4", "darkgreen", "darkcyan", "darkblue", "darkviolet")
  names(color_scale) <- c('superkingdom', 'phylum', "class", "order", "family", "genus", "species")
  
  p <- p + geom_label2(aes(subset=!isTip, label=name, size=rank), fill='grey95') +
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
               panel = 'Paralogs') +
    scale_fill_viridis_c(option="turbo") +
    scale_size_manual(values=size_scale) +
    scale_color_manual(values=color_scale)
  
  p <- p + geom_tiplab(aes(label=name)) + xlim_tree(6)
  
  p <- p + theme(
    legend.key.size = unit(1, "cm"),
    legend.key.width = unit(1, "cm"),
    legend.text = element_text(size = 12), 
    legend.position=legend) + scale_x_ggtree()
  
  print('themed')
  
  print('heatmapped')
  
  facet_widths(p, widths=c(4,1))
  
  print('printed')
  
  ggsave(filename, width = width, height = height, units = "cm", limitsize = F)
  
  print('saved')
}
