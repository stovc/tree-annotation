library(ggnewscale)

plot_tree <- function(tree, 
                      layout="circular",
                      branch_length=T,
                      tips=NULL,
                      highlight=NA,
                      label=NA,
                      label_out=NA,
                      bars=NULL,
                      color="taxonomy",
                      label_nodes=T,
                      bootstrap=T,
                      tax_start='all',
                      tax_expand=NULL,
                      link=NULL,
                      domains=NULL,
                      context=NULL,
                      circle=NA,
                      label_ancectors=NA,
                      collapse=NA, clade_labels=NA,
                      width=100,
                      height=100,
                      legend="left",
                      filename="tree.pdf",
                      title="") {
  
  COLOR_VECT = c("red", "darkgreen", "blue",              # red green blue
                      "orange", "#808000", "navy",      # orange olive lime
                      "#42d4f4", "#911eb4", "#f032e6",       # cyan purple magenta
                      "#800000", "#469990", "green",          # maroon teal navy
                      "#9A6324", "#eedd00", "#fabed4",       # brown yellow pink
                      "#ffd8b1", "#fffac8", "#aaffc3",       # apricot, beige, mint
                      "#dcbeff", "cornflowerblue", "#bfef45", # lawander --- lime
                      "maroon", "steelblue")                 # lavender, cfb
  
  # title construction constants
  BRANCH_LENGTH_FALSE="-branch_length=F"
  
  # branch length
  if (branch_length == T) {
    branch_length = "branch.length"
  } 
  else {
    title <- paste(title, BRANCH_LENGTH_FALSE)
  }
  
  # reshaping tree
  if (color == "taxonomy") {
    tree <- calibrate_taxonomy(
      tree, 
      start=tax_start, 
      expand=tax_expand, 
      threshold=1)
    
    tab_tree <- as_tibble(tree)  # tabular form of the tree
    
    # construct color list
    tax <- split(tab_tree$node, tab_tree$taxonomy)
    tax <- tax[order(sapply(tax,length),decreasing=T)]
    
    u_tax <- names(tax)
    u_tax <- u_tax[u_tax != "None"]
    u_tax <- u_tax[u_tax != "Other"]
    u_tax <- c(u_tax, "Other", "None")
    
    n_tax <- length(u_tax)
    
    #colors = c(COLOR_VECT, rep("gray", n_tax - 1 - length(COLOR_VECT)), "black")
    named_colors = c(COLOR_VECT[1:(n_tax - 2)], 'gray', 'black')
    
    names(named_colors) = u_tax
  }
  
  if (!is.null(tips)) {
    label  <- 
      print(label)
    
    assign_df <- tab_tree %>% 
      mutate(tip_lab = 
               do.call(paste, c(tab_tree[tips], sep = "|"))
      ) %>% 
      select("node", "tip_lab")
    
    tree <- full_join(tree, assign_df, by="node")
  }
  
  tab_tree <- as_tibble(tree)  # tabular form of the tree
  
  # plot tree
  if (!is.null(link)) {
    xlim = 170
    hjust = -1
  }
  else {
    xlim = NULL
    hjust = 0
  }
  
  p <- ggtree(tree, aes_string(color=color), layout=layout, branch.length=branch_length, xlim=xlim)
  if (color == "taxonomy") {
    p <- p + scale_color_manual(values=named_colors)
    p <- p + guides(colour = guide_legend(ncol = 1)) + 
      guides(colour = guide_legend(ncol = 1))  # guide_legend(override.aes = list(size = 20), ncol = 1)) - like that to make legend elements bigger
    print('MET')
  } else {
    print('NOT MET')
    p <- p + scale_color_manual(values=COLOR_VECT)
  }
  
  # Label tips
  if (!is.null(tips)) {
    p <- p + geom_tiplab(aes(label=tip_lab), hjust=hjust, align=T, size=2)
  }
  
  # Highlight
  if (!is.na(highlight)) {
    data=tab_tree[!is.na(tab_tree[,highlight]), c("node", highlight)]
    p <- p + geom_hilight(data=data, aes(node=node), fill="grey95", color="black", alpha=.5)
  }
  
  # Label
  if (!is.na(label)) {
    p <- p + geom_label2(aes(label=paralog_anc), color="black", fill='grey95')
  }
  
  # Label_out
  if (!is.na(label_out)) {
    data=tab_tree[!is.na(tab_tree[,label_out]), c("node", label_out)]
    p <- p + geom_cladelab(data=data, mapping = aes(node = node, label=paralog_anc), angle="auto")
  }
  
  # Bars
  if (!is.null(bars)) {
    # update title
    str_to_add <- paste0("-bars=[", paste(bars, collapse=", "), "]")
    title <- paste(title, str_to_add)
    
    # iterate collumns to be added as barplots
    for (bar in bars) {
      data=tab_tree[!is.na(tab_tree[,bar]), c("label", bar)]
      colnames(data) <- c("label","value")
      p <- p + geom_fruit(data=data,
                          geom=geom_col,
                          mapping=aes(y=label, x=value),
                          axis.params=list(
                            axis       = "x",
                            text.size  = 1.8,
                            hjust      = 1,
                            vjust      = 0.5,
                            nbreak     = 4,
                            title=bar
                          ),
                          grid.params=list()
      ) + new_scale_fill()
    }
  }
  
  # Bootstrap
  if (bootstrap == T) {
    p <- p + geom_point2(aes(subset=!isTip, 
                             fill=cut(UFboot, c(0, 70, 90, 100))), 
                         shape=21, size=2, colour="black") +
      scale_fill_manual(values=c("white", "grey", "black"), guide='legend', 
                        name='Bootstrap', 
                        breaks=c('(90,100]', '(70,90]', '(0,70]'), 
                        labels=expression(BS>=90,70 <= BS * " < 90", BS < 70))
  }
  
  # Length
  #if (circle == "length") {
  #  p <- p + geom_point2(aes(subset=isTip, fill=length), size=4, colour="black", shape=21) +
  #    scale_fill_viridis_c()
  #}
  
  #p <- p + geom_tippoint(
  #  aes(
  #    color = taxonomy, size=factor(length)
  #    )
  #  )
  #geom_label(aes(label = UFboot, fill = UFboot), size=3, color="black")
  
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
  if (!is.null(link)) {
    # prepare taxa link
    genomes <- split(tab_tree$label, tab_tree[, link])
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
    label_taxon <- tab_tree %>% select(label, taxonomy)
    names(label_taxon)[1] <- "from"
    taxalink_map <- merge(taxalink_map, label_taxon, by = "from")
    
    # add taxalink layer
    p <- p + geom_taxalink(data = taxalink_map, 
                           mapping=aes(taxa1=from, taxa2=to, color=taxonomy), 
                           alpha=0.8, size=0.6, curvature=1, ncp=10, outward=F)
  }
  
  # Domain architecture
  if (!is.null(domains)) {
    p <- p + geom_facet(data = domains,
                        mapping = aes(xmin = start, xmax = end, 
                                      fill = domain),
                        geom = geom_motif,
                        panel = 'Domains',
                        on = '.', label = 'domain', align = 'left',
                        arrowhead_width = grid::unit(0, "mm"),
                        arrowhead_height = grid::unit(3, "mm")) +
      xlim_tree(100)
  }
  
  # genome context
  if (!is.null(context)) {
    p <- p + geom_facet(data = context,
                        mapping = aes(xmin = start, xmax = end, 
                                      fill = motif, forward = strand),
                        geom = geom_motif, 
                        panel = 'Genome context',
                        on = '.', label = 'motif', align = 'left',
                        arrowhead_width = grid::unit(3, "mm"),
                        arrow_body_height = grid::unit(3, "mm"),
                        arrowhead_height = grid::unit(3, "mm")) + 
      xlim_tree(50)
  }
  
  # Legend
  #p <- p + theme(                This works when you have ~5000 tips
  #  legend.key.size = unit(2, "cm"),
  #  legend.key.width = unit(1.5,"cm"),
  #  legend.text = element_text(size = 16), 
  #  legend.position=legend)
  
  p <- p + theme(legend.position=legend)
  
  # Title
  p <- p + ggtitle(title)
  
  if (is.null(domains) & is.null(context)) {
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
                          width='auto', 
                          height='auto',
                          legend='right',
                          filename='org_tree.svg') {
  
  RANKS = c('superkingdom', 'phylum', "class", "order", "family", "genus", "species")
  
  # prepare heatmap data
  tab_tree <- as_tibble(tree)
  
  data_heatmap <- tibble(
    tab_tree[, c(-3, -5:-6)]
  )  # parent, node, label, c(paralogs)
  
  data_heatmap <- data_heatmap[!(data_heatmap$node %in% tab_tree$parent), ]
  
  data_heatmap <- data_heatmap[, c(-1, -2)]
  
  data_heatmap <- pivot_longer(data_heatmap, -1, names_to='paralog')
  data_heatmap <- data_heatmap[is.na(data_heatmap$value) == F,]
  
  data_heatmap$x = unclass(factor(data_heatmap$paralog, levels=levels))
  mapping <- data_heatmap[duplicated(data_heatmap$paralog) == F, ]
  
  # plot tree
  p <- ggtree(tree, aes_string(color='rank'), layout='rectangular', branch.length='none')
  
  y <- max(p$data$y)  # no clue what it is...
  
  # add annotations
  p <- p + geom_label2(aes(subset=!isTip, label=name, size=rank), fill='grey95') +
    geom_facet(data = data_heatmap, geom = geom_tile,
               mapping = aes(x=x, fill=value),
               colour='black',
               panel = 'Paralogs') +
    geom_facet(data = data_heatmap, geom = geom_text,
               mapping = aes(x=x, label=round(value,1)),
               colour='white',
               panel = 'Paralogs') +
    geom_facet(data = mapping, geom = geom_text,
               mapping = aes(x=x, y=0, label=paralog),
               colour='black',
               panel = 'Paralogs')
  
  p <- p + geom_tiplab(aes(label=name)) + xlim_tree(6)
  
  # apply scales
  size_scale = c(14, 11, 9, 7, 5, 3, 2)
  names(size_scale) <- RANKS
  
  color_scale = c("darkred", "darkorange", "yellow4", "darkgreen", "darkcyan", "darkblue", "darkviolet")
  names(color_scale) <- c('superkingdom', 'phylum', "class", "order", "family", "genus", "species")
  
  p <- p +
    scale_fill_viridis_c(option="turbo") +
    scale_size_manual(values=size_scale) +
    scale_color_manual(values=color_scale)
  
  # Apply theme
  p <- p + theme(
    legend.key.size = unit(1, "cm"),
    legend.key.width = unit(1, "cm"),
    legend.text = element_text(size = 12), 
    legend.position=legend) + scale_x_ggtree()
  
  # calculate height and width
  if (width == 'auto') {
    n_ranks = length(intersect(unique(tab_tree$rank), RANKS))
    n_paralogs = nrow(mapping)
    
    tree_width = n_ranks * 22
    heatmap_width = n_paralogs * 2
    
    width = tree_width + heatmap_width
  }
  
  if (height == 'auto') {
    n_tips = sum(isTip(tab_tree, tab_tree$node))
    height = n_tips
  }
  
  # print and save
  facet_widths(p, widths=c(tree_width, heatmap_width))
  ggsave(filename, width = width, height = height, units = "cm", limitsize = F)
  
}
