# MUTATIONS BY LOCUS / VCS ANALYSIS ############################################
################################################################################

# helper function to count occurrences of the target character in a list 
count_character <- function(list_element, target) {
  sum(unlist(lapply(list_element, function(x) sum(x == target))))
}


# extraction function 
extract_mutations <- function(df, gene_length, substitution_column){
  print(substitution_column)
  # Initialize lists to store unique mutations per locus
  mutations_list <- vector('list', gene_length)
  mutations_list_unique <- vector('list', gene_length)
  
  
  # Iterate through the substitutions column
  # for .aa files it would 'aa_mutationss'
  # for .dna files it would be 'substitutions'
  for (i in 1:nrow(df)) {
    subs <- strsplit(as.character(df[[substitution_column]][i]), " ")[[1]]
    
    for (sub in subs) {
      if (!is.na(sub) && sub != ""){
        # extract numeric elements from sub (the locus)
        locus <- as.numeric(gsub("[^0-9]", "", sub))
        if (!is.na(locus) && locus >= 1 && locus <= gene_length) {
          mutation <- substr(sub, nchar(sub), nchar(sub))
          mutations_list_unique[[locus]] <- unique(c(mutations_list[[locus]], mutation))
          mutations_list[[locus]] <- c(mutations_list[[locus]], mutation)
        }
      }
    }
  }
  # Count the number of mutations and create a string of unique mutations for each locus
  mutation_count <- sapply(mutations_list, length)
  earlystop_count <- sapply(mutations_list, count_character, target = '*')
  mutation_count_unique <- sapply(mutations_list_unique, length)
  unique_mutations <- sapply(mutations_list, function(x) paste(unique(x), collapse = " "))
  
  return(list(mutation_count = mutation_count, earlystop_count = earlystop_count,
              mutation_count_unique = mutation_count_unique,
              unique_mutations = unique_mutations, 
              mutations_list = mutations_list))
  
}



# parsing function - put it all together, return data frames
# run separate for codons and aa ('aa_substitions' vs 'substitutions' column) 
count_mutations <- function(df, gene_length, substitution_column, idx1=NULL, idx2=NULL) {
  print(substitution_column)
  if (substitution_column == 'aa_mutations'){
    
    # Ignore indels and wt for mutational saturation purposes
    df <- df %>% filter(number_of_indels == 0 & aa_mutations != 'wt' & aa_mutations != 'synonymous' & n_aa_substitutions > 0)
    
    mutation_count_aa <- extract_mutations(df, gene_length, 'aa_mutations')
    
    data.frame(
      locus = 1:gene_length,
      mutation_count_aa = mutation_count_aa$mutation_count,
      earlystop_count = mutation_count_aa$earlystop_count,
      mutation_count_unique_aa = mutation_count_aa$mutation_count_unique,
      unique_mutations_aa = mutation_count_aa$unique_mutations, 
      mutations_list_aa = I(mutation_count_aa$mutations_list))
    
    
    # Legacy from parsing the CONSENSUS output directly from alignparse
    # in the .dna tables the column is 'substitutions'
  } else if (substitution_column == 'substitutions') {
    mutation_count_codons <- extract_mutations(df, gene_length, 'substitutions')
    
    x <- data.frame(
      locus = 1:gene_length,
      mutation_count_codons = mutation_count_codons$mutation_count,
      mutation_count_unique_codons = mutation_count_codons$mutation_count_unique,
      unique_mutations_codons = mutation_count_codons$unique_mutations,
      mutations_list_codons = I(mutation_count_codons$mutations_list))
    
    # masking primer binding regions if they exist (index provided in targets manifest)
    if (!is.null(idx1) & !is.null(idx2)) {
      x <- x[idx1:idx2, ]
    }
    
    x
  }
  
}


# VARIANT CALL SUPPORT ANALYSIS DNA ############################################


# does not keep mutbylocus (will do so in aa version)
vcs_analysis <- function(min_vcs, targets_manifest){
  numbarcode <- list()
  outputs <- list()
  
  for (i in 1:nrow(targets_manifest)) {
    target <- get(targets_manifest$dataframe[[i]])
    target <- subset(target, variant_call_support >= min_vcs)
    numbarcode[[i]] <- nrow(target)
    
    mutbylocus <- paste0(targets_manifest$dataframe[[i]], '_mutbylocusDNA_',min_vcs)
    
    
    # count mutations by locus
    # if there are primer binding masking region indicated by maskidx1, maskidx2,
    # count_mutations will return the sliced by locus dataframe
    if (is.na(targets_manifest$maskidx1[[i]])) {
      outputs[[mutbylocus]] <- count_mutations(target, targets_manifest$Length[[i]], 'substitutions')
    } else {
      outputs[[mutbylocus]] <- count_mutations(target, targets_manifest$Length[[i]], 'substitutions', 
                                               targets_manifest$maskidx1[[i]], 
                                               targets_manifest$maskidx2[[i]])
    }
  }
  
  
  Percentloci <- sapply(outputs, function(x) (sum(x[["mutation_count_codons"]] !=0)/nrow(x))*100)
  MeanNucSubsPerLocus <- sapply(outputs, function(x) mean(x[["mutation_count_codons"]]))
  MeanNucSubsPerLocus_Unique <- sapply(outputs, function(x) mean(x[["mutation_count_unique_codons"]]))
  MeanMutSaturation <- sapply(outputs, function(x) mean(x[["mutation_count_unique_codons"]]/3))
  
  data.frame(target = targets_manifest$dataframe, NumBarcodes = unlist(numbarcode),
             Percentloci = Percentloci,
             MeanNucSubsPerLocus = MeanNucSubsPerLocus,
             MeanNucSubsPerLocus_Unique = MeanNucSubsPerLocus_Unique,
             MeanMutSaturation = MeanMutSaturation,
             VCS_thresh = rep(min_vcs, nrow(targets_manifest)))
  
}

################################################################################

# VISUALIZATION
integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}





vcs_analysis_aa <- function(min_vcs, targets_manifest){
  numbarcode <- list()
  outputs <- list()
  
  for (i in 1:nrow(targets_manifest)) {
    target <- get(targets_manifest$dataframe[[i]])
    target <- subset(target, variant_call_support >= min_vcs)
    numbarcode[[i]] <- nrow(target)
    
    mutbylocus <- paste0(targets_manifest$dataframe[[i]], '_mutbylocusAA_',min_vcs)
    mutbylocus_out <- count_mutations(target, targets_manifest$Length[[i]]/3, 'aa_substitutions')
    
    # Fix the self-counting for non 'M' amino acids'
    mutspace_df <- get(targets_manifest$mutspace[[i]])
    mutspace_df$num_possible_mutations <- ifelse(mutspace_df$original_aa != 'M',
                                                 mutspace_df$num_possible_mutations - 1, 
                                                 mutspace_df$num_possible_mutations)
    
    
    mutbylocus_out <- merge(mutbylocus_out, mutspace_df, all.x=TRUE)
    
    outputs[[mutbylocus]] <- mutbylocus_out
    
  }
  
  Percentloci <- sapply(outputs, function(x) (sum(x[["mutation_count_aa"]] !=0)/nrow(x))*100)
  MeanAASubsPerLocus <- sapply(outputs, function(x) mean(x[["mutation_count_aa"]]))
  MeanAASubsPerLocus_Unique <- sapply(outputs, function(x) mean(x[["mutation_count_unique_aa"]]))
  MeanMutSaturationAA <- sapply(outputs, function(x) mean(x[["mutation_count_unique_aa"]]/x[["num_possible_mutations"]]))
  
  data.frame(target = targets_manifest$dataframe, NumBarcodes = unlist(numbarcode),
             Percentloci = Percentloci,
             MeanAASubsPerLocus = MeanAASubsPerLocus,
             MeanAASubsPerLocus_Unique = MeanAASubsPerLocus_Unique,
             MeanMutSaturationAA = MeanMutSaturationAA,
             VCS_thresh = rep(min_vcs, nrow(targets_manifest)))
  
}


# Function to compute Shannon diversity index (H) and equitability (E)
compute_shannon_equitability <- function(mutations_list, num_possible_mutations) {
  # Check if num_possible_mutations is positive
  if (num_possible_mutations <= 0 || is.na(num_possible_mutations)) {
    return(list(H = NA, E = NA))
  }
  
  # Compute counts of each unique mutation
  counts <- table(mutations_list)
  
  # Total number of mutations
  N <- sum(counts)
  
  # Compute probabilities (p_i)
  p_i <- counts / N
  
  # Compute Shannon diversity index (H) using log2
  H <- -sum(p_i * log2(p_i))
  
  # Compute maximum possible Shannon diversity index (H_max)
  H_max <- log2(num_possible_mutations)
  
  # Compute Shannon equitability metric (E)
  E <- H / H_max
  
  # Return H and E
  return(list(H = H, E = E))
}



# 
# fix_mutspace <- function(df) {
#   library(stringr)
#   df <- df %>%
#     rename(possible_mutations = any_of('possible_aa_mutations')) %>%
#     mutate(
#       # Remove the original_aa from possible_mutations
#       possible_mutations = str_remove_all(possible_mutations, paste0("\\b", original_aa, "\\b")),
#       # Clean up any extra whitespace
#       possible_mutations = str_trim(str_squish(possible_mutations)),
#       # Count the number of possible mutations
#       num_possible_mutations = ifelse(possible_mutations == "", 0, str_count(possible_mutations, "\\S+"))
#     ) %>%
#     dplyr::select(-any_of('num_possible_aa_mutations'))
#   return(df)
# }


fix_mutspace <- function(df) {
  # If df might be an S4 DataFrame, coerce to tibble to keep dplyr semantics
  df <- tibble::as_tibble(df)
  
  df <- df %>%
    dplyr::rename_with(~ "possible_mutations", tidyselect::any_of("possible_aa_mutations")) %>%
    dplyr::mutate(
      possible_mutations = stringr::str_remove_all(.data$possible_mutations,
                                                   paste0("\\b", .data$original_aa, "\\b")),
      possible_mutations = stringr::str_squish(possible_mutations),
      num_possible_mutations = dplyr::if_else(
        is.na(possible_mutations) | possible_mutations == "", 0L,
        stringr::str_count(possible_mutations, "\\S+")
      )
    ) %>%
    dplyr::select(-tidyselect::any_of("num_possible_aa_mutations"))
  
  df
}


# This is to get mut by locus
parse_mutations <- function(min_vcs, targets_manifest){
  numbarcode <- list()
  outputs <- list()
  
  for (i in 1:nrow(targets_manifest)) {
    print(paste("Processing:", targets_manifest$dataframe[[i]]))
    
    target <- get(targets_manifest$dataframe[[i]])
    target <- subset(target, variant_call_support >= min_vcs)
    numbarcode[[i]] <- nrow(target)
    print(numbarcode[[i]])
    
    mutbylocus <- paste0(targets_manifest$dataframe[[i]], '_mutbylocusAA_',min_vcs)
    mutbylocus_out <- count_mutations(target, targets_manifest$Length[[i]]/3, 'aa_mutations')
    print(head(mutbylocus_out))
    
    # Fix the self-counting for non 'M' amino acids'
    mutspace_df <- get(targets_manifest$mutspace[[i]])
    
    # NO longer needed with PREPROCESSED VARDICTS (whitelisted barcodes)
    # Or alternatively, after applying fix_mutations()
    # mutspace_df$num_possible_mutations <- ifelse(mutspace_df$original_aa != 'M',
    #                                              mutspace_df$num_possible_mutations - 1, 
    #                                              mutspace_df$num_possible_mutations)
    
    
    mutbylocus_out <- merge(mutbylocus_out, mutspace_df, all.x=TRUE)
    
    
    # Calculate the Shannon equitability metric
    shannon <- mapply(
      compute_shannon_equitability,
      mutbylocus_out$mutations_list,
      mutbylocus_out$num_possible_mutations,
      SIMPLIFY = FALSE
    )
    
    mutbylocus_out$H <- sapply(shannon, function(x) x$H)
    mutbylocus_out$E <- sapply(shannon, function(x) x$E)
    
    print(head(mutbylocus_out))  # Debugging
    
    outputs[[mutbylocus]] <- mutbylocus_out
    
  }
  outputs
}



# for plotting by gene segment, helper function to add GeneSegment column to mutbylocus.
# Requires vector of break points eg. (1, 500, 700). Left inclusive right exclusive

# add_gene_segment <- function(df, breakpoints) {
#   # Ensure breakpoints are sorted
#   breakpoints <- sort(breakpoints)
#   
#   # Create labels for segments
#   labels <- paste0("Seg", seq_along(breakpoints) + 1)
#   
#   # Add labels for below the first breakpoint
#   labels <- c("Seg1", labels)
#   
#   # Use case_when to create the new GeneSegment column
#   df <- df %>%
#     mutate(GeneSegment = case_when(
#       locus < breakpoints[1] ~ "Seg1",
#       locus >= breakpoints[length(breakpoints)] ~ paste0("Seg", length(breakpoints) + 1),
#       TRUE ~ labels[findInterval(locus, breakpoints)]
#     ))
#   
#   return(df)
# }


add_gene_segment <- function(df, breakpoints) {
  # Ensure breakpoints are sorted
  breakpoints <- sort(breakpoints)
  
  # Create labels for segments
  labels <- paste0("Seg", seq_len(length(breakpoints)))
  
  
  # Use findInterval to create the new GeneSegment column
  df <- df %>%
    mutate(GeneSegment = labels[findInterval(locus, breakpoints, rightmost.closed = TRUE)])
  
  return(df)
}




# plotting function for mut by locus
plot_mutbylocus_all <- function(df, save_path, gene, segment.breaks, num.mut.max, mut.sat.max) {
  
  print(paste0('Working on gene ', gene))
  
  df_mutbylocus <- df %>% pivot_longer(cols = mutation_count_aa:earlystop_count,
                                       names_to = 'count_var',
                                       values_to = 'count')
  print(colnames(df_mutbylocus))
  
  p_mutbylocus <- ggplot(df_mutbylocus, aes(x=locus, y=count, color=count_var))+
    geom_point()+
    geom_smooth(method='loess')+
    theme_classic()+
    scale_color_brewer(palette = 'Set2')+
    theme(axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          axis.text.y = element_text(size=15),
          axis.text.x = element_text(size=15),
          panel.border = element_rect(fill=NA, colour="black", linewidth=1),
          legend.position = 'none')+
    scale_y_continuous(trans='log10', limits=c(NA, num.mut.max))+
    ylab('Number of Mutations')+
    xlab('AA Locus')+
    ggtitle(gene)

  
  
  
  df_mutsat <- df %>% 
    mutate(MutSat = (mutation_count_unique_aa/(num_possible_mutations))*100)
  
  df_mutsat <- add_gene_segment(df_mutsat, segment.breaks)
  
  print(paste0(gene, " Number of gene segments: ", length(unique(df_mutsat$GeneSegment))))
  print(paste0(gene, " Segment breaks: ", segment.breaks))
  
  print("Head of df_mutsat")
  print(head(df_mutsat))
  
  p_mutsat <- ggplot(df_mutsat, aes(x=locus, y=MutSat))+
    geom_hline(yintercept = 100, color='lightgrey')+
    geom_hline(yintercept = 80, color='lightgrey')+
    geom_line(linewidth = 0.7)+
    scale_y_continuous(limits=c(0, mut.sat.max))+
    theme_classic()+
    theme(axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          axis.text.y = element_text(size=15),
          axis.text.x = element_text(size=15),
          panel.border = element_rect(fill=NA, colour="black", linewidth=1))+
    ylab('Mutational Space Saturation (%)')+
    xlab('AA Locus')+
    ggtitle(gene)
  
  
  p_mutsat_hist <- ggplot(df_mutsat, aes(x=MutSat))+
    geom_histogram(binwidth = 1, color = 'black') + 
    #geom_density(color='black')+
    scale_y_continuous(trans='log10', limits=c(NA,1000))+
    scale_x_continuous(limits=c(0,mut.sat.max))+
    theme_classic()+
    theme(axis.title.x = element_text(size=13),
          axis.title.y = element_text(size=16),
          axis.text.y = element_text(size=15),
          axis.text.x = element_text(size=15),
          panel.border = element_rect(fill=NA, colour="black", linewidth=1))+
    xlab('Mutational Space Saturation (%)')+
    ylab('Number of Loci')+
    ggtitle(gene)
  
  p_mutsat_dens <- ggplot(df_mutsat, aes(x=MutSat))+
    #geom_histogram(binwidth = 1, color = 'black') + 
    geom_density(color='black')+
    #scale_y_continuous(trans='log10', limits=c(NA,1000))+
    scale_x_continuous(limits=c(0,mut.sat.max))+
    theme_classic()+
    theme(axis.title.x = element_text(size=13),
          axis.title.y = element_text(size=16),
          axis.text.y = element_text(size=15),
          axis.text.x = element_text(size=15),
          panel.border = element_rect(fill=NA, colour="black", linewidth=1))+
    xlab('Mutational Space Saturation (%)')+
    ylab('Density')+
    ggtitle(gene)
  
  
  p_mutsat_seg <- ggplot(df_mutsat, aes(x = GeneSegment, y = MutSat))+
    geom_boxplot(outlier.shape = NA, color = 'black') + 
    geom_jitter(width = 0.3, color = 'black')+
    theme_classic()+
    scale_y_continuous(limits=c(0, mut.sat.max))+
    theme(axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15),
          axis.text.y = element_text(size=15),
          axis.text.x = element_text(size=12),
          panel.border = element_rect(fill=NA, colour="black", linewidth=1))+
    xlab('Gene Segment')+
    ylab('Mutational Space Saturation (%)')+
    ggtitle(gene)
  
  p_combo <- ggarrange(p_mutbylocus, p_mutsat, 
                       ggarrange(p_mutsat_hist, p_mutsat_dens, p_mutsat_seg, ncol =3),
                       nrow = 3, heights=c(1, 1, 0.75))
  
  ggsave(paste0(save_path, '.png'),
         plot = p_combo, 
         device = 'png',
         width = 25, height = 35, units = 'cm')
  
}


plot_shannonE_all <- function(df, save_path, gene, segment.breaks) {

  
  
  df <- add_gene_segment(df, segment.breaks)
  
  print("Head of df shannonE")
  print(head(df))
  
  p_shannonE <- ggplot(df, aes(x=locus, y=E))+
    geom_hline(yintercept = 1, color='lightgrey')+
    geom_hline(yintercept = 0.8, color='lightgrey')+
    geom_line(linewidth = 0.7)+
    scale_y_continuous(limits=c(0, 1.5))+
    theme_classic()+
    theme(axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          axis.text.y = element_text(size=15),
          axis.text.x = element_text(size=15),
          panel.border = element_rect(fill=NA, colour="black", linewidth=1))+
    ylab('Shannon Equitability')+
    xlab('AA Locus')+
    ggtitle(gene)
  
  
  p_shannonE_hist <- ggplot(df, aes(x=E))+
    geom_histogram(binwidth = 0.05, color = 'black') + 
    #geom_density(color='black')+
    scale_y_continuous(trans='log10', limits=c(NA,1000))+
    scale_x_continuous(limits=c(0,1.5))+
    theme_classic()+
    theme(axis.title.x = element_text(size=13),
          axis.title.y = element_text(size=16),
          axis.text.y = element_text(size=15),
          axis.text.x = element_text(size=15),
          panel.border = element_rect(fill=NA, colour="black", linewidth=1))+
    xlab('Shannon Equitability')+
    ylab('Number of Loci')+
    ggtitle(gene)
  
  p_shannonE_dens <- ggplot(df, aes(x=E))+
    #geom_histogram(binwidth = 1, color = 'black') + 
    geom_density(color='black')+
    #scale_y_continuous(trans='log10', limits=c(NA,1000))+
    scale_x_continuous(limits=c(0,1.5))+
    theme_classic()+
    theme(axis.title.x = element_text(size=13),
          axis.title.y = element_text(size=16),
          axis.text.y = element_text(size=15),
          axis.text.x = element_text(size=15),
          panel.border = element_rect(fill=NA, colour="black", linewidth=1))+
    xlab('Shannon Equitability')+
    ylab('Density')+
    ggtitle(gene)
  
  
  p_shannonE_seg <- ggplot(df, aes(x = GeneSegment, y =E))+
    geom_boxplot(outlier.shape = NA, color = 'black') + 
    geom_jitter(width = 0.3, color = 'black')+
    theme_classic()+
    scale_y_continuous(limits=c(0, 1.5))+
    theme(axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15),
          axis.text.y = element_text(size=15),
          axis.text.x = element_text(size=12),
          panel.border = element_rect(fill=NA, colour="black", linewidth=1))+
    xlab('Gene Segment')+
    ylab('Shannon Equitability')+
    ggtitle(gene)
  
  p_combo <- ggarrange(p_shannonE, 
                       ggarrange(p_shannonE_hist, p_shannonE_dens, p_shannonE_seg, ncol =3),
                       nrow = 3, heights=c(1, 1, 0.75))
  
  ggsave(paste0(save_path, '.png'),
         plot = p_combo, 
         device = 'png',
         width = 25, height = 23, units = 'cm')
  
}



# parse vardict variants and return annotated dataframe
parse_variants <- function(varcounts.df, early.locus, late.locus) {
  
  
  # Proceed with the rest of the function
  varcounts.df <- varcounts.df %>%
    mutate(
      is_wt = aa_mutations == 'wt',
      is_special = grepl("fs|ins|del", aa_mutations),
      
      # Extract wildtype amino acid (if present at the start)
      wildtype = ifelse(
        is_wt, "",
        ifelse(is_special, NA, str_extract(aa_mutations, "^[A-Z\\*]"))
      ),
      
      # Extract position (digits in the middle)
      position = ifelse(
        is_wt, NA, as.numeric(str_extract(aa_mutations, "[0-9]+"))
      ),
      
      # Extract mutation based on case
      mutation = case_when(
        is_wt ~ "",
        is_special ~ str_extract(aa_mutations, "fs|ins|del"),
        TRUE ~ str_extract(aa_mutations, "[A-Z\\*]+$")
      )
    ) %>%
    select(-is_wt, -is_special)
  
  # Extract unique positions and wildtypes
  unique_wildtypes <- varcounts.df %>%
    filter(!is.na(wildtype) & wildtype != "") %>%
    select(position, wildtype) %>%
    distinct()
  
  # Merge to fill missing wildtype values for special cases
  varcounts.df <- varcounts.df %>%
    left_join(unique_wildtypes, by = "position", suffix = c("", ".from_unique")) %>%
    mutate(
      wildtype = coalesce(wildtype, wildtype.from_unique)
    ) %>%
    select(-wildtype.from_unique)
  
  # Updated func_map function to handle NA values
  func_map <- function(wt, mut, pos) {
    if (is.na(wt) || nchar(wt) == 0) {
      return(NA)
    }
    
    if (wt == mut) {
      return("synonymous")
    } else if (mut == "del") {
      return("deletion") 
    } else if (mut == "ins") {
      return("insertion")
    } else if (mut == "fs") {
      return("frameshift")
    } else if (mut == "*" & pos < early.locus) {
      return("early stop")
    } else if (mut == "*" & pos >= early.locus & pos < late.locus) {
      return("medium stop")
    } else if (mut == "*" & pos >= late.locus) {
      return("late stop")
    } else {
      return("missense")
    }
  }
  
  varcounts.df <- varcounts.df %>%
    rowwise() %>%
    mutate(type = func_map(wildtype, mutation, position)) %>%
    ungroup() %>%
    mutate(type = replace_na(type, 'wildtype'))
  
  print(head(varcounts.df))
  
  return(varcounts.df)
}

