# A selection of scripts written by Jack Shaw
# (and then stolen from Dunhill in prep. for modification)
# 1. pfim network
# 3. helper scripts.

# network, node and motif metric functions ----
pfim <- function(data,
                 col_taxon = "taxon",
                 col_num_size = NULL,
                 cat_combo_list = NULL,
                 cat_trait_types = NULL,
                 num_size_rule = function(res_size, con_size) {
                   ifelse(res_size <= con_size, 1, 0)
                 },
                 certainty_req = "all",
                 return_full_matrix = FALSE,
                 print_dropped_taxa = FALSE,
                 hide_printout = FALSE,
                 ...) {
  
  # Hide "no visible binding for global variable" comment
  . <- Var1 <- Var2 <- pres_sum <- size_consumer <- size_resource <- taxon_consumer <- taxon_resource <- trait <- trait_type <- trait_type_consumer <- trait_type_interaction <- trait_type_pres <- trait_value_pres <- NULL
  
  if (is.null(cat_combo_list)) {
    cat_combo_list <- pfwim::example_data_traitrules
  } else {
    
  }
  
  names(data)[names(data) == col_taxon] <- "taxon"
  col_taxon <- "taxon"
  
  trait_cats <- tolower(unique(c(cat_combo_list[, c("trait_type_resource")], cat_combo_list[, c("trait_type_consumer")])))
  trait_cats_cons <- tolower(unique(c(cat_combo_list[, c("trait_type_consumer")])))
  col_taxon <- tolower(col_taxon)
  
  if (is.null(cat_trait_types)) {
  } else {
    trait_cats <- cat_trait_types
    trait_cats_cons <- cat_trait_types
  }
  
  fd <- as.data.frame(data)
  colnames(fd) <- tolower(colnames(fd))
  
  
  # Run for all trait types
  fw_match_traits <- c()
  for (i in trait_cats_cons) {
    col_combo_list_mini <- cat_combo_list %>% dplyr::filter(trait_type_consumer == i)
    
    res_taxa <- fd %>% dplyr::select(c(col_taxon, unique(col_combo_list_mini[, c("trait_type_resource")])))
    con_taxa <- fd %>%
      dplyr::select(c(col_taxon, i)) %>%
      filter_all(any_vars(!is.na(.)))
    con_taxa <- con_taxa[Reduce(`&`, lapply(con_taxa, function(x) !(is.na(x) | x == ""))), ]
    
    res_taxa_longer <- res_taxa %>%
      dplyr::mutate_all(., as.character) %>%
      tidyr::pivot_longer(-dplyr::all_of(col_taxon), names_to = "trait_type", values_to = "trait") %>%
      dplyr::filter(trait_type %in% col_combo_list_mini[, c("trait_type_resource")])
    
    con_taxa_longer <- con_taxa %>%
      dplyr::mutate_all(., as.character) %>%
      tidyr::pivot_longer(-dplyr::all_of(col_taxon), names_to = "trait_type", values_to = "trait") %>%
      dplyr::filter(trait != "primary")
    
    crossing_taxa <- tidyr::crossing(as.data.frame(res_taxa_longer) %>% 
                                       setNames(paste0(names(.), "_resource")), con_taxa_longer %>% 
                                       setNames(paste0(names(.), "_consumer"))) %>%
      dplyr::left_join(col_combo_list_mini[, c("trait_type_resource", "trait_type_consumer", "trait_consumer")] %>% 
                         dplyr::mutate(trait_type_pres = 1), by = c("trait_type_resource", "trait_type_consumer", "trait_consumer"))
    crossing_taxa <- crossing_taxa %>%
      dplyr::filter(trait_type_pres == 1) %>%
      dplyr::distinct() %>%
      dplyr::left_join(col_combo_list_mini %>% 
                         dplyr::mutate(trait_value_pres = 1), by = c("trait_type_resource", "trait_type_consumer", "trait_consumer", "trait_resource"))
    crossing_taxa <- crossing_taxa %>%
      dplyr::filter(trait_value_pres == 1) %>%
      dplyr::mutate(trait_type_interaction = i)
    crossing_taxa <- crossing_taxa[, c("taxon_resource", "taxon_consumer", "trait_type_interaction", "trait_value_pres")]
    
    fw_match_traits <- rbind(fw_match_traits, crossing_taxa)
    if (hide_printout == FALSE) print(i)
  }
  
  # If feasibility of interaction also defined by a numerical size rule, then utilize
  if (is.null(col_num_size)) {
    
  } else {
    fw_size_pastes <- unique(paste(fd[, c(col_taxon)], fd[, c(col_num_size)], sep = "SEP"))
    fw_size_expand <- as.data.frame(expand.grid(fw_size_pastes, fw_size_pastes)) %>%
      tidyr::separate(Var1, c("taxon_resource", "size_resource"), sep = "SEP", convert = TRUE) %>%
      tidyr::separate(Var2, c("taxon_consumer", "size_consumer"), sep = "SEP", convert = TRUE) %>%
      dplyr::mutate(trait_value_pres = num_size_rule(size_resource, size_consumer), trait_type_interaction = col_num_size) %>%
      dplyr::filter(trait_value_pres == 1) %>%
      dplyr::select(taxon_resource, taxon_consumer, trait_value_pres, trait_type_interaction)
    fw_match_traits <- bind_rows(fw_match_traits, fw_size_expand)
    if (hide_printout == FALSE) print(col_num_size)
  }
  
  
  fw_match_traits2 <- fw_match_traits %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(names_from = trait_type_interaction, values_from = trait_value_pres)
  fw_match_traits2$pres_sum <- fw_match_traits2 %>%
    dplyr::select(all_of(unique(fw_match_traits$trait_type_interaction))) %>%
    rowSums(na.rm = TRUE)
  
  # How many trait combinations need to be matched in order for an interaction to be viable
  if (certainty_req == "all") {
    certainty_val <- length(unique(fw_match_traits$trait_type_interaction))
  } else {
    certainty_val <- certainty_req
    certainty_val <- ifelse(certainty_val > length(unique(fw_match_traits$trait_type_interaction)), length(unique(fw_match_traits$trait_type_interaction)), certainty_val)
  }
  
  
  # Return full matrix with trait feasibility indicated or not
  if (return_full_matrix == TRUE) {
    fw_match_traits3 <- fw_match_traits2 %>% dplyr::distinct()
    # print(paste(length(setdiff(unique(fd$taxon), unique(c(fw_match_traits3$taxon_resource, fw_match_traits3$taxon_consumer)))), "taxa dropped from web"))
    
    
    
    return(as.matrix(fw_match_traits3))
  } else {
    fw_match_traits3 <- fw_match_traits2 %>%
      dplyr::filter(pres_sum >= certainty_val) %>%
      dplyr::distinct()
    if (hide_printout == FALSE) print(paste(length(setdiff(unique(fd$taxon), unique(c(fw_match_traits3$taxon_resource, fw_match_traits3$taxon_consumer)))), "taxa dropped from web"))
    
    if (print_dropped_taxa == TRUE) {
      if (hide_printout == FALSE) print(setdiff(unique(fd$taxon), unique(c(fw_match_traits3$taxon_resource, fw_match_traits3$taxon_consumer))))
    } else {
    }
    
    edgelist <- as.matrix(fw_match_traits3[, c("taxon_resource", "taxon_consumer")])
    
    # Create an adjacency matrix called benguela.AM, containing only zeros
    spp_list <- edgelist %>%
      as.tibble() %>%
      pivot_longer(
        everything(),
        names_to = c("name"),
        values_to = c("species")
      ) %>%
      distinct(species)
    
    spp_list <- as.vector(spp_list$species)
    
    
    int_mat <- matrix(0, length(spp_list), length(spp_list),
                      dimnames = list(spp_list, spp_list))
    
    for (i in 1:nrow(edgelist)) {
      resource = edgelist[i,1]
      consumer = edgelist[i,2]
      
      int_mat[consumer,resource] = 1
    }
    
    return(int_mat)
  }
}


sample_pdf <- function(M = 100,
                       y = 2.5,
                       func = function(r, M, y) exp(-r / (exp((y - 1) * (log(M) / (y))))),
                       n_samp = 100) {
  row <- c()
  for (i in 1:M) {
    row <- c(row, func(i, M, y))
  }
  row2 <- as.data.frame(row)
  # ggplot(row2)+geom_point(aes(x=1:M,y=row))
  
  ar <- sum(row)
  row <- row / ar
  
  return(sample(1:M, n_samp, replace = T, prob = row))
}


powerlaw_prey <- function(el,
                          n_samp = 50,
                          func = function(r, M, y) exp(-r / (exp((y - 1) * (log(M) / (y)))))) {
  con_node_node_name_inferred <- NULL
  
  edgelist <- as.data.frame(el)
  colnames(edgelist) <- c("res_node_node_name_inferred", "con_node_node_name_inferred")
  
  web_list <- lapply(1:n_samp, matrix, data = NA, nrow = 0, ncol = 2)
  
  for (i in unique(edgelist$con_node_node_name_inferred)) {
    min <- edgelist %>% dplyr::filter(con_node_node_name_inferred == i)
    min <- as.data.frame(min)
    rich <- as.numeric(as.character(length(unique(min[, c("res_node_node_name_inferred")]))))
    
    t <- sample_pdf(M = rich, n_samp = n_samp)
    
    names(t) <- rep(i, length(t))
    
    for (j in 1:length(t)) {
      
      # Currently taking upper bound
      min2 <- min %>% sample_n(min(rich, t[[j]]))
      
      
      
      web_list[[j]] <- rbind(web_list[[j]], min2)
    }
  }
  
  web_list <- lapply(web_list, as.matrix)
  web_list <- lapply(web_list, na.omit)
  
  
  
  return(web_list)
}

#### UTILS----

mat2list <- function(mat) {
  mat <- as.matrix(mat)
  list <- data.frame(X1 = colnames(mat)[col(mat)], X2 = rownames(mat)[row(mat)], dist = c(mat))
  colnames(list) <- c("col", "row")
  return(list)
}


clean_date <- function(drop_punctuation = TRUE) {
  t <- format(Sys.time(), format = "%Y-%m-%dT%H:%M:%S%Z")
  
  t <- ifelse(drop_punctuation == FALSE, t, gsub("[^0-9A-Za-z///' ]", "", t, ignore.case = TRUE))
  
  return(t)
}


jlookup <- function(old_values,
                    replacement_df,
                    matching_col,
                    new_values) {
  temp_df <- as.data.frame(old_values) %>%
    dplyr::left_join(replacement_df[, c(matching_col, new_values)], by = c("old_values" = matching_col))
  
  return(temp_df[, c(new_values)])
}
