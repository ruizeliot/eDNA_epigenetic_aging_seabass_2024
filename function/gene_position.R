
gene_position = function(position_gb, genes, pattern_list, type_list, col_name = "V1"){ 
  
  get_gene_position = function(data, search_pattern, type, name, col_name){
    
    # Create an empty list to store the orientation of the DNA sequence (5' to 3' = "normal" // 3' to 5' = "complement")
    orientation = list()
    
    # For each sequence splitted in lapply, search the row number in which there is the required pattern
    position = lapply(data, function(x) grep(search_pattern, x[[col_name]])) 
    
    # For each sequence splitted in lapply, search the row number in which there is the required type of gene
    index_gene = lapply(data, function(x) grep(type, x[[col_name]]))
    
    # For each sequence, takes the closest position of the row containing the pattern "type" (index_gene) 
    # to the position of the row containing the pattern giving the name of the gene (position)
    for(i in 1:length(index_gene)){ index_gene[[i]] = tail(index_gene[[i]][index_gene[[i]] < tail(position[[i]], 1)], 1)}
    
    # For each sequence, extract the row containing the position based on the index of the row found above
    for(i in 1:length(data)){ position[[i]] = data[[i]][index_gene[[i]],] }
    
    # For each sequence, searches if there is the mention "complement" in the row containing the positions
    for(j in 1:length(data)){ 
      
      if(length(grep("complement", position[[j]])) != 0) orientation[j] = "complement"
      
      else orientation[j] = "normal" 
      
    }
    
    # Once the orientation infos has been extracted, removes all mention about the orientation
    position = lapply(position, function(x) gsub("join\\(|,|complement\\(|\\)|>|<", "", x))
    
    # Convert the positions stored as list elements in a vector. Unknown positions (no element) are assigned NA here
    position = unlist(lapply(position, function(x) ifelse(length(x) == 0, NA, x)))
    
    # Extract the starting position by looking at the last word (first = type), and then the numbers before ".."
    start_position = as.numeric(stringr::word(stringr::word(position, -1), 1, sep = fixed("..")))
    
    # Extract the starting position by looking at the last word (first = type), and then the numbers after ".."
    end_position = as.numeric(stringr::word(stringr::word(position, -1), -1, sep = fixed("..")))
    
    # Combine all informations in a dataframe
    output = tibble::tibble(data.frame(start_position, end_position, unlist(orientation)))
    
    # Name the column of the dataframe by pasting the type of the infos with the name specified in the argument of the function
    colnames(output) = paste(c("START", "END", "ORIENTATION"), name, sep = "_")
    
    return(output)
    
  }
  
  # Split the dataframe after each occurence of the word "LOCUS" which is the first category (left column on NCBI) in each descriptions  
  position_splited = split(position_gb, findInterval(1:nrow(position_gb), which(grepl("LOCUS", position_gb[[col_name]])) + 1))[-1]
  
  # Launch the function get_gene_position for all possible names of the gene or mispellings 
  position_list = lapply(1:length(genes), function(x) 
    suppressWarnings(get_gene_position(position_splited, pattern_list[[genes[x]]], type_list[[genes[x]]], genes[x], col_name = col_name)))
  
  # To get the accession number which is always the last element of the row named "VERSION"
  accession = unlist(lapply(position_splited, function(x) x[which(grepl("VERSION", x[[col_name]])), ]))
  
  # To remove the version number from the accession number
  accession = tibble::tibble(data.frame(ACCESSION = stringr::word(stringr::word(accession, -2, sep = fixed(".")), -1)))
  
  # To combine all 4 dataframes per genes (possible because different column names each times)
  final_df = tibble::tibble(cbind(accession, do.call(cbind, position_list)))
  
  # To return the results in the form of a list separating complete mitogenomes in which all genes position were found, and other sequences missing those infos or not being complete mitogenome
  final_list = list(COMPLETE_POSITIONS = final_df[rowSums(is.na(final_df)) == 0, ],
                    MISSING_POSITIONS = final_df[rowSums(is.na(final_df)) %in% c(1:7), ])
  if(nrow(final_list[[1]]) == 0) return(final_list[[2]])
  else if(nrow(final_list[[2]]) == 0) return(final_list[[1]])
  else return(final_list)
  
}
