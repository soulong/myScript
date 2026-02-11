
library(readr)
library(webchem)


get_pubchem <- function(query, from="smiles") {
  
  query <- as.character(na.omit(lib$smiles))
  
  res <- get_cid(
    query=query,
    from=from,
    domain="compound") %>% 
    transmute(Query=query, CID=cid)
  
  prop <- pc_prop(
    cid=res$CID, 
    properties=c("IUPACName","CanonicalSMILES","InChIKey","MolecularWeight")) %>% 
    as_tibble()
  
  if(is.na(prop)) {
    print("retrieve props failed, return CID only")
    return(res)
    
  } else {
    res <- left_join(res, 
                     mutate(prop, CID=as.character(CID)), 
                     by=join_by(CID))
    return(res)
  }
  
}



if(F) {
  lib <- "/media/hao/Data1/compound_screen/cmpd_library/shulab_cmpd_library_merged.csv" %>% 
    read_csv() %>% glimpse()
  
  data <- lib$smiles %>%
    # .[1:10] %>% 
    get_pubchem()
  write_csv(data, "/media/hao/Data1/compound_screen/cmpd_library/shulab_cmpd_library_merged_pubchem.csv")
}



