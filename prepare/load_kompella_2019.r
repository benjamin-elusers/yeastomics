library(tabulizer)
library(tidyverse)
kompella_et_al = tabulizer::extract_tables(file="~/Downloads/Data_Sheet_1_Definition of the Minimal Contents for the Molecular Simulation of the Yeast Cytoplasm.pdf",
                                           pages = 4:7, output='data.frame', columns=list(5,5,5,5))

SI.1 = kompella_et_al[[1]] %>%
       dplyr::rename(rk='X',orf='X.1',gene='X.2',desc=4,pdb=5) %>%
       dplyr::filter(!is.na(rk))
SI.2 = kompella_et_al[[2]] %>%
        dplyr::rename(gene='X.1',desc=3,pdb='X.2') %>%
        dplyr::filter(X != '') %>%
        mutate(rk = row_number()+nrow(SI.1), orf=str_extract(X,SGD.nomenclature())) %>%
        relocate(rk,orf) %>%
        dplyr::select(-X)
SI.3 = kompella_et_al[[3]] %>%
       dplyr::rename(gene='X.1',desc=3,pdb='X.2') %>%
       dplyr::filter(X!='') %>%
       mutate(rk = row_number()+nrow(SI.1)+nrow(SI.2), orf=str_extract(X,SGD.nomenclature())) %>%
       relocate(rk,orf) %>%
       dplyr::select(-X)
SI.4 = kompella_et_al[[4]] %>%
        dplyr::rename(rk='X',orf='X.1',gene='X.2',desc=4,pdb='X.3') %>%
        dplyr::filter(!is.na(rk))

SI = bind_rows(SI.1,SI.2,SI.3,SI.4) %>%
      dplyr::filter(!is.na(rk)) %>%
      mutate(orf = str_replace_all(orf," ",""),
             gene = str_replace_all(gene," ",""),
             pdb = str_extract(str_replace_all(pdb," ",""), "\\d\\w{3}_?\\w?")) %>%
      dplyr::select(-desc)

