yeastomics="https://raw.githubusercontent.com/benjamin-elusers/yeastomics/main/"
source(file.path(yeastomics,"src/utils.r"))
source(file.path(yeastomics,"src/function_annotation.r"))
source(file.path(yeastomics,"src/function_sequence.r"))
source(file.path(yeastomics,"src/function_phylogenetic.r"))
source(file.path(yeastomics,"src/function_analysis.r"))
source(file.path(yeastomics,"src/function_datalocal.r"))
source(file.path(yeastomics,"src/function_datapub.r"))
library(tidyverse)
library(hablar)

##### 1. UNIPROT-TO-ORF MAPPING #####
sc_uni = readRDS("data/uniprot-features.rds") %>%
  dplyr::select("UNIPROTKB","SGD","EXISTENCE","SCORE","FAMILIES","PNAME")
sc_sgd = load.sgd.features() %>% dplyr::select("sgdid","type","qual",ORF="name","gname","chr","strand")
uni2sgd = left_join(sc_uni,sc_sgd,by=c('SGD'='sgdid'))

##### 2. PAXDB ABUNDANCE RANKED #####
sc_pax = get.paxdb(4932,'integrated') %>%
  group_by(taxid) %>%
  mutate(rk = rank(desc(ppm_int)),
         rk.pc = percent_rank(desc(ppm_int)),
         pc= 100 * ppm_int / sum_(ppm_int)) %>%
  arrange(rk.pc) %>% mutate(pc.csum = cumsum(pc))

##### 3. PAXDB MAPPED TO UNIPROT-TO-ORF MAPPING #####
sc_pax2orf = left_join(sc_pax,uni2sgd, by=c('protid'='ORF')) %>%
  relocate(UNIPROTKB,EXISTENCE,SCORE,FAMILIES,
           ORF=protid,gname,SGD,type,qual,chr,strand)

##### 4. PDB TO ORF MAPPING VIA UNIPROT #####
#sc_complex = load.3dcomplex.yeast() # residue level data
sc_pdb = get.mapping.3dcomplex.yeast()
sc_pdb_nr = sc_pdb %>%
            arrange(seqid,code,chain_name,desc(ident),desc(overlap),resol,desc(length_atom)) %>%
            group_by(seqid) %>%
            mutate( pdb_rk=row_number() ) %>% #, length_atom, length_full) ) %>%
            dplyr::filter(pdb_rk==1)# REQUIRES ACCESS TO FILESERVER 'mata.weizmann.ac.il'

##### 5. RANKED PROTEIN ABUNDANCE MAPPED TO PDB #####
options(dplyr.width=Inf)
final = left_join(sc_pax2orf,sc_pdb_nr, by=c('UNIPROTKB'='seqid'))

top50 = final %>%
  dplyr::filter(rk < 51) %>%
  dplyr::select(UNIPROTKB,ORF,gname, ppm_int,rk, rk.pc, pc, pc.csum,
                code, chain_name, ident, overlap, resol)
top50 %>%
  dplyr::select(rk, gname, UNIPROTKB, ppm_int, pc.csum, code, pid=ident, ol=overlap, R=resol) %>%
  distinct() %>%
  print(n=30)

homology = tribble(
  ~name,        ~pdb,   ~chain,       ~oligomer,    ~pid,  ~is_homology,
  'FBA1',   '6lnk.1',     'A',    'homo-2-mer',   73.43,          TRUE,
  'ENO2',   '2one.1',     'A',    'homo-2-mer',   95.41,          TRUE,
 'CDC19',     '1a3x',  'ABCD',    'homo-4-mer',   100.0,         FALSE,
  'SSA2',   '4fl9.1',     'A',       'monomer',   80.18,          TRUE,
  'PMA1',   '6lly.1',     'A',       'monomer',   25.50,          TRUE,
  'ILV5',   '6vo2.1',     'A',    'homo-2-mer',   34.38,          TRUE,
  'ALD6',   '5fhz.1',     'A',    'homo-4-mer',   48.53,          TRUE,
  'SSA1',   '4fl9.1',     'A',       'monomer',   80.36,          TRUE,
'RPL21B',   '5h4p.1',     'U',       'monomer',   98.75,          TRUE,
 'RPP2B',  '4v5z.23',     'A',       'monomer',   35.09,         FALSE,
  'MET6',   '3ppc.1',     'A',       'monomer',   75.52,          TRUE,
 'RPL7A',     '6c0f',     'I', 'hetero-40-mer',   100.0,         FALSE
)

top20 = top50 %>% left_join(homology, by=c('gname'='name')) %>%
  mutate( pdb_code = ifelse(is.na(code), no=code, yes=pdb),
          pdb_chain = ifelse(is.na(chain_name), no=chain_name, yes=chain),
          ident     = ifelse(is.na(ident), no=ident, pid)
        )  %>%
  dplyr::select(-c(code, chain_name, ident)) %>% dplyr::filter(rk<21)

sprintf("fetch %s",paste(str_sub(top20$pdb_code,1,4),collapse=" "))



