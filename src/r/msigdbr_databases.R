# GO Biological Processes (GO:BP)
db_go <- msigdbr(species = "Mus musculus") %>%
  filter(gs_subcollection == "GO:BP") %>%   # NOTE: some msigdbr versions use gs_subcat; others use gs_subcollection
  transmute(term = gs_name, gene = gene_symbol) %>%
  distinct()

# KEGG
db_kegg <- msigdbr(species = "Mus musculus") %>%
  filter(gs_subcollection == "CP:KEGG") %>%
  transmute(term = gs_name, gene = gene_symbol) %>%
  distinct()

# KEGG Legacy
db_kegg <- msigdbr(species = "Mus musculus", category = "C2", subcollection = "CP:KEGG_LEGACY") %>%
  transmute(term = gs_name, gene = gene_symbol) %>%
  distinct()
# Hallmark (H)
db_hallmark <- msigdbr(species = "Mus musculus", category = "H") %>%
  transmute(term = gs_name, gene = gene_symbol) %>%
  distinct()
# Reactome (CP:REACTOME)
db_reactome <- msigdbr(species = "Mus musculus") %>%
  filter(gs_subcollection == "CP:REACTOME") %>%   # NOTE: some msigdbr versions use gs_subcat; others use gs_subcollection
  transmute(term = gs_name, gene = gene_symbol) %>%
  distinct()
