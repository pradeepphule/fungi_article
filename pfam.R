ldomai
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# Enrichment analysis for Pfam domains
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

## define the input data
data <- as.character(domains$Pfam)

### load Pfam domain informtion (as 'InfoDataFrame' object)
Pfam <- dcRDataLoader('Pfam')
Pfam
## 1) GOBP enrichment analysis, producing an object of S4 class 'Eoutput'
### By default, using all annotatable domains as the background
eoutput <- dcEnrichment(data, domain="Pfam", ontology="GOBP")
eoutput
#### look at Pfam domains annotated by the most signficant term
tmp <- as.character(view(eoutput, top_num=1, sortBy="pvalue", details=T)$members)
tmp <- unlist(strsplit(tmp,","))
Data(Pfam)[match(tmp,rowNames(Pfam)),]
## 2) GOMF enrichment analysis, producing an object of S4 class 'Eoutput'
eoutput <- dcEnrichment(data, domain="Pfam", ontology="GOMF")