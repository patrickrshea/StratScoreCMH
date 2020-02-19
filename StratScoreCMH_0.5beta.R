#!/usr/bin/env Rscript
### load libraries -----------------------------------------------------------------

library("optparse")
library("plyr")


### parse options from command line ------------------------------------------
option_list = list(
  make_option(c("-m", "--matrix"), type="character", default=NULL, 
              help="collapsing matrix file", metavar="character"),
  make_option(c("-p", "--pheno"), type="character", default=NULL, 
              help="phenotype file", metavar="character"),
  make_option(c("-c", "--cluster"), type="character", default=NULL, 
              help="stratification score file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="results.csv", 
              help="output file name [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$matrix)){
  print_help(opt_parser)
  stop("At least three arguments must be supplied (collapsing matrix, phenotype, and cluster score files)./n", call.=FALSE)
}

### functions -----------------------------------------------------------------

doCMH <- function(temp) {
	contingency_table<-xtabs(~ Gene + Group + Strata, data=temp)
	if (is.null(contingency_table)) {
		model.results$p.value <- "NA"
	} else if (sum(temp$Gene)==0)  {   ##CATCH ANY INVARIANT GENES AND ASSIGN A P-VALUE OF 1##
		model.results$p.value <- "1"
	} else {	
		gene.model<-mantelhaen.test(contingency_table, exact = TRUE)
		model.results$p.value <- gene.model$p.value
		}	
return(model.results)
}


### main ---------------------------------------------------------------------


###READ COLLAPSING MATRIX IN FROM FILE
df<-read.table(opt$matrix, header=T, row.names=1)
df<- as.data.frame(t(df))


###READ AND EXTRACT PHENOS
df1<-read.table(opt$pheno, header=T)
phenos<-subset(df1, select=c("Sample", "Group"))


##MERGE PHENOS
df3<-merge(df, phenos, by.x="row.names", by.y="Sample")
rownames(df3)<-df3$Row.names
df3<-df3[-1]


### READ STRATIFICATION SCORE CLUSTERS FROM PROPENSITY MATCHING
#
#   input format: sampleID clusterID
############################################
strata<-read.table(opt$cluster, row.names=1)
colnames(strata)<-c("strataID")

cmh_array<-merge(df3, strata, by="row.names")
rownames(cmh_array)<-cmh_array$Row.names
collapsed<-cmh_array[-1]

###REMOVE ANY STRATA WITH N=1 (caused by lack of overlap between clusters and genotype matrix )
strata_counts<-count(collapsed, vars="strataID")
bad_strata<-strata_counts[which(strata_counts$freq<2),1]
collapsed<-collapsed[! collapsed$strataID %in% bad_strata, ]

##INITIALIZE STORAGE VARIABLES
results=NULL
model.results=NULL
temp=NULL

genes <- colnames(collapsed)
genes<-genes[-((length(genes) - 2):length(genes))]

##LOOP THROUGH GENES AND PERFORM CMH TEST

for (gene in genes) {

	genotypes <- collapsed[gene]
	genotypes <- ifelse(genotypes > 0, 1, 0) # convert from ATAV format
	colnames(genotypes) <- 'gene'

	###EXTRACT GENOTYPES, CASE_CONTROL STATUS, AND STRATA INTO SEPARATE DATAFRAME
	temp<-data.frame(genotypes, collapsed$Group, collapsed$strataID)
	colnames(temp)<-c("Gene", "Group", "Strata") 


	##CREATE 3D CONTINGENCY TABLE AND CALCULATE CMH TEST
	model_results<-doCMH(temp)

       
	   results = rbind(results, data.frame(gene, model_results))
	}

##OUTPUT RESULTS
write.csv(results, opt$out)
