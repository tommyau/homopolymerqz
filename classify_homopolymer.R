# classify_homopolymer.R
#   Compare pyrosequencing signal intensity values against reference ranges
#
# Chun Hang AU (chau@hksh.com)
#   Hong Kong Sanatorium and Hospital
#
suppressPackageStartupMessages(suppressWarnings(library(RSQLite)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(sqldf)))
suppressPackageStartupMessages(suppressWarnings(library(plyr)))

args <- commandArgs(trailingOnly = TRUE)
# use "-" as argument to print to STDOUT

flowstatdbfile<-args[1]
flowstattable<-args[2]
refrangetable<-args[3]
fileout<-args[4]
refseq<-args[5]
pos<-args[6]
gapcol<-args[7]
genotype<-args[8]

sqlite    <- dbDriver("SQLite")
flowstatdb <- dbConnect(sqlite,flowstatdbfile)

# support print to STDOUT
if (fileout == "-") fileout<-stdout()

out<-dbGetQuery(flowstatdb, paste("SELECT a.autoanalysis, a.processing, a.batch, a.sample, a.refseq, a.pos, a.gapcol, a.genotype, a.flowvalue_group, SUM(a.flowvalue_readcnt) as flowvalue_group_sum, (SELECT SUM(b.flowvalue_readcnt) FROM flowstat b WHERE a.processing = b.processing AND a.batch = b.batch AND a.sample = b.sample AND a.refseq = b.refseq AND a.pos = b.pos AND a.gapcol = b.gapcol AND a.genotype = b.genotype GROUP BY b.processing, b.batch, b.sample, b.refseq, b.pos, b.gapcol, b.genotype) AS pos_coverage, r.processing, r.refseq, r.pos, r.gapcol, r.genotype, r.flowvalue_group, r.refN, r.refmean, r.refsd, r.refse, r.refiqr, r.refq1, r.refq3, r.refq1outlier, r.refq3outlier
FROM ", flowstattable," a, ", refrangetable," r
WHERE a.processing = r.processing AND a.refseq = r.refseq AND a.pos = r.pos AND a.gapcol = r.gapcol AND a.genotype = r.genotype AND a.flowvalue_group = r.flowvalue_group
 AND a.processing='ampliconprocessing' AND a.refseq = '", refseq,"' AND a.pos = ",pos," AND a.gapcol = ", gapcol," AND a.genotype = '", genotype,"'
GROUP BY a.processing, a.batch, a.sample, a.refseq, a.pos, a.gapcol, a.genotype, a.flowvalue_group;", sep=""))
if (nrow(out) >= 1) {
 out$flowvalue_group_sum_frac<-out$flowvalue_group_sum / out$pos_coverage
 out$classq3outlier <- out$flowvalue_group_sum_frac >= out$refq3outlier
 out$zscore <- ( out$flowvalue_group_sum_frac - out$refmean ) / out$refsd
 
 # since append=TRUE and col.names=TRUE, this leads to the following warning message
 #Warning message:
 #In write.table(subset(out, classq3outlier == TRUE & flowvalue_group_sum_frac >=  :
 #  appending column names to file
 # we are going to suppress such wanring
 suppressWarnings(write.table(subset(out, classq3outlier==TRUE & flowvalue_group_sum_frac >= 0.05), file=fileout, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t", append=TRUE))
}