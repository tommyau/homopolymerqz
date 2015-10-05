# HomopolymerQZ

## Requirements, as tested on 64-bit CentOS 5.5
* Perl (tested on version 5.8.8)
* Perl module Tie::IxHash (tested on version 1.22)
* SQLite (tested on version 3.7.17)
* R (tested on version 2.14.1)
* R library RSQLite (tested on version 0.11.1)
* R library ggplot2 (tested on version 0.9.0)
* R library sqldf (tested on version 0.4-6.1)
* R library plyr (tested on version 1.7.1)


## Procedures

### Setup reference range from control samples
Extract pyrosequencing signal intensity data and import them into the HomopolymerQZ database
```
AVAPROJECT="/path/to/AVAproject"
ALIGNFOLDER="$AVAPROJECT/Amplicons/Results/Align"
FLOWSTAT="flowstatref.tab"
FLOWSTATDB="homopolymerqz.db"
FLOWSTATTABLE="flowstat"
REFRANGETABLE="refrange"
PROJECTNAME="controlsamples"
BATCHNAME="controlsamples"
PROCESSINGMETHOD="ampliconprocessing"

# init db
echo "create table $FLOWSTATTABLE (sample text, refseq text, pos integer, gapcol integer, direction text, direction_readcnt integer, gaps integer, ns integer, column integer, numgroups integer, genotype  text, genotype_readcnt integer, flowvalue real, flowvalue_readcnt integer, autoanalysis text, batch text, processing text, flowvalue_group real);" | sqlite3 $FLOWSTATDB

# import values
for f in `find $ALIGNFOLDER -name "*.signalDist.txt"`
do
    perl parse_signalDist_v2.pl --projectdef $AVAPROJECT/Amplicons/ProjectDef/ampliconsProject.txt --define autoanalysis=${PROJECTNAME} --define batch=${BATCHNAME} --define processing=${PROCESSINGMETHOD} --input $f | awk -f flowstat_assign_flowvalue_group.awk >> $FLOWSTAT
done

# create database index
echo "CREATE INDEX flowstat_eight_field_index ON $FLOWSTATTABLE(processing, refseq, pos, gapcol, genotype, flowvalue_group, batch, sample);" | sqlite3 $FLOWSTATDB
```

Build reference ranges
```
R
suppressPackageStartupMessages(suppressWarnings(library(RSQLite)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(sqldf)))
suppressPackageStartupMessages(suppressWarnings(library(plyr)))
flowstatdbfile<-"homopolymerqz.db"
sqlite    <- dbDriver("SQLite")
flowstatdb <- dbConnect(sqlite,flowstatdbfile)

oridirmerged<-dbGetQuery(flowstatdb,
"SELECT processing, batch, sample, refseq, pos, gapcol, genotype, flowvalue_group, SUM(flowvalue_readcnt) as flowvalue_group_sum, (SELECT SUM(b.flowvalue_readcnt) FROM flowstat b WHERE a.processing = b.processing AND a.batch = b.batch AND a.sample = b.sample AND a.refseq = b.refseq AND a.pos = b.pos AND a.gapcol = b.gapcol AND a.genotype = b.genotype GROUP BY b.processing, b.batch, b.sample, b.refseq, b.pos, b.gapcol, b.genotype) AS pos_coverage
FROM flowstat a
GROUP BY a.processing, a.batch, a.sample, a.refseq, a.pos, a.gapcol, a.genotype, a.flowvalue_group;")
oridirmerged$flowvalue_group_sum_frac<-oridirmerged$flowvalue_group_sum / oridirmerged$pos_coverage

tasklist<-dbGetQuery(flowstatdb,
"SELECT DISTINCT processing, refseq, pos FROM flowstat;"
)

process_task<-function(x){
    sqlstatement<-paste("SELECT processing, batch, sample, refseq, pos, gapcol, genotype, flowvalue_group, SUM(flowvalue_readcnt) as flowvalue_group_sum, (SELECT SUM(b.flowvalue_readcnt) FROM flowstat b WHERE a.processing = b.processing AND a.batch = b.batch AND a.sample = b.sample AND a.refseq = b.refseq AND a.pos = b.pos AND a.gapcol = b.gapcol AND a.genotype = b.genotype GROUP BY b.processing, b.batch, b.sample, b.refseq, b.pos, b.gapcol, b.genotype) AS pos_coverage
    FROM flowstat a
    WHERE
    a.processing = '", x$processing, "'
    AND a.refseq = '", x$refseq, "'
    AND a.pos = '", x$pos, "'
    GROUP BY a.processing, a.batch, a.sample, a.refseq, a.pos, a.gapcol, a.genotype, a.flowvalue_group;", sep="")
    task_oridirmerged<-dbGetQuery(flowstatdb,sqlstatement)
    task_oridirmerged$flowvalue_group_sum_frac<-task_oridirmerged$flowvalue_group_sum / task_oridirmerged$pos_coverage
    task_oridirmerged_refrange <- ddply(task_oridirmerged, .(processing, refseq, pos, gapcol, genotype, flowvalue_group), summarise, 
		   refN    = length(flowvalue_group_sum_frac),
		   refmean = mean(flowvalue_group_sum_frac),
		   refsd   = sd(flowvalue_group_sum_frac),
		   refse   = sd(flowvalue_group_sum_frac) / sqrt(length(flowvalue_group_sum_frac)),
		   refiqr  = IQR(flowvalue_group_sum_frac),
		   refq1   = quantile(flowvalue_group_sum_frac,1/4),
		   refq3   = quantile(flowvalue_group_sum_frac,3/4)
		   )
    task_oridirmerged_refrange$refq1outlier<-task_oridirmerged_refrange$refq1 - 1.5 * task_oridirmerged_refrange$refiqr
    task_oridirmerged_refrange$refq3outlier<-task_oridirmerged_refrange$refq3 + 1.5 * task_oridirmerged_refrange$refiqr
    return(task_oridirmerged_refrange)
}
refrangetable<-"refrange"
for(i in 1:nrow(tasklist)) {
    out<-process_task(tasklist[i,])
    dbWriteTable(flowstatdb, refrangetable, out, append = T)
}
dbSendQuery(flowstatdb, paste("CREATE INDEX refrange_six_field_index ON ",refrangetable, "(processing, refseq, pos, gapcol, genotype, flowvalue_group)", sep=""))
q("no")
```

### Analyze testing DNA samples

Extract pyrosequencing signal intensity data and import them into the HomopolymerQZ database
```
AVAPROJECT="/path/to/AVAproject"
ALIGNFOLDER="$AVAPROJECT/Amplicons/Results/Align"
FLOWSTAT="flowstat.tab"
FLOWSTATDB="homopolymerqz.db"
FLOWSTATTABLE="flowstat"
REFRANGETABLE="refrange"
PROJECTNAME="testingsamples"
BATCHNAME="testingsamples"
PROCESSINGMETHOD="ampliconprocessing"

# import values
for f in `find $ALIGNFOLDER -name "*.signalDist.txt"`
do
    perl parse_signalDist_v2.pl --projectdef $AVAPROJECT/Amplicons/ProjectDef/ampliconsProject.txt --define autoanalysis=${PROJECTNAME} --define batch=${BATCHNAME} --define processing=${PROCESSINGMETHOD} --input $f | awk -f flowstat_assign_flowvalue_group.awk >> $FLOWSTAT
done

echo -e "DROP TABLE IF EXISTS $FLOWSTATTABLE;\nCREATE TABLE $FLOWSTATTABLE (sample text, refseq text, pos integer, gapcol integer, direction text, direction_readcnt integer, gaps integer, ns integer, column integer, numgroups integer, genotype  text, genotype_readcnt integer, flowvalue real, flowvalue_readcnt integer, autoanalysis text, batch text, processing text, flowvalue_group real);\n.mode tabs\n.import $FLOWSTAT $FLOWSTATTABLE" | sqlite3 $FLOWSTATDB
echo "CREATE INDEX flowstat_eight_field_index ON $FLOWSTATTABLE(processing, refseq, pos, gapcol, genotype, flowvalue_group, batch, sample);" | sqlite3 $FLOWSTATDB
```

Compare a homopolymer locus against reference range, with coordinates according to the reference sequences used in Amplicon Variant Analyzer projects
```
Rscript --no-save --no-restore --slave classify_homopolymer.R $FLOWSTATDB $FLOWSTATTABLE $REFRANGETABLE - <refseq> <start_position> <0> <base>
```
For example
```
Rscript --no-save --no-restore --slave classify_homopolymer.R $FLOWSTATDB $FLOWSTATTABLE $REFRANGETABLE - BRCA1_120329_revcomp 1398 0 G
```