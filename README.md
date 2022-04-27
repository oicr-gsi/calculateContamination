# calculateContamination

QC workflow to determine contamination metrics on tumor bam files.

## Overview

![calculateContamination workflow diagram](./calculateContamination.svg?raw=true "calculateContamination workflow diagram")

## Dependencies

* [gatk 4.2.0.0](https://gatk.broadinstitute.org)
* [hg38-gatk-gnomad 2.0](https://gnomad.broadinstitute.org/)


## Usage

### Cromwell
```
java -jar cromwell.jar run calculateContamination.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`inputType`|String|Either 'bam' or 'fastq'
`bwaMem.runBwaMem_bwaRef`|String|The reference genome to align the sample with by BWA
`bwaMem.runBwaMem_modules`|String|Required environment modules


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`inputGroups`|Array[InputGroup]?|None|Array of fastq structs containing reads and readgroups
`bamFiles`|Array[BamInputs]?|None|Array of bam structs containing bam files and indices


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`bwaMem.adapterTrimmingLog_timeout`|Int|48|Hours before task timeout
`bwaMem.adapterTrimmingLog_jobMemory`|Int|12|Memory allocated indexing job
`bwaMem.indexBam_timeout`|Int|48|Hours before task timeout
`bwaMem.indexBam_modules`|String|"samtools/1.9"|Modules for running indexing job
`bwaMem.indexBam_jobMemory`|Int|12|Memory allocated indexing job
`bwaMem.bamMerge_timeout`|Int|72|Hours before task timeout
`bwaMem.bamMerge_modules`|String|"samtools/1.9"|Required environment modules
`bwaMem.bamMerge_jobMemory`|Int|32|Memory allocated indexing job
`bwaMem.runBwaMem_timeout`|Int|96|Hours before task timeout
`bwaMem.runBwaMem_jobMemory`|Int|32|Memory allocated for this job
`bwaMem.runBwaMem_threads`|Int|8|Requested CPU threads
`bwaMem.runBwaMem_addParam`|String?|None|Additional BWA parameters
`bwaMem.adapterTrimming_timeout`|Int|48|Hours before task timeout
`bwaMem.adapterTrimming_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.adapterTrimming_addParam`|String?|None|Additional cutadapt parameters
`bwaMem.adapterTrimming_modules`|String|"cutadapt/1.8.3"|Required environment modules
`bwaMem.slicerR2_timeout`|Int|48|Hours before task timeout
`bwaMem.slicerR2_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.slicerR2_modules`|String|"slicer/0.3.0"|Required environment modules
`bwaMem.slicerR1_timeout`|Int|48|Hours before task timeout
`bwaMem.slicerR1_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.slicerR1_modules`|String|"slicer/0.3.0"|Required environment modules
`bwaMem.countChunkSize_timeout`|Int|48|Hours before task timeout
`bwaMem.countChunkSize_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.outputFileNamePrefix`|String|"output"|Prefix for output file
`bwaMem.numChunk`|Int|1|number of chunks to split fastq file [1, no splitting]
`bwaMem.doTrim`|Boolean|false|if true, adapters will be trimmed before alignment
`bwaMem.trimMinLength`|Int|1|minimum length of reads to keep [1]
`bwaMem.trimMinQuality`|Int|0|minimum quality of read ends to keep [0]
`bwaMem.adapter1`|String|"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"|adapter sequence to trim from read 1 [AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC]
`bwaMem.adapter2`|String|"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"|adapter sequence to trim from read 2 [AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT]
`tumorOnlyMetrics.modules`|String|"gatk/4.2.0.0 hg38-gatk-gnomad/2.0"|Required environment modules
`tumorOnlyMetrics.refVCF`|String|"$HG38_GATK_GNOMAD_ROOT/small_exac_common_3.hg38.vcf.gz"|Path the reference VCF required by GATK
`tumorOnlyMetrics.memory`|Int|24|Memory allocated for this job
`tumorOnlyMetrics.timeout`|Int|12|Time in hours before task timeout
`getMetrics.refVCF`|String|"$HG38_GATK_GNOMAD_ROOT/small_exac_common_3.hg38.vcf.gz"|Path the reference VCF required by GATK
`getMetrics.modules`|String|"gatk/4.2.0.0 hg38-gatk-gnomad/2.0"|Required environment modules
`getMetrics.memory`|Int|24|Memory allocated for this job
`getMetrics.timeout`|Int|12|Time in hours before task timeout


### Outputs

Output | Type | Description
---|---|---
`contaminationMetrics`|File|Metrics about contamination for inputs bams/fastqs


## Commands
 This section lists command(s) run by calculateContamination workflow
 
 * Running calculateContamination
 
 A WDL workflow that runs GATK4 CalculateContamination from FastQ or bam files.
 
 ### Create pileups for tumor/normal pair and compare.
 
 ```
 module load ~{modules}
 
 mv ~{tumorBamFile} ./tumorBamFile.bam
 mv ~{tumorBaiFile} ./tumorBamFile.bam.bai
 
 mv ~{normalBamFile} ./normalBamFile.bam
 mv ~{normalBaiFile} ./normalBamFile.bam.bai
 
 gatk GetPileupSummaries \
 -I tumorBamFile.bam \
 -V ~{refVCF} \
 -L ~{refVCF} \
 -O tumor.summaries.table
 
 gatk GetPileupSummaries \
 -I normalBamFile.bam \
 -V ~{refVCF} \
 -L ~{refVCF} \
 -O normal.summaries.table
 
 gatk CalculateContamination \
 -I tumor.summaries.table \
 -matched normal.summaries.table \
 -O contamination.table
 
 ```
 ### Create pileup for tumor bam and get contamination metrics.
 
 ```
 module load ~{modules}
 
 mv ~{tumorBamFile} ./tumorBamFile.bam
 mv ~{tumorBaiFile} ./tumorBamFile.bam.bai
 
 gatk GetPileupSummaries \
 -I tumorBamFile.bam \
 -V ~{refVCF} \
 -L ~{refVCF} \
 -O tumor.summaries.table
 
 gatk CalculateContamination \
 -I tumor.summaries.table \
 -O contamination.table
 
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
