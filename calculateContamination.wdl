version 1.0

import "imports/pull_bwaMem.wdl" as bwaMem

struct InputGroup {
    File fastq1
    File fastq2
    String readGroups
}

struct BamInputs {
    File bamFile
    File baiFile
}

workflow calculateContamination {
    input {
        Array[InputGroup]? inputGroups
        Array[BamInputs]? bamFiles
        String inputType
        String refVCF
        String modules 
    }

    parameter_meta {
        inputGroups: "Array of fastq structs containing reads and readgroups"
        bamFiles: "Array of bam structs containing bam files and indices"
        inputType: "Either 'bam' or 'fastq'"
        refVCF: "Path the reference VCF required by GATK"
        modules: "Required environment modules"
    }
 
    meta {
        author: "Murto Hilali"
        email: "mhilali@oicr.on.ca"
        description: "QC workflow to determine contamination metrics on tumor bam files."
        dependencies: [
            {
                name: "gatk/4.2.0.0",
                url: "https://gatk.broadinstitute.org"
            },
            {
                name: "hg38-gatk-gnomad/2.0",
                url: "https://gnomad.broadinstitute.org/"
            }
        ]
        output_meta: {
            contaminationMetrics: "Metrics about contamination for inputs bams/fastqs"

        }
    }

# =======================================================
#   Accept fastqs and align them into bam files.
#   Bam and index file(s) collected into a new array.
# =======================================================

    if ( inputType=="fastq" && defined(inputGroups) ){
        Array[InputGroup] inputGroups_ = select_first([inputGroups])
        scatter (ig in inputGroups_) {
            call bwaMem.bwaMem {
                input:
                    fastqR1 = ig.fastq1,
                    fastqR2 = ig.fastq2,
                    readGroups = ig.readGroups
            }
            BamInputs indexedBamFiles = {
                "bamFile":bwaMem.bwaMemBam,
                "baiFile":bwaMem.bwaMemIndex
            }
        }
    }

# =======================================================
#   Check to see if bam files array has 1 or 2 files.
#   Determines whether we run tumor/normal or tumor-only.
# =======================================================

    Array[BamInputs] bamFiles_ = select_first([bamFiles, indexedBamFiles])    
    
    if ( length(bamFiles_)==1 ) {
        call tumorOnlyMetrics {
            input:
                tumorBamFile = bamFiles_[0].bamFile,
                tumorBaiFile = bamFiles_[0].baiFile,
                refVCF = refVCF,
                modules =modules
        }
    }
    if ( length(bamFiles_)==2 ) {
        call getMetrics {
            input:
                tumorBamFile = bamFiles_[0].bamFile,
                tumorBaiFile = bamFiles_[0].baiFile,
                normalBamFile = bamFiles_[1].bamFile,
                normalBaiFile = bamFiles_[1].baiFile,
                refVCF = refVCF,
                modules =modules
        }
    }

    output {
        File contaminationMetrics = select_first([tumorOnlyMetrics.tumorContaminationTable, getMetrics.pairContaminationTable])
    }
}


task getMetrics{
    input {
        File normalBamFile
        File normalBaiFile
        File tumorBamFile
        File tumorBaiFile
        String refVCF
        String modules
        Int memory = 24
        Int timeout = 12
    }

    parameter_meta {
        normalBamFile: "Reference bam file"
        normalBaiFile: "Index of reference bam file"
        tumorBamFile: "Tumor bam file"
        tumorBaiFile: "Index of tumor bam file"
        memory: "Memory allocated for this job"
        timeout: "Time in hours before task timeout"
    }

    command <<<

ln -s ~{tumorBamFile} ./tumorBamFile.bam
ln -s ~{tumorBaiFile} ./tumorBamFile.bam.bai
ln ~{normalBamFile} ./normalBamFile.bam
ln ~{normalBaiFile} ./normalBamFile.bam.bai

$GATK_ROOT/bin/gatk GetPileupSummaries \
-I ./tumorBamFile.bam \
-V ~{refVCF} \
-L ~{refVCF} \
-O tumor.summaries.table

$GATK_ROOT/bin/gatk GetPileupSummaries \
-I ./normalBamFile.bam \
-V ~{refVCF} \
-L ~{refVCF} \
-O normal.summaries.table

$GATK_ROOT/bin/gatk CalculateContamination \
-I tumor.summaries.table \
-matched normal.summaries.table \
-O contamination.table

    >>>

    runtime {
        modules: "~{modules}"
        memory: "~{memory}G"
        timeout: "~{timeout}"
    }

    output {
        File pairContaminationTable = "contamination.table"
    }

    meta {
        output_meta: {
            pairContaminationTable: "Table containing contamination metrics for T/N pair"

        }
    }
}

task tumorOnlyMetrics{
    input {
        File tumorBamFile
        File tumorBaiFile
        String modules
        String refVCF
        Int memory = 24
        Int timeout = 12
    }

    parameter_meta {
        tumorBamFile: "Tumor bam file"
        tumorBaiFile: "Index of tumor bam file"
        memory: "Memory allocated for this job"
        timeout: "Time in hours before task timeout"
    }

    command <<<

ln -s ~{tumorBamFile} ./tumorBamFile.bam
ln -s ~{tumorBaiFile} ./tumorBamFile.bam.bai

$GATK_ROOT/bin/gatk GetPileupSummaries \
-I ./tumorBamFile.bam \
-V ~{refVCF} \
-L ~{refVCF} \
-O tumor.summaries.table

$GATK_ROOT/bin/gatk CalculateContamination \
-I tumor.summaries.table \
-O contamination.table

    >>>

        runtime {
            modules: "~{modules}"
            memory: "~{memory}G"
            timeout: "~{timeout}"
        }

        output {
            File tumorContaminationTable = "contamination.table"
        }

        meta {
            output_meta: {
                tumorContaminationTable: "Table containing tumor bam contamination metrics"

            }
        }
}
