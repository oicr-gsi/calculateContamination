version 1.0

import "imports/pull_bwaMem.wdl" as bwaMem

workflow calculateContamination {
    input {
        File fastq1
        File fastq2
        File? bamFile
    }

    parameter_meta {
    }

    meta {
        author: "Murto Hilali"
        email: "mhilali@oicr.on.ca"
        description: ""
        dependencies: [
        ]
        output_meta: {
        }
    }

    if ( defined(bamFile) == false ){
        call bwaMem.bwaMem {
            input:
                fastqR1 = fastq1,
                fastqR2 = fastq2,
                runBwaMem_bwaRef = "foo",
                runBwaMem_modules = "foo",
                readGroups = "foo"
        }
        File bamFile = bwaMem.bwaMem.outputBam

    }



    call hello { 
        input:
    }

    call hello as taskName {
        input:
    }

    output {
    } 
}

task hello {
        input {
            String modules = ""
            Int memory = 24
            Int timeout = 12
        }

        parameter_meta {
            modules: "Required environment modules"
            memory: "Memory allocated for this job"
            timeout: "Time in hours before task timeout"
        }

        command <<<

        >>>

        runtime {
            modules: "~{modules}"
            memory: "~{memory}G"
            timeout: "~{timeout}"
        }

        output {

        }

        meta {
            output_meta: {
            }
        }
}