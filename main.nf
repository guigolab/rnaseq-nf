params.reads = "$baseDir/data/*_{1,2}.fastq.gz"
params.genome = "$baseDir/data/genome.fa"
params.annotation = "$baseDir/data/annotation.gtf"
params.outdir = 'results'
params.units = [ 'expected_count', 'TPM' ]

genomeFile = file(params.genome)
annotationFile = file(params.annotation)
units = [ params.units ].flatten()

Channel
    .fromFilePairs(params.reads, checkExists:true)
    .ifEmpty{ error("No input files found") }
    .set { readsChannel }

process genomeIndex {
    input:
    file genome from genomeFile
    file annotation from annotationFile

    output:
    file 'genomeIndex' into genomeIndexChannel

    """
    mkdir -p genomeIndex
    STAR --runThreadN ${task.cpus} \\
         --runMode genomeGenerate \\
         --genomeDir genomeIndex \\
         --genomeFastaFiles $genome \\
         --sjdbGTFfile $annotation \\
         --genomeSAindexNbases 11
    """
}

process transcriptomeIndex {
    input:
    file genome from genomeFile
    file annotation from annotationFile

    output:
    file 'transcriptomeIndex' into transcriptomeIndexChannel

    """
    mkdir -p transcriptomeIndex
    rsem-prepare-reference --gtf $annotation \\
                           $genome \\
                           transcriptomeIndex/RSEMref
    """
}

process mapping {
    tag { prefix }
    publishDir { "${params.outdir}/${task.process}" }

    input:
    file index from genomeIndexChannel
    set prefix, file(reads) from readsChannel

    output:
    set prefix, file('*.sortedByCoord.out.bam') into genomeAlignmentsChannel
    set prefix, file('*.toTranscriptome.out.bam') into transcriptomeAlignmentsChannel

    """
    STAR --runThreadN ${task.cpus} \\
         --genomeDir $index \\
         --readFilesIn $reads \\
         --outSAMunmapped Within \\
         --outFilterType BySJout \\
         --outSAMattributes NH HI AS NM MD \\
         --readFilesCommand pigz -p${task.cpus} -dc \\
         --outSAMtype BAM SortedByCoordinate \\
         --quantMode TranscriptomeSAM \\
         --outFileNamePrefix ${prefix}_
    """
}

process quantification {
    tag { prefix }
    publishDir { "${params.outdir}/${task.process}" }

    input:
    file transcriptomeIndex from transcriptomeIndexChannel
    set prefix, file(transcriptomeAlignments) from transcriptomeAlignmentsChannel

    output:
    file "${prefix}.genes.results" into geneQuantificationChannel
    file "${prefix}.isoforms.results" into isoformQuantificationChannel

    """
    rsem-calculate-expression --num-threads ${task.cpus} \\
                              --bam \\
                              --paired-end \\
                              --estimate-rspd \\
                              --forward-prob 0 \\
                              --no-bam-output \\
                              --seed 12345 \\
                              $transcriptomeAlignments \\
                              $transcriptomeIndex/RSEMref \\
                              $prefix
    """
}

process matrix {
    tag { unit }
    publishDir { "${params.outdir}/${task.process}" }

    input:
    each unit from units
    file quantification from geneQuantificationChannel.collect().sort{ it.simpleName }

    output:
    file outputMatrix

    script:
    files = quantification.collect { "'${it.name}'" }.join(',')
    outputMatrix = "mouse.gene.matrix.${unit}.tsv"
    """
    #!/usr/bin/env python
    import csv
    import os

    d = {}
    samples = set()

    for f in [ ${files} ]:
        with open(f) as fd:
            sample = os.path.basename(f).split('.')[0]
            samples.add(sample)
            csvfile = csv.DictReader(fd, delimiter='\t')
            for line in csvfile:
                element_id = line['gene_id']
                element_value_list = []
                v = line['${unit}']
                element_value = float(v) if v != "NA" else "NA"
                element_value_list += [str(element_value)]
                d.setdefault(element_id, {}).setdefault(sample, ",".join(element_value_list))

    f = open('${outputMatrix}','w')

    f.write('\\t'.join(sorted(samples, key=lambda x: x.lower()))+'\\n')
    for element in sorted(d.iterkeys()):
        values = d[element]
        f.write(element+'\\t')
        missing_value = 'NA'
        f.write('\\t'.join(str(values.get(s, missing_value)) for s in sorted(samples, key=lambda x: x.lower()))+'\\n')
    f.close()
    """

}