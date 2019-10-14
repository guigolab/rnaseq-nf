params.reads = "$baseDir/data/*_{1,2}.fastq.gz"
params.genome = "$baseDir/data/genome.fa"
params.annotation = "$baseDir/data/annotation.gtf"

genomeFile = file(params.genome)
annotationFile = file(params.annotation)

Channel
    .fromFilePairs(params.reads)
    .set { readsChannel }

process genomeIndex {
    input:
    file genome from genomeFile
    file annotation from annotationFile

    output:
    file 'genomeIndex' into genomeIndexChannel

    """
    mkdir -p genomeIndex
    STAR --runThreadN ${task.cpus} \
         --runMode genomeGenerate \
         --genomeDir genomeIndex \
         --genomeFastaFiles $genome \
         --sjdbGTFfile $annotation \
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
    rsem-prepare-reference --gtf $annotation \
                           $genome \
                           transcriptomeIndex/RSEMref
    """
}

process mapping {
    tag { prefix }

    input:
    file index from genomeIndexChannel
    set prefix, file(reads) from readsChannel

    output:
    set prefix, file('*.sortedByCoord.out.bam') into genomeAlignmentsChannel
    set prefix, file('*.toTranscriptome.out.bam') into transcriptomeAlignmentsChannel

    """
    STAR --runThreadN ${task.cpus} \
         --genomeDir $index \
         --readFilesIn $reads \
         --outSAMunmapped Within \
         --outFilterType BySJout \
         --outSAMattributes NH HI AS NM MD \
         --readFilesCommand pigz -p${task.cpus} -dc \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM \
         --outFileNamePrefix ${prefix}_
    """
}

process quantification {
    input:
    file transcriptomeIndex from transcriptomeIndexChannel
    set prefix, file(transcriptomeAlignments) from transcriptomeAlignmentsChannel

    """
    rsem-calculate-expression --num-threads ${task.cpus} \
                              --bam \
                              --paired-end \
                              --estimate-rspd \
                              --forward-prob 0 \
                              --no-bam-output \
                              --seed 12345 \
                              $transcriptomeAlignments \
                              $transcriptomeIndex/RSEMref \
                              $prefix
    """
}