#!/usr/bin/env nextflow

genome_file = file(params.genome_file)
name        = params.name

OUTDIR = params.outdir


// Check if paired or unpaired analysis
mode = "paired"
fastq = Channel.create()
if(params.fastq_N_R1 && params.fastq_N_R2) {	
    fastq = [['T', file(params.fastq_T_R1), file(params.fastq_T_R2)],
	     ['N', file(params.fastq_N_R1), file(params.fastq_N_R2)]]
}
else {
    fastq = [['T', file(params.fastq_T_R1), file(params.fastq_T_R2)]]
    mode = "unpaired"
}    


// Split bed file in to smaller parts to be used for parallel variant calling
Channel
    .fromPath("${params.regions_bed}")
    .ifEmpty { exit 1, "Regions bed file not found: ${params.regions_bed}" }
    .splitText( by: 150, file: 'bedpart.bed' )
    .into { beds_mutect; beds_freebayes; beds_tnscope; beds_vardict }

// Pindel bed file
if(params.pindel) {
  Channel
      .fromPath("${params.pindelbed}")
      .ifEmpty { exit 1, "Pindel regions bed file not found: ${params.pindelbed}" }
      .set { bed_pindel }
}





process bwa_align {
    cpus 20
    
    input: 
	set val(type), file(r1), file(r2) from fastq

    output:
	set val(type), file("${type}_bwa.sort.bam"), file("${type}_bwa.sort.bam.bai") into bwa_bam

    script:
    if( params.sentieon_bwa ) {
	"""
        sentieon bwa mem -M -R '@RG\\tID:${name}_${type}\\tSM:${name}_${type}\\tPL:illumina' -t ${task.cpus} $genome_file $r1 $r2 | sentieon util sort -r $genome_file -o ${type}_bwa.sort.bam -t ${task.cpus} --sam2bam -i -
        """
    }
    else {
	"""
        bwa mem -R '@RG\\tID:${name}_${type}\\tSM:${name}_${type}\\tPL:illumina' -M -t ${task.cpus} $genome_file $r1 $r2 | samtools view -Sb - | samtools sort -o ${type}_bwa.sort.bam -
        samtools index ${type}_bwa.sort.bam
        """
    }
}




process markdup {
    publishDir "${OUTDIR}/bam/myeloid", mode: 'copy', overwrite: true
    cpus 16
    memory '64 GB'
    
    input:
	set val(type), file(bam), file(bai) from bwa_bam

    output:
	set val(type), file("${name}.${type}.markdup.bam"), file("${name}.${type}.markdup.bam.bai") into bams
	set val(type), file(bam), file(bai), file("dedup_metrics.txt") into bams_qc

    """
    sentieon driver -t ${task.cpus} -i $bam --algo LocusCollector --fun score_info score.gz
    sentieon driver -t ${task.cpus} -i $bam --algo Dedup --score_info score.gz --metrics dedup_metrics.txt ${name}.${type}.markdup.bam
    """
}


// Split tumor and normal bams into different channels
bamT = Channel.create()
bamN = Channel.create()
bams.choice(bamT, bamN) {it[0] == "T" ? 0 : 1}

// Send them to the different variant callers
(bamT_freebayes, bamT_mutect, bamT_tnscope, bamT_vardict, bamT_pindel) = bamT.into(5)
bamN_freebayes = Channel.from( ["N", file("NO_FILE"), file("NO_FILE")] )
bamN_mutect    = Channel.from( ["N", file("NO_FILE"), file("NO_FILE")] )
bamN_tnscope   = Channel.from( ["N", file("NO_FILE"), file("NO_FILE")] )
bamN_vardict   = Channel.from( ["N", file("NO_FILE"), file("NO_FILE")] )
bamN_pindel    = Channel.from( ["N", file("NO_FILE"), file("NO_FILE")] )
if( mode == "paired" ) {
    (bamN_freebayes, bamN_mutect, bamN_tnscope, bamN_vardict, bamN_pindel) = bamN.into(5)
}




process vardict {
    cpus 1
    
    input:
	set val(typeT), file(bamT), file(baiT) from bamT_vardict
        set val(typeN), file(bamN), file(baiN) from bamN_vardict
        each file(bed) from beds_vardict

    output:
	set val("vardict"), file("vardict_${bed}.vcf") into vcfparts_vardict

    when:
	params.vardict
    
    script:
    if( mode == "paired" ) {
   	"""
	vardict -G $genome_file -f 0.03 -N ${name}_T -b "$bamT|$bamN" -c 1 -S 2 -E 3 -g 4 $bed | testsomatic.R | var2vcf_paired.pl -N "${name}_T|${name}_N" -f 0.03 > vardict_${bed}.vcf
        """
    }
    else if( mode == "unpaired" ) {
   	"""
	vardict -G $genome_file -f 0.03 -N ${name}_T -b $bamT -c 1 -S 2 -E 3 -g 4 $bed | teststrandbias.R | var2vcf_valid.pl -N ${name}_T -E -f 0.03 > vardict_${bed}.vcf
        """
    }
}



process freebayes {
    cpus 1
    
    input:
	set val(typeT), file(bamT), file(baiT) from bamT_freebayes
        set val(typeN), file(bamN), file(baiN) from bamN_freebayes
        each file(bed) from beds_freebayes

    output:
	set val("freebayes"), file("freebayes_${bed}.vcf") into vcfparts_freebayes

    when:
	params.freebayes

    //#cat freebayes_${bed}.filt1.vcf | awk -f /data/bnf/scripts/freebayes_somatic.awk | vcfuniq | vcffilter -A -F 'LowCov' -g "DP > 100" | awk '$10 != "." ' > freebayes_${bed}.vcf

    script:
    if( mode == "paired" ) {
   	"""
        freebayes -f $genome_file -t $bed --pooled-continuous --pooled-discrete --min-repeat-entropy 1 -F 0.03 $bamT $bamN > freebayes_${bed}.vcf.raw
        vcffilter -F LowCov -f "DP > 500" -f "QA > 1500" freebayes_${bed}.vcf.raw | vcffilter -F LowFrq -o -f "AB > 0.05" -f "AB = 0" | vcfglxgt > freebayes_${bed}.filt1.vcf
        filter_freebayes_somatic.pl freebayes_${bed}.filt1.vcf > freebayes_${bed}.vcf
        """
    }
    else if( mode == "unpaired" ) {
   	"""
        freebayes -f $genome_file -t $bed --pooled-continuous --pooled-discrete --min-repeat-entropy 1 -F 0.03 $bamT > freebayes_${bed}.vcf
        """
    }
}


process mutect {
    cpus 1
    
    input:
	set val(typeT), file(bamT), file(baiT) from bamT_mutect
        set val(typeN), file(bamN), file(baiN) from bamN_mutect
        each file(bed) from beds_mutect

    output:
	set val("mutect"), file("mutect_${bed}.vcf") into vcfparts_mutect

    
    when:
	params.mutect

    script:
    if( mode == "paired" ) {
	"""
	gatk --java-options "-Xmx2g" Mutect2 -R $genome_file -I $bamT -I $bamN -tumor ${name}_T -normal ${name}_N -L $bed -O mutect_${bed}.vcf
        """
    }
    else if( mode == "unpaired" ) {
	"""
	gatk --java-options "-Xmx2g" Mutect2 -R $genome_file -I $bamT -tumor ${name}_T -L $bed -O mutect_${bed}.vcf
        """
    }
}



process sentieon_preprocess_bam {
    cpus 10

    input:
	set val(typeT), file(bamT), file(baiT) from bamT_tnscope
        set val(typeN), file(bamN), file(baiN) from bamN_tnscope

    output:
	set file(bamT), file(baiT), file("T_recal.table") into processed_bamT_tnscope
	set file(bamN), file(baiN), file("N_recal.table") into processed_bamN_tnscope

    when:
	params.tnscope

    script:
    if( mode == "paired" ) {
	"""
        sentieon driver -t ${task.cpus} -r $genome_file -i ${bamT} --algo QualCal T_recal.table
        sentieon driver -t ${task.cpus} -r $genome_file -i ${bamN} --algo QualCal N_recal.table
        """
    }
    else if( mode == "unpaired" ) {
	"""
        sentieon driver -t ${task.cpus} -r $genome_file -i ${bamT} --algo QualCal T_recal.table
        touch ${bamN} 
        touch ${baiN} 
        touch N_recal.table
        """
    }
}


process sentieon_qc {
    cpus 40
    memory '64 GB'
    publishDir "${OUTDIR}/postmap/myeloid/", mode: 'copy', overwrite: 'true'

    input:
	set type, file(bam), file(bai), file(dedup) from bams_qc

    """
        sentieon driver \\
                --interval $params.regions_bed \\
                -r $genome_file \\
                -t ${task.cpus} -i ${bam} \\
                --algo MeanQualityByCycle mq_metrics.txt \\
                --algo QualDistribution qd_metrics.txt \\
                --algo GCBias --summary gc_summary.txt gc_metrics.txt \\
                --algo AlignmentStat aln_metrics.txt \\
                --algo InsertSizeMetricAlgo is_metrics.txt \\
                --algo CoverageMetrics --cov_thresh 1 --cov_thresh 10 --cov_thresh 30 --cov_thresh 100 --cov_thresh 250 --cov_thresh 500 cov_metrics.txt
        sentieon driver \\
                -r $genome_file \\
                -t ${task.cpus} -i ${bam} \\
                --algo HsMetricAlgo --targets_list $params.interval_list --baits_list $params.interval_list hs_metrics.txt
        qc_sentieon.pl ${name}_${type} panel > ${name}_${type}.QC
        """

}




process sentieon_tnscope {
    cpus 4
    
    input:
	set file(bamT), file(baiT), file(recal_tableT) from processed_bamT_tnscope
        set file(bamN), file(baiN), file(recal_tableN) from processed_bamN_tnscope
        each file(bed) from beds_tnscope

    output:
	set val("tnscope"), file("tnscope_${bed}.vcf") into vcfparts_tnscope

    script:
    if( mode == "paired" ) {
	"""
        sentieon driver -t ${task.cpus} -r $genome_file -i $bamT -q $recal_tableT -i $bamN -q $recal_tableN --interval $bed --algo TNscope --tumor_sample ${name}_T --normal_sample ${name}_N --clip_by_minbq 1 --max_error_per_read 3 --min_init_tumor_lod 2.0 --min_base_qual 10 --min_base_qual_asm 10 --min_tumor_allele_frac 0.00005 tnscope_${bed}.tmp.vcf
        sentieon driver -t ${task.cpus} -r $genome_file --algo TNModelApply --model $params.tnscope_model -v tnscope_${bed}.tmp.vcf tnscope_${bed}.tmp2.vcf
        bcftools filter -s "ML_FAIL" -i "INFO/ML_PROB > 0.81" tnscope_${bed}.tmp2.vcf -m x -o tnscope_${bed}.vcf
        """
    }
    else {
	"""
        sentieon driver -t ${task.cpus} -r $genome_file -i $bamT -q $recal_tableT --interval $bed --algo TNscope --tumor_sample ${name}_T --clip_by_minbq 1 --max_error_per_read 3 --min_init_tumor_lod 2.0 --min_base_qual 10 --min_base_qual_asm 10 --min_tumor_allele_frac 0.00005 tnscope_${bed}.vcf
        """ 
    }
}


// Prepare vcf parts for concatenation
vcfparts_freebayes = vcfparts_freebayes.groupTuple()
vcfparts_tnscope   = vcfparts_tnscope.groupTuple()
vcfparts_mutect    = vcfparts_mutect.groupTuple()
vcfparts_vardict   = vcfparts_vardict.groupTuple()
vcfs_to_concat = vcfparts_freebayes.mix(vcfparts_mutect, vcfparts_tnscope, vcfparts_vardict)

process concatenate_vcfs {
    publishDir "${OUTDIR}/vcf/myeloid", mode: 'copy', overwrite: true
    
    input:
	set vc, file(vcfs) from vcfs_to_concat

    output:
	set val("sample"), file("${name}_${vc}.vcf.gz") into concatenated_vcfs

    """
    vcf-concat $vcfs | vcf-sort -c | gzip -c > ${vc}.concat.vcf.gz
    vt decompose ${vc}.concat.vcf.gz -o ${vc}.decomposed.vcf.gz
    vt normalize ${vc}.decomposed.vcf.gz -r $genome_file | vt uniq - -o ${name}_${vc}.vcf.gz
    """
}



process aggregate_vcfs {
    input:
        set sample, file(vcfs) from concatenated_vcfs.groupTuple()
    
    output:
	file("${name}.agg.vcf") into agg_vcf

    """
    aggregate_vcf.pl --vcf ${vcfs.join(",")} |vcf-sort -c > ${name}.agg.vcf
    """
}



process annotate_vep {
    container = '/fs1/resources/containers/container_VEP.sif'
    publishDir "${OUTDIR}/vcf/myeloid", mode: 'copy', overwrite: true
    cpus 16
    
    input:
	file(vcf) from agg_vcf

    output:
	file("${name}.vep.vcf") into vep

    """
    vep -i ${vcf} -o ${name}.vep.vcf \\
    --offline --merged --everything --vcf --no_stats \\
    --fork ${task.cpus} \\
    --force_overwrite \\
    --plugin CADD $params.CADD --plugin LoFtool \\
    --fasta $params.VEP_FASTA \\
    --dir_cache $params.VEP_CACHE --dir_plugins $params.VEP_CACHE/Plugins \\
    --distance 200 \\
    -cache -custom $params.GNOMAD \\
    """
}
