#!/usr/bin/dev nextflow

params.ref = '/fs1/resources/ref/hg19/fasta/human_g1k_v37_decoy.fasta'
//params.candidates_bcf = '/data/bnf/dev/sima/varloci/candidates.vcf'
params.purity = 0.75
//params.smpl_id ='KTC-HD829-C2'
smpl_id ='KTC-HD829-C2'
params.normal_bam='/fs1/results/myeloid/bam/KTC-HD829-C2.normal.dedup.bam'
params.tumor_bam ='/fs1/results/myeloid/bam/KTC-HD829-C2.tumor.dedup.bam'
params.outdir ='/data/bnf/dev/sima/varloci'
OUTDIR = params.outdir
params.vcfs_path = '/fs1/results/myeloid/vcf/'



process merge_vcfs{

	publishDir "$OUTDIR/vcf", mode :'copy'

	output:
		file "CandidateVariants.vcf" into candidate_vcf
	
	"""
	python /data/bnf/dev/sima/varloci/bin/merge_for_varlociraptor.py --callers 'vardict,tnscope,freebayes' --dir ${params.vcfs_path} --output CandidateVariants.vcf
	"""
	}

process  varloci_prepro {

	publishDir "$OUTDIR/vcf", mode :'copy'
	tag "${smpl_id}"
	
	input:
		file candidates from candidate_vcf
		//file(normal.bam)		
		//file(tumor_bam)
	
	output:
		set val(smpl_id), file("${smpl_id}.tumor.observations.bcf"), file("${smpl_id}.normal.observations.bcf") into tum_normal_obs

	 
	
	"""
	 varlociraptor preprocess variants ${params.ref} --bam ${params.normal_bam} --output ${smpl_id}.normal.observations.bcf < ${candidates}
	 varlociraptor preprocess variants ${params.ref} --bam ${params.tumor_bam} --output ${smpl_id}.tumor.observations.bcf < ${candidates}
	"""
}




process varloci_calling_TN {

	publishDir "$OUTDIR/vcf", mode :'copy'
	tag "${smpl_id}"
	//when:
	//	params.TN
	input:
		set val(smpl_id), file(tumor_bcf), file(normal_bcf) from tum_normal_obs
	output:
		set val(smpl_id), file("${smpl_id}.varloci.calls.bcf") into calls_bcf, calls_bcf2
	script:
	
	"""
	varlociraptor call variants tumor-normal --purity ${params.purity} --tumor ${tumor_bcf} --normal ${normal_bcf} --output ${smpl_id}.varloci.calls.bcf
	"""
}
process varlociCalling_generic {

	publishDir "$OUTDIR/vcf", mode :'copy'
	tag "${smpl_id}"
	//when:
	//	params.generic
	input:
		set val(smpl_id), file(tumor_bcf) from tum_obs
	output:
		set val(smpl_id), file("${smpl_id}.varloci.calls.bcf") into calls_bcf, calls_bcf2
	script:
	//tumor is the name of sample given in the scenario.
	"""
	varlociraptor call variants  generic --scenario scenario.yaml --obs tumor=${tumor_bcf} > ${smpl_id}.varloci.calls.bcf
	"""
}

//based on the given event in scenario you change the events param in filtering processes down here:
process  FDR_filtering {

	publishDir "$OUTDIR/vcf", mode :'copy'
	tag "${smpl_id}"
	 
	when:
		params.fdr
	input:
		set val(smpl_id), file(calls_file) from calls_bcf
	output:
		set val(smpl_id), file ("${smpl_id}.calls.fdrFilter.bcf") into callsfiltered
	script:
	"""
		varlociraptor filter-calls control-fdr ${calls_file} --events SOMATIC_TUMOR --fdr 0.05 --var SNV > ${smpl_id}.calls.fdrFilter.bcf
	"""
	}

process posteriorOdds_filtering {
	publishDir "$OUTDIR/vcf", mode :'copy'
	tag "${smpl_id}"
	
	when:
		params.posteriorOdds
		
	input:
                set val(smpl_id), file(calls_file) from calls_bcf2
	output:
                set val(smpl_id), file ("${smpl_id}.calls.posteriorFilter.bcf") into callsfiltered_post
    script:
	"""
	varlociraptor filter-calls posterior-odds --events SOMATIC_TUMOR --odds strong < ${calls_file} > ${smpl_id}.calls.posteriorFilter.bcf  
	"""
	 }




