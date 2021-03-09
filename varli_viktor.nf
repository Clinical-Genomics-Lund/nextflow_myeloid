#!/usr/bin/env nextflow

genome_file = file(params.genome_file)

OUTDIR = params.outdir+'/'+params.subdir
CRONDIR = params.crondir

csv = file(params.csv)
mode = csv.countLines() > 2 ? "paired" : "unpaired"
println(csv)
println(mode)

Channel
    .fromPath(params.csv).splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, row.type, file(row.bam), file(row.bai)) }
    .into { bam_varli }




process varlo_merge_prepro{
	cpus 5
	time '5h'


	//input:
	//	set group, file(vcf) from vcf_done

	when:
		params.varlo


	output:
		file("CandidateVariants.vcf") into candidate

	"""
	source activate py3-env
	python3 /trannel/proj/varlociraptor/bin/merge_for_varlociraptor.py --callers 'vardict,tnscope,freebayes' --dir $params.vcfs_path --output CandidateVariants.vcf
	"""
	}

process split_candidates {
	cpus 1
	time '20m'
	tag "$group"
    container = '/fs1/resources/containers/wgs_active.sif'

	input:
		file(vcf) from candidate

	output:
		file("*.parts") into candidate_parts

	'''
	grep -v '^#' CandidateVariants.vcf | split -l 1000 - --filter='sh -c "{ grep ^# CandidateVariants.vcf; cat; } > $FILE.parts"'
	'''

}
candidate_parts.flatMap().set{ cp }

process preprocess {
	cpus 1
	time '20m'
	tag "$id"
	memory '1GB'

	input:
		set group, id, type, file(bam), file(bais), file(part) from bam_varli.combine(cp)
		//each file(part) from candidate_parts

	output:
		//set group, id, type, slice, file("${id}.${type}.observations.${part}.sort.bcf.gz"), file("${id}.${type}.observations.${part}.sort.bcf.gz.csi") into vcfparts_varlo
		set group, id, type, val(slice), file("${id}.${type}.observations.${part}.bcf") into varlo_call_vcf

	script:
		pattern = part =~ /(\w+)\.(parts)/
		slice = pattern[0][1]
	"""
	source activate py3-env
    varlociraptor preprocess variants $genome_file --bam ${bam} --output ${id}.${type}.observations.${part}.bcf < $part
	bcftools sort ${id}.${type}.observations.${part}.bcf > ${id}.${type}.observations.${part}.sort.bcf
	bgzip ${id}.${type}.observations.${part}.sort.bcf
	bcftools index ${id}.${type}.observations.${part}.sort.bcf.gz
	"""

}


process varloci_calling {
	cpus 1
	time '20m' 
	//publishDir "$OUTDIR/vcf", mode :'copy'
	tag "$group"

	input:
		set group, id, type, slice, file(bcfs) from varlo_call_vcf.groupTuple(by: [0,3])

	output:
		set group, file("${group}.varloci.${slice}.calls.sort.bcf.gz"), file("${group}.varloci.${slice}.calls.sort.bcf.gz.csi") into combine_bcf

	script:
		tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }
		normal_idx = type.findIndexOf{ it == 'normal' || it == 'N' }
	
		if (mode == 'paired') {
			"""
			source activate py3-env
			varlociraptor call variants tumor-normal --purity 0.50 --tumor ${bcfs[tumor_idx]} --normal ${bcfs[normal_idx]} > ${group}.varloci.${slice}.calls.bcf
			bcftools sort ${group}.varloci.${slice}.calls.bcf > ${group}.varloci.${slice}.calls.sort.bcf
			bgzip ${group}.varloci.${slice}.calls.sort.bcf
			bcftools index ${group}.varloci.${slice}.calls.sort.bcf.gz
			"""
		}
		else {
			"""
			source activate py3-env
			varlociraptor call variants generic --scenario $params.varlo --obs tumor=$bcfs > ${group}.varloci.${slice}.calls.bcf
			bcftools sort ${group}.varloci.${slice}.calls.bcf > ${group}.varloci.${slice}.calls.sort.bcf
			bgzip ${group}.varloci.${slice}.calls.sort.bcf
			bcftools index ${group}.varloci.${slice}.calls.sort.bcf.gz
			"""
		}
}

process concatenate_vcfs_varlo {
	cpus 1
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true
	time '20m'    
	tag "$id"

	input:
		set group, file(vcfs), file(csis) from combine_bcf.groupTuple()

	output:
		set group, file("${group}.concat.bcf") into filter_varlo

	"""
	bcftools concat -a $vcfs > tmp.merged.bcf
	bcftools sort tmp.merged.bcf -O u > ${group}.concat.bcf
	"""
}


// process  FDR_filtering {

// 	publishDir "$OUTDIR/vcf", mode :'copy'
// 	tag "${smpl_id}"
	 
// 	when:
// 		params.fdr
// 	input:
// 		set group, val(id), file(calls_file) from calls_bcf
// 	output:
// 		set val(smpl_id), file ("${smpl_id}.calls.fdrFilter.bcf") into callsfiltered
// 	script:
// 	"""
// 	varlociraptor filter-calls control-fdr ${calls_file} --events SOMATIC_TUMOR --fdr 0.05 --var SNV > ${smpl_id}.calls.fdrFilter.bcf
// 	"""
// }

// process posteriorOdds_filtering {
// 	publishDir "$OUTDIR/vcf", mode :'copy'
// 	tag "${smpl_id}"
	
// 	when:
// 		params.posteriorOdds
		
// 	input:
//                 set val(smpl_id), file(calls_file) from calls_bcf2
// 	output:
//                 set val(smpl_id), file ("${smpl_id}.calls.posteriorFilter.bcf") into callsfiltered_post
//     script:
// 	"""
// 	varlociraptor filter-calls posterior-odds --events SOMATIC_TUMOR --odds strong < ${calls_file} > ${smpl_id}.calls.posteriorFilter.bcf  
// 	"""
// }