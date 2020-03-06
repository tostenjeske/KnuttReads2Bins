##
## Snakefile - The main executable Snakefile for KnuttReads2Bins
##
## Knutt.org/Knutt2Reads2Bins


# It is just a proxy for the Snakefile_KnuttReads2Bins. This is done to
# make the workflow callable without speciying the Snakefile location
# and also allow the "include"ing of this workflow without conflicts
#  on the "all" rule




include: "Snakefile_0KnuttReads2Bins"

rule paper:
   input:
      rules.rawQC.input,
      lambda w: [fastQC_for_file(f) for f in expand(rules.cutadapt_paired_reads.output,
                                                    sample=sample_names)],
      rules.merge.input,
      lambda w: [fastQC_for_file(f) for f in expand(merge_res_file+"{mergefile}.fastq.gz",
                                                    mergefile=["merged","unmgd_R1",
                                                               "unmgd_R2","merged_qtr",
                                                               "unmgd_R1_qtr"],
                                                    sample=sample_names,trimmed=adpt_poss)],
      lambda w: [fastQC_for_file(f) for f in rules.annotationReads.input],
      rules.classifySSUreads.input,
      rules.classifyReads.input,
      rules.annotateReads.input,
      rules.assemble.input,
      rules.metaQUAST.input,
      rules.coverageContigs.input,
      rules.checkM.input,
      rules.adapterTrimmingData.output,
      expand(rules.mergeData.output,trimmed=adpt_poss[0] if config["adaptertrim"] else adpt_poss[1]),
      expand(rules.mergeInsertData.output,trimmed=adpt_poss[0] if config["adaptertrim"] else adpt_poss[1]),
      expand(rules.qulitytrim_data.output,trimmed=adpt_poss[0] if config["adaptertrim"] else adpt_poss[1]),
      rules.readclass_SSU.output,
      rules.kaiju_data.output,
      expand(rules.annodata.output,db=customdbs+integrateddbs),
      rules.metaquastdata.output,
      rules.covdata.output,
      rules.covdata_details.output,
      rules.jgi_depthdata.output,
      rules.checkm_lineage.output,
      rules.checkmprofile.output,
      expand(rules.krona.output,krona_type=("readclass_SSU","readclass_kaiju")),
      expand(rules.kronawithdb.output,dbswithkrona=integrateddbs+customdbs),
