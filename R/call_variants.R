#' Call Somatic/Germline structural variants using Manta

#' This function takes a pair of matched samples or single | multiple germline samples
#' and calls variants on them. Variant calling mode is established based on wheather tumor_bam is supplied or not.
#'
#' @param bin_path [REQUIRED] Path to strelka binary. Somatic or Germline.
#' @param tumor_bam [REQUIRED] Path to tumor  bam file.
#' @param normal_bam [OPTIONAL] Path to normal samples bam files.
#' @param ref_genome [REQUIRED] Path to reference genome.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param patient_id [OPTIONAL] Patient ID to name files after. If not given tumor_bam file name will be used, and if this is not given then normal_bam file will be used.
#' @param targeted [OPTIONAL] If exome/capture method. Default FALSE
#' @param threads [OPTIONAL] Number of threads per job. Default 3
#' @param verbose [DEFAULT==FALSE] Enables progress messages.
#' @export

call_sv_manta=function(bin_path=system.file("tools/manta/bin", "configManta.py", package = "RandomCode"),tumor_bam="",
normal_bam="",ref_genome="",output_dir="",patient_id="",verbose=FALSE,targeted=FALSE,threads=3){

  if(bin_path==""){
    stop("Perhaps forgot to call RandomCode::install_tools() first?")
  }

  sep="/"
  if(output_dir==""){
    sep=""
  }

  if (patient_id=="" & tumor_bam!=""){
    sample_name=ULPwgs::get_sample_name(tumor_bam)
  }else if(patient_id=="" & tumor_bam==""){
    sample_name=ULPwgs::get_sample_name(normal_bam)
  }else{
    sample_name=patient_id
  }

  if (tumor_bam!=""){
    tumor_bam=paste(" --tumorBam ",tumor_bam)
    normal_bam=paste(" --normalBam ",normal_bam)
    out_file_dir=paste0(output_dir,sep,sample_name,"_MANTA_SV_SOMATIC")
  }else{
    normal_bam=paste0(" --bam ", paste0(normal_bam,collapse=" --bam "))
    out_file_dir=paste0(output_dir,sep,sample_name,"_MANTA_SV_GERMLINE")
  }

  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  exome=""
  if (targeted){
    exome=" --exome "
  }

  if(verbose){
    print(paste("python2.7 ",bin_path, tumor_bam, normal_bam, " --referenceFasta ", ref_genome," --runDir ",out_file_dir, exome))
  }
  system(paste("python2.7 ",bin_path, tumor_bam, normal_bam, " --referenceFasta ", ref_genome," --runDir ",out_file_dir, exome))


  if(verbose){
    print(paste0("python2.7 ",out_file_dir,"/runWorkflow.py -j ",threads))
  }
  system(paste0("python2.7 ",out_file_dir,"/runWorkflow.py  -j ",threads))

  out_file_dir_sv=paste0(out_file_dir,"/results/variants/SVs")

  if (!dir.exists(out_file_dir_sv)){
      dir.create(out_file_dir_sv,recursive=TRUE)
  }

  if (tumor_bam==""){
    system(paste("cp",paste0(out_file_dir,"/results/variants/diploidSV.vcf.gz"),out_file_dir_sv))
    system(paste("gunzip",paste0(out_file_dir_sv,"/*")))
    system(paste("cp",paste0(out_file_dir,"/results/variants/diploidSV.vcf.gz"),out_file_dir_sv))
    system(paste("cp",paste0(out_file_dir,"/results/variants/diploidSV.vcf.gz.tbi"),out_file_dir_sv))
  }else{
    system(paste("cp",paste0(out_file_dir,"/results/variants/somaticSV.vcf.gz"),out_file_dir_sv))
    system(paste("gunzip",paste0(out_file_dir_sv,"/*")))
    system(paste("cp",paste0(out_file_dir,"/results/variants/somaticSV.vcf.gz"),out_file_dir_sv))
    system(paste("cp",paste0(out_file_dir,"/results/variants/somaticSV.vcf.gz.tbi"),out_file_dir_sv))
  }

}


#' Call Somatic/Germline single nucleotide variants using Strelka

#' This function takes a pair of matched samples or single | multiple germline variants
#' and calls variants on them. To specify the variant calling mode supply the right strelka binary.
#' It is recommended to supply a vcf with indel candidates. This can be generated using MANTA workflow
#'
#' @param bin_path [REQUIRED] Path to strelka binary. Somatic or Germline.
#' @param bin_path2 [REQUIRED] Path to manta binary. Somatic or Germline.
#' @param bin_path3 [REQUIRED] Path to bcftools binary.
#' @param bin_path4 [REQUIRED] Path to bgzip binary.
#' @param bin_path5 [REQUIRED] Path to tabix binary.
#' @param tumor_bam [REQUIRED] Path to tumor  bam file.
#' @param normal_bam [OPTIONAL] Path to normal samples bam files.
#' @param ref_genome [REQUIRED] Path to reference genome.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param patient_id [OPTIONAL] Patient ID to name files after. If not given tumor_bam file name will be used, and if this is not given then normal_bam file will be used.
#' @param targeted [OPTIONAL] If exome/capture method. Default FALSE
#' @param threads [OPTIONAL] Number of threads per job. Default 3
#' @param exec_options [OPTIONAL] Type of execution. local (Single node) / sge (multiple nodes). Default local.
#' @param verbose [DEFAULT==FALSE] Enables progress messages.
#' @export

call_variants_strelka=function(bin_path=system.file("tools/strelka/bin",
"configureStrelkaSomaticWorkflow.py", package = "RandomCode"),
bin_path2=system.file("tools/manta/bin", "configManta.py", package = "RandomCode"),
bin_path3=system.file("tools/bcftools/bin", "bcftools", package = "RandomCode"),
bin_path4=system.file("tools/htslib/bin", "bgzip", package = "RandomCode"),
bin_path5=system.file("tools/htslib/bin", "tabix", package = "RandomCode"),tumor_bam="",normal_bam="",ref_genome="",output_dir="",
patient_id="",verbose=FALSE,targeted=FALSE,threads=3,exec_options="local"){


  if(bin_path==""){
    stop("Perhaps forgot to call RandomCode::install_tools() first?")
  }

  sep="/"
  if(output_dir==""){
    sep=""
  }

  if (patient_id=="" & tumor_bam!=""){
    sample_name=ULPwgs::get_sample_name(tumor_bam)
  }else if(patient_id=="" & tumor_bam=="") {
    sample_name=ULPwgs::get_sample_name(normal_bam)
  }else{
    sample_name=patient_id
  }
  if (tumor_bam!=""){
    out_file_dir=paste0(output_dir,sep,sample_name,"_STRELKA_SNV_SOMATIC")
  }else{
    out_file_dir=paste0(output_dir,sep,sample_name,"_STRELKA_SNV_GERMLINE")
  }
  if (!dir.exists(out_file_dir)){
        dir.create(out_file_dir,recursive=TRUE)
  }

  call_sv_manta(bin_path=bin_path2,tumor_bam=tumor_bam,normal_bam=normal_bam,ref_genome=ref_genome,output_dir=output_dir,verbose=verbose,targeted=targeted,threads=threads,patient_id=patient_id)

  if (tumor_bam!=""){
    tumor_bam=paste(" --tumorBam ",tumor_bam)
    normal_bam=paste(" --normalBam ",normal_bam)
  }else{
    normal_bam=paste0(" --bam ", paste0(normal_bam,collapse=" --bam "))
  }

  exome=""
  if (targeted){
    exome=" --exome "
  }

  indel_candidates=""
  if (tumor_bam!=""){
    indel_candidates=paste0(" --indelCandidates ",paste0(output_dir,sep,sample_name,"_MANTA_SV_SOMATIC"),"/results/variants/candidateSmallIndels.vcf.gz")
  }


  if(verbose){
    print(paste("python2.7 ",bin_path, tumor_bam, normal_bam, " --referenceFasta ", ref_genome,indel_candidates," --runDir ",out_file_dir, exome))
  }
  system(paste("python2.7 ",bin_path, tumor_bam, normal_bam, " --referenceFasta ", ref_genome,indel_candidates," --runDir ",out_file_dir, exome))


  if(verbose){
    print(paste0("python2.7 ",out_file_dir,"/runWorkflow.py -m ",exec_options," -j ",threads))
  }
  system(paste0("python2.7 ",out_file_dir,"/runWorkflow.py -m ",exec_options," -j ",threads))

  if (tumor_bam==""){
    vcf_filter_variants(unfil_vcf=paste0(out_file_dir,"/results/variants/variants.vcf.gz"),
    bin_path=bin_path3,bin_path2=bin_path4,bin_path3=bin_path5,qual="",mq="",type="snp",verbose=verbose,output_dir=paste0(out_file_dir,"/results/variants/SNPs"))
    vcf_filter_variants(unfil_vcf=paste0(out_file_dir,"/results/variants/variants.vcf.gz"),
    bin_path=bin_path3,bin_path2=bin_path4,bin_path3=bin_path5,qual="",mq="",type="indel",verbose=verbose,output_dir=paste0(out_file_dir,"/results/variants/INDELs"))
  }else{
    if (!dir.exists(paste0(out_file_dir,"/results/variants/SNPs"))){
          dir.create(paste0(out_file_dir,"/results/variants/SNPs"),recursive=TRUE)
    }

    if (!dir.exists(paste0(out_file_dir,"/results/variants/INDELs"))){
          dir.create(paste0(out_file_dir,"/results/variants/INDELs"),recursive=TRUE)
    }
    system(paste("mv",paste0(out_file_dir,"/results/variants/somatic.snvs*"),paste0(out_file_dir,"/results/variants/SNPs")))
    system(paste("mv",paste0(out_file_dir,"/results/variants/somatic.indel*"),paste0(out_file_dir,"/results/variants/INDELs")))
  }
}
