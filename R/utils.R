#' Get the name of a sample from input file name
#'
#' This function takes the absolute/relative path to a file and
#' returns its base name without the file extension suffix.
#'
#' @param file_path Path to the input file
#' @return A string with the name of the file
#' @export


get_sample_name=function(file_path=""){
  sample_name=unlist(strsplit(basename(file_path),"\\."))[1]
  return(sample_name)
}


#' Get the extension of a file
#'
#' This function takes the absolute/relative path to a file and
#' returns the file extension suffix.
#'
#' @param file_path Path to the input file
#' @return A string with the extension of the file
#' @export
get_file_extension=function(file_path=""){
    ext = strsplit(basename(file_path), split="\\.")[[1]]
    ext = paste(ext[-1],collapse=".")
    return(ext)
}

#' Get sample name from two sample replicates
#'
#' This function takes the absolute/relative path to two files
#' and returns the longest common string among their basenames
#'
#' @param file_path Path to the input file
#' @param file_path2 Path to the second input file
#' @return A string with the longest common basename
#' @export

intersect_sample_name=function(file_path="",file_path2=""){
  tmp_name=get_sample_name(file_path2)
  sample_name=sapply(sapply(c(0:(nchar(tmp_name)-1)),function (i) substr(tmp_name,1,nchar(tmp_name)-i)),function (x) grepl(x,file_path))
  sample_name=names(which(sample_name)[1])
  sample_name=sub("(.*)[_.-].*","\\1",sample_name)
  return(sample_name)

}





#' VCF filtering using bcftools
#'
#' This function filters VCF calls using bcftools
#'
#' @param bin_path Path to gatk binary. Default tools/gatk/gatk.
#' @param bin_path2 Path to bgzip binary. Default tools/htslib/bgzip.
#' @param bin_path3 Path to tabix binary. Default tools/htslib/tabix.
#' @param unfil_vcf Path to unfiltered vcf file.
#' @param qual Quality filter. Default 30.
#' @param mq Mapping quality filter. Default 40.
#' @param ref [Optional] Filter variants with ID.
#' @param state [Optional] Variant state to select. Options: het/homo
#' @param type [Optional] Variant type to include. Options: snp/indel. Default includes both
#' @param filter [Optional] Filter to include Options: PASS/.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export


vcf_filter_variants=function(unfil_vcf="",bin_path=system.file("inst/tools/bcftools/bin","tabix",
package = "RandomCode"),bin_path2=system.file("inst/tools/htslib/bin","bgzip",
package = "RandomCode"),bin_path3=system.file("inst/tools/htslib/bin","tabix",
package = "RandomCode"),qual=30,mq=40,state="",ref="",type="",filter="",verbose=FALSE,output_dir=""){

  sep="/"
  if(output_dir==""){
    sep=""
  }
  sample_name=ULPwgs::get_sample_name(unfil_vcf)
  out_file_dir=paste0(output_dir,sep,sample_name,"_FILTERED")
  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir,recursive=TRUE)
  }

  out_file=paste0(out_file_dir,"/",sample_name,".FILTERED.vcf")

  if (qual!=""){
    qual=paste0("\'%QUAL>",qual)
  }else{
    if (mq!=""){
      mq=paste0("\'%MQ>",mq)
    }else{
        if (filter!=""){
          filter=paste0("\'%FILTER=\"",filter,"\" ")
        }else{
          if (type!=""){
            type=paste0("\'%TYPE=\"",type,"\" ")
          }else{
            if (ref!=""){
              ref=paste0("\'%ID!=\"",ref,"\" ")
            }else{
              if (state!=""){
                state=paste0("\'GT[0]=\"",state,"\" ")
              }
            }
          }
        }
      }
    }


    if (state!=""&!grepl("GT",state)){
      state=paste0(" & GT[0]=\"",state,"\" ")
    }

    if (ref!=""&!grepl("%",ref)){
      ref=paste0(" & ID!=\"",ref,"\" ")
    }

    if (type!=""&!grepl("%",type)){
      type=paste0(" & TYPE=\"",type,"\" ")
    }
    if (filter!=""&!grepl("%",filter)){
      filter=paste0(" & FILTER=\"",filter,"\" ")
    }
    if (mq!=""&!grepl("%",mq)){
      mq=paste0(" & MQ>",mq)
    }

  if(verbose){
    print(paste(bin_path,"view  -i ",qual,filter,state,type,ref,mq,"\'",unfil_vcf,">",out_file))
  }
  system(paste(bin_path,"view  -i ",qual,filter,state,type,ref,mq,"\'",unfil_vcf,">",out_file))
  system(paste("cp", out_file, paste0(out_file,".tmp")))
  bgzip(bin_path=bin_path2,file=out_file)
  tab_indx(bin_path=bin_path3,file=paste0(out_file,".gz"))
  system(paste("cp", paste0(out_file,".tmp"), out_file))
  system(paste("rm -rf", paste0(out_file,".tmp")))
}



#' Index tab separated genomic regions
#'
#' This function takes a .vcf.bgz tab separated genomic region file and generates an index for it
#'
#' @param bin_path [REQUIRED] Path to bcftools binary. Default tools/htslib/tabix.
#' @param file [REQUIRED] Path to VCF file to index.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @export

tab_indx=function(bin_path=system.file("inst/tools/htslib/bin","tabix",
package = "RandomCode"),file="",verbose=FALSE){

  if (verbose){
    print(paste0(bin_path," -fp vcf",file))
  }
  system(paste0(bin_path," -fp vcf ",file))
}



#' Bgzips a VCF file
#'
#' This function takes a VCF file and Bgzips it
#'
#' @param bin_path [Required] Path to bcftools binary. Default tools/htslib/tabix.
#' @param file [Required] Path to VCF file to bgzip.
#' @param verbose [Optional] Enables progress messages. Default False.
#' @export



bgzip=function(bin_path=system.file("inst/tools/htslib/bin","bgzip",
package = "RandomCode"),file="",verbose=FALSE){
  if (verbose){
    print(paste0(bin_path," -f ",file))
  }
  system(paste0(bin_path," -f ",file))
}


#' Install all tools locally
#'
#' This function will install all required tools
#'
#' @export



install_tools=function(){

  strelka=system.file("inst/tools","strelka-2.9.10.centos6_x86_64.tar.bz2", package = "RandomCode")
  manta=system.file("inst/tools","manta-1.6.0.centos6_x86_64.tar.bz2", package = "RandomCode")
  bcftools=system.file("inst/tools","bcftools-1.14.tar.bz2", package = "RandomCode")
  htslib=system.file("inst/tools","htslib-1.14.tar.bz2", package = "RandomCode")
  system(paste0("tar -xzvf ",strelka, " && mv ",strelka, strsplit(strelka,"-")[[1]][1]))
  system(paste0("tar -xzvf ",manta, " && mv ",manta, strsplit(manta,"-")[[1]][1]))
  system(paste0("tar -xzvf ",bcftools, " && mv ",bcftools, strsplit(bcftools,"-")[[1]][1]))
  system(paste0("tar -xzvf ",htslib, " && mv ",htslib, strsplit(htslib,"-")[[1]][1]))

  system(paste0("cd ",strsplit(bcftools,"-")[[1]][1],"; ./configure --prefix=",
  strsplit(bcftools,"-")[[1]][1],"/bin",";make;make install"))

  system(paste0("cd ",strsplit(htslib,"-")[[1]][1],"; ./configure --prefix=",
  strsplit(htslib,"-")[[1]][1],"/bin",";make;make install"))

}
