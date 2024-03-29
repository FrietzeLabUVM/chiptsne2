#### Example data ####

#' exampleQueryGR
#'
#' @return GRanges example
#' @export
#'
#' @importFrom GenomicRanges mcols
#' @examples
#' exampleQueryGR()
exampleQueryGR = function(){
    query_gr = seqsetvis::CTCF_in_10a_overlaps_gr
    colnames(GenomicRanges::mcols(query_gr)) = paste0(
        "peak_",
        colnames(GenomicRanges::mcols(query_gr))
    )
    query_gr
}
#' exampleProfDT
#'
#' @return data.table example
#' @export
#'
#' @examples
#' exampleProfDT()
exampleProfDT = function(){
    seqsetvis::CTCF_in_10a_profiles_dt
}

#' query_gr
#'
#' @return ChIPtsne object for testing
#' @export
#'
#' @examples
#' exampleChIPtsne2()
exampleChIPtsne2 = function(){
    query_gr = exampleQueryGR()
    prof_dt = exampleProfDT()

    ChIPtsne2.from_tidy(prof_dt, query_gr)
}

#' exampleChIPtsne2.with_meta
#'
#' @return ChIPtsne object for testing, includes meta data
#' @export
#' @importFrom tidyr separate
#'
#' @examples
#' exampleChIPtsne2.with_meta()
exampleChIPtsne2.with_meta = function(){
    query_gr = exampleQueryGR()
    prof_dt = exampleProfDT()
    meta_dt = prof_dt %>%
        dplyr::select(sample) %>%
        unique %>%
        tidyr::separate(sample, c("cell", "mark"), sep = "_", remove = FALSE)

    ChIPtsne2.from_tidy(prof_dt, query_gr, sample_metadata = meta_dt)
}

#' exampleBamFiles
#'
#' @return character file paths to example bam files.
#' @export
#'
#' @examples
#' exampleBamFiles()
exampleBamFiles = function(){
    dir(system.file("extdata", package = "chiptsne2"), pattern = "MCF.+bam$", full.names = TRUE)
}

#' exampleBigWigFiles
#'
#' @return character file paths to example bigWig files.
#' @export
#'
#' @examples
#' exampleBigWigFiles()
exampleBigWigFiles = function(){
    dir(system.file("extdata", package = "seqsetvis"), pattern = "bw$", full.names = TRUE)
}

#' exampleBamConfigFile
#'
#' @return character file path to example bam config file.
#' @export
#'
#' @examples
#' exampleBamConfigFile()
exampleBamConfigFile = function(){
    system.file("extdata/bam_config.csv", package = "chiptsne2", mustWork = TRUE)
}

#' exampleBigWigConfigFile
#'
#' @return character file path to example bigWig config file.
#' @export
#'
#' @examples
#' exampleBigWigConfigFile()
exampleBigWigConfigFile = function(){
    system.file("extdata/bigwig_config.csv", package = "chiptsne2", mustWork = TRUE)
}

#' exampleBam_data.frame
#'
#' @return data.frame with example bam files and metadata.
#' @export
#'
#' @examples
#' exampleBam_data.frame()
exampleBam_data.frame = function(){
    bam_config_file = exampleBamConfigFile()
    bam_config_df = .parse_config_body(bam_config_file)
    pkg_dir = system.file(package = "chiptsne2", "extdata", mustWork = TRUE)
    bam_config_df$file = file.path(pkg_dir, bam_config_df$file)
    bam_config_df


}

#' exampleBigWig_data.frame
#'
#' @return data.frame with example bigWig files and metadata.
#' @export
#'
#' @examples
#' exampleBigWig_data.frame()
exampleBigWig_data.frame = function(){
    bigwig_config_file = exampleBigWigConfigFile()
    bigwig_config_df = .parse_config_body(bigwig_config_file)
    pkg_dir = system.file(package = "chiptsne2", "extdata", mustWork = TRUE)
    bigwig_config_df$file = file.path(pkg_dir, bigwig_config_df$file)
    bigwig_config_df
}
