#' ChIPtsne2.from_tidy
#'
#' @param prof_dt Profile data.table, as returned from seqsetvis::ssvFetch* functions
#' @param query_gr The query GRanges object used to fetch prof_dt.
#' @param name_VAR Variable name that contains sample ids/names. Links prof_dt to meta_dt. Default is "sample".
#' @param position_VAR Variable name that contains positional information in prof_dt. Default is "x".
#' @param value_VAR Variable name that contains signal value information in prof_dt. Default is "y".
#' @param region_VAR Variable name that contains region ID information in prof_dt. Default is "id".
#' @param sample_metadata Metadata for entries in prof_dt's name_VAR, must include name_VAR
#' @param region_metadata Metadata to append to rowRanges, mcols of query_gr will also be used.
#' @param auto_sample_metadata If true, additional attributes in prof_dt will used for metadata (minus certain region related attributes such as seqnames, start, end, etc.)
#' @param obj_history Existing history for object, may be used to describe origin on prof_dt or query_gr.
#' @param fetch_config A FetchConfig object. Optional. Allows signal to be fetched as needed.
#' @param init If TRUE, initialize history with birthday, session_info, and chiptsne2_version
#'
#' @return ChIPtsne2 object with tsne results
#' @importFrom utils sessionInfo
#' @importFrom dplyr any_of all_of
#' @export
#'
#' @examples
#' query_gr = seqsetvis::CTCF_in_10a_overlaps_gr
#' prof_dt = seqsetvis::CTCF_in_10a_profiles_dt
#' ct2 = ChIPtsne2.from_tidy(prof_dt, query_gr)
#' ct2
ChIPtsne2.from_tidy = function(prof_dt,
                               query_gr,
                               sample_metadata = NULL,
                               region_metadata = NULL,
                               name_VAR = "sample",
                               position_VAR = "x",
                               value_VAR = "y",
                               region_VAR = "id",
                               auto_sample_metadata = TRUE,
                               obj_history = list(),
                               fetch_config = FetchConfig.null(),
                               init = TRUE
){
    if(init){
        init_history = list(birthday = date(), session_info = utils::sessionInfo(), chiptsne2_version = utils::packageDescription("chiptsne2")$Version)
        obj_history = c(init_history, obj_history)
    }

    if(is(prof_dt, "GRanges")){
        prof_dt = data.table::as.data.table(prof_dt)
    }
    #basic VAR checks
    if(!all(c(name_VAR, position_VAR, value_VAR, region_VAR) %in% colnames(prof_dt))){
        missed = !c(name_VAR, position_VAR, value_VAR, region_VAR) %in% colnames(prof_dt)
        missed_msg = paste(c("name_VAR", "position_VAR", "value_VAR", "region_VAR")[missed],
                           c(name_VAR, position_VAR, value_VAR, region_VAR)[missed], sep = ": ")
        stop(paste(c("Missing required VAR in prof_dt:", missed_msg), collapse = "\n"))
    }
    if(!is.null(sample_metadata)){
        if(is.null(rownames(sample_metadata))){
            if(!name_VAR %in% colnames(sample_metadata)){
                stop("name_VAR: ", name_VAR, " not found in sample_metadata.")
            }
        }else{
            if(!name_VAR %in% colnames(sample_metadata)){
                sample_metadata[[name_VAR]] = rownames(sample_metadata)
            }
        }
    }
    #determine colname order
    cn = NULL
    # use sample_metadata over prof_dt if provided and levels of factor if set, otherwise current order of character
    if(!is.null(sample_metadata)){
        if(is.factor(sample_metadata[[name_VAR]])){
            sample_metadata[[name_VAR]] = droplevels(sample_metadata[[name_VAR]])
            cn = levels(sample_metadata[[name_VAR]])
        }else if(is.character(sample_metadata[[name_VAR]])){
            cn = unique(sample_metadata[[name_VAR]])
        }
    }else{
        if(is.factor(prof_dt[[name_VAR]])){
            prof_dt[[name_VAR]] = droplevels(prof_dt[[name_VAR]])
            cn = levels(prof_dt[[name_VAR]])
        }else if(is.character(prof_dt[[name_VAR]])){
            cn = unique(prof_dt[[name_VAR]])
        }
    }
    if(is.null(cn)){
        stop("Could not determine column order from prof_dt (or sample_metadata if provided) from name_VAR: ", name_VAR, "\nIs ", name_VAR, " present and a character or factor?")
    }

    #determine rownames order
    rn = NULL
    # use prof_dt levels of factor if set, otherwise current order of character
    if(!is.null(prof_dt)){
        if(is.factor(prof_dt[[region_VAR]])){
            prof_dt[[region_VAR]] = droplevels(prof_dt[[region_VAR]])
            rn = levels(prof_dt[[region_VAR]])
        }else if(is.character(prof_dt[[region_VAR]])){
            rn = unique(prof_dt[[region_VAR]])
        }
    }
    if(is.null(rn)){
        stop("Could not determine row order from prof_dt from region_VAR: ", region_VAR, "\nIs ", region_VAR, " present and a character or factor?")
    }
    if(!is.null(region_metadata)){
        if(region_VAR %in% colnames(region_metadata)){
            rownames(region_metadata) = region_metadata[[region_VAR]]
            region_metadata[[region_VAR]] = NULL
        }
        if(is.null(rownames(region_metadata))){
            stop("region_metadata must have rownames or contain region_VAR: ", region_VAR)
        }
        if(!all(rownames(region_metadata) %in% rn)){
            stop("rownames of region_metadata are not consistent with profile region_VAR: ", region_VAR)
        }
    }

    #create wide profile matrix
    tmp_wide = tidyr::pivot_wider(
        prof_dt,
        names_from = dplyr::all_of(c(name_VAR, position_VAR)),
        values_from = dplyr::all_of(c(value_VAR)),
        id_cols = dplyr::all_of(c(region_VAR))
    )
    new_rowToRowMat = as.matrix(tmp_wide[, -1])
    rownames(new_rowToRowMat) = tmp_wide[[region_VAR]]
    new_rowToRowMat = new_rowToRowMat[rn,, drop = FALSE]

    tmp = c("value_VAR")
    names(tmp) = value_VAR

    if(is.null(sample_metadata)){
        if(auto_sample_metadata){
            drop_vars = c("seqnames",
                          "start",
                          "end",
                          "width",
                          "strand",
                          "id",
                          "y",
                          "x",
                          "cluster_id",
                          position_VAR,
                          value_VAR)
            sample_metadata = prof_dt %>%
                # dplyr::select(all_of(c(name_VAR))) %>%
                dplyr::select(!dplyr::any_of(c(drop_vars))) %>%
                unique
            if(nrow(sample_metadata) != length(cn)){
                stop("Something has gone wrong attempting to automatically derive sample_metadata. Either supply explicitly or disable with auto_sample_metadata = FALSE")
            }
        }else{
            sample_metadata = data.frame(V1 = cn, check.names = FALSE)
            colnames(sample_metadata) = name_VAR
        }
    }
    sample_metadata = data.frame(sample_metadata, check.names = FALSE)
    rownames(sample_metadata) = sample_metadata[[name_VAR]]
    sample_metadata[[name_VAR]] = NULL

    map_dt = prof_dt %>%
        dplyr::select(dplyr::all_of(c(name_VAR, position_VAR))) %>%
        unique
    map_dt = dplyr::mutate(map_dt, cn = paste(get(name_VAR), get(position_VAR), sep = "_"))
    stopifnot(map_dt$cn == colnames(new_rowToRowMat))
    map_dt = dplyr::mutate(map_dt, nr = seq(nrow(map_dt)))
    new_colToRowMatCols = split(map_dt$cn, map_dt[[name_VAR]])

    #impose cn and rn

    new_rowToRowMat = new_rowToRowMat[rn, , drop = FALSE]
    new_colToRowMatCols = new_colToRowMatCols[cn]
    sample_metadata = sample_metadata[cn, , drop = FALSE]

    new_prof_max_mat = .recalculateMax(new_rowToRowMat, new_colToRowMatCols)
    new_prof_max_mat = new_prof_max_mat[rn, cn, drop = FALSE]

    if(is.null(query_gr)){
        region_metadata = region_metadata[rn,]
        ChIPtsne2_no_rowRanges(
            assay = list(max = new_prof_max_mat),
            rowToRowMat = new_rowToRowMat,
            colToRowMatCols = new_colToRowMatCols,
            colData = sample_metadata,
            rowData = region_metadata,
            name_VAR = name_VAR,
            position_VAR = position_VAR,
            value_VAR = value_VAR,
            region_VAR = region_VAR,
            fetch_config = fetch_config,
            metadata = obj_history)
    }else{
        if(is.null(names(query_gr))){
            stop("names() must be set on query_gr. Maybe call seqsetvis::prepare_fetch_GRanges_names on query_gr before fetching prof_dt?")
        }
        if(!all(unique(as.character(prof_dt[[region_VAR]])) %in% names(query_gr))){
            stop("region_VAR: ", region_VAR, " in prof_dt is not consistent with names of query_gr")
        }
        if(!is.null(region_metadata)){
            #merge region_metadata into query_gr
            #region_metadata is not used after
            query_gr = .add_region_metadata(query_gr, region_metadata, region_VAR, overwrite = TRUE)
        }
        if(!setequal(names(query_gr), rn)){
            stop("names(query_gr) is not consistent with region_VAR: ", region_VAR, " in prof_dt")
        }
        query_gr = query_gr[rn]
        ChIPtsne2(
            assay = list(max = new_prof_max_mat),
            rowRanges = query_gr,
            rowToRowMat = new_rowToRowMat,
            colToRowMatCols = new_colToRowMatCols,
            colData = sample_metadata,
            name_VAR = name_VAR,
            position_VAR = position_VAR,
            value_VAR = value_VAR,
            region_VAR = region_VAR,
            fetch_config = fetch_config,
            metadata = obj_history)
    }
}

#' ChIPtsne2.from_FetchConfig
#'
#'
#' @param fetch_config An object of class FetchConfig
#' @param query_gr The query GRanges object used to fetch prof_dt.
#' @param region_metadata Metadata to append to rowRanges, mcols of query_gr will also be used.
#' @param obj_history Existing history for object, may be used to describe origin on prof_dt or query_gr.
#' @param init If TRUE, initialize history with birthday, session_info, and chiptsne2_version
#'
#' @return A ChIPtsne2 created using profiles fetched used provided FetchConfig
#' @export
#'
#' @examples
#' bam_cfg_f = exampleBamConfigFile()
#' fetch_config = FetchConfig.load_config(bam_cfg_f)
#' query_gr = seqsetvis::CTCF_in_10a_overlaps_gr
#' ct2 = ChIPtsne2.from_FetchConfig(fetch_config, query_gr)
#' ct2
ChIPtsne2.from_FetchConfig = function(fetch_config,
                                      query_gr,
                                      region_metadata = NULL,
                                      obj_history = list(),
                                      init = TRUE
){
    name_VAR = fetch_config@name_VAR
    sample_metadata = fetch_config@meta_data
    query_gr = seqsetvis::prepare_fetch_GRanges_names(query_gr)
    query_gr = GenomicRanges::resize(query_gr, width = fetch_config@view_size, fix = "center")
    query_gr = seqsetvis::prepare_fetch_GRanges_width(query_gr, win_size = fetch_config$window_size)

    fetch_res = fetch_signal_at_features(fetch_config, query_gr)
    prof_dt = fetch_res$prof_dt

    ct2 = ChIPtsne2.from_tidy(prof_dt = prof_dt,
                              name_VAR = name_VAR,
                              query_gr = query_gr,
                              sample_metadata = sample_metadata,
                              region_metadata = region_metadata,
                              obj_history = obj_history,
                              fetch_config = fetch_config,
                              init = init)
    ct2
}

#' ChIPtsne2.history
#'
#' @param ct2 `r doc_ct2_nrr()`
#'
#' @return list of history items
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' ChIPtsne2.history(ct2)
#' ct2 = setNameVariable(ct2, "name2")
#' # the previous operation is now in the history
#' ChIPtsne2.history(ct2)
ChIPtsne2.history = function(ct2){
    ct2@metadata
}

.cloneChIPtsne2 = function(ct2, new_rowToRowMat = NULL, new_colToRowMatCols = NULL, new_name_VAR = NULL, new_position_VAR = NULL, new_value_VAR = NULL, new_region_VAR = NULL, new_fetch_config = NULL, new_rowRanges = NULL, new_colData = NULL, new_assays = NULL, new_metadata = NULL){
    if(is.null(new_rowToRowMat)) new_rowToRowMat = rowToRowMat(ct2)
    if(is.null(new_colToRowMatCols)) new_colToRowMatCols = colToRowMatCols(ct2)
    if(is.null(new_name_VAR)) new_name_VAR = ct2@name_VAR
    if(is.null(new_position_VAR)) new_position_VAR = ct2@position_VAR
    if(is.null(new_value_VAR)) new_value_VAR = ct2@value_VAR
    if(is.null(new_region_VAR)) new_region_VAR = ct2@region_VAR
    if(is.null(new_fetch_config)) new_fetch_config = ct2@fetch_config
    if(is.null(new_rowRanges)) new_rowRanges = rowRanges(ct2)
    if(is.null(new_colData)) new_colData = colData(ct2)
    if(is.null(new_assays)) new_assays = as.list(ct2@assays@data)
    if(is.null(new_metadata)) new_metadata = ct2@metadata


    ChIPtsne2(rowToRowMat = new_rowToRowMat,
              colToRowMatCols = new_colToRowMatCols,
              name_VAR = new_name_VAR,
              position_VAR = new_position_VAR,
              value_VAR = new_value_VAR,
              region_VAR = new_region_VAR,
              fetch_config = new_fetch_config,
              rowRanges = new_rowRanges,
              colData = new_colData,
              assays = new_assays,
              metadata = new_metadata)
}

.cloneChIPtsne2_no_rowData = function(ct2, new_rowToRowMat = NULL, new_colToRowMatCols = NULL, new_name_VAR = NULL, new_position_VAR = NULL, new_value_VAR = NULL, new_region_VAR = NULL, new_fetch_config = NULL, new_rowData = NULL, new_colData = NULL, new_assays = NULL, new_metadata = NULL){
    if(is.null(new_rowToRowMat)) new_rowToRowMat = rowToRowMat(ct2)
    if(is.null(new_colToRowMatCols)) new_colToRowMatCols = colToRowMatCols(ct2)
    if(is.null(new_name_VAR)) new_name_VAR = ct2@name_VAR
    if(is.null(new_position_VAR)) new_position_VAR = ct2@position_VAR
    if(is.null(new_value_VAR)) new_value_VAR = ct2@value_VAR
    if(is.null(new_region_VAR)) new_region_VAR = ct2@region_VAR
    if(is.null(new_fetch_config)) new_fetch_config = ct2@fetch_config
    if(is.null(new_rowData)) new_rowData = rowData(ct2)
    if(is.null(new_colData)) new_colData = colData(ct2)
    if(is.null(new_assays)) new_assays = as.list(ct2@assays@data)
    if(is.null(new_metadata)) new_metadata = ct2@metadata


    ChIPtsne2_no_rowRanges(rowToRowMat = new_rowToRowMat,
              colToRowMatCols = new_colToRowMatCols,
              name_VAR = new_name_VAR,
              position_VAR = new_position_VAR,
              value_VAR = new_value_VAR,
              region_VAR = new_region_VAR,
              fetch_config = new_fetch_config,
              rowData = new_rowData,
              colData = new_colData,
              assays = new_assays,
              metadata = new_metadata)
}

#' cloneChIPtsne2
#'
#' Create a clone (copy) of a ChIPtsne2 object replacing any specified components.
#'
#' @param ct2 `r doc_ct2()` or `r doc_ct2_nrr()`
#'
#' @param new_rowToRowMat rowToRowMat slot override
#' @param new_colToRowMatCols colToRowMatCols slot override
#' @param new_name_VAR name_VAR slot override
#' @param new_position_VAR position_VAR slot override
#' @param new_value_VAR value_VAR slot override
#' @param new_region_VAR region_VAR slot override
#' @param new_fetch_config slot override
#' @param new_rowRanges rowRanges slot override
#' @param new_colData colData slot override
#' @param new_assays assays slot override
#' @param new_metadata metadata slot override
#'
#' @export
#' @return Clone of input `r doc_ct2()` or `r doc_ct2_nrr()` with specified slots modified.
#' @rdname ct2-clone
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' getNameVariable(ct2)
#' ct2 = cloneChIPtsne2(ct2, new_name_VAR = "new_name")
#' getNameVariable(ct2)
#'
#'
setGeneric("cloneChIPtsne2", function(ct2, new_rowToRowMat = NULL, new_colToRowMatCols = NULL, new_name_VAR = NULL, new_position_VAR = NULL, new_value_VAR = NULL, new_region_VAR = NULL, new_fetch_config = NULL, new_rowRanges = NULL, new_colData = NULL, new_assays = NULL, new_metadata = NULL) standardGeneric("cloneChIPtsne2"))

#' @export
#' @rdname ct2-clone
setMethod("cloneChIPtsne2", c("ChIPtsne2_no_rowRanges"), .cloneChIPtsne2)

.cloneChIPtsne2_fromTidy = function(ct2, new_prof_dt = NULL, new_obj_history = NULL, new_name_VAR = NULL, new_position_VAR = NULL, new_value_VAR = NULL, new_region_VAR = NULL, new_fetch_config = NULL, new_query_gr = NULL, init = FALSE, new_sample_metadata = NULL, new_region_metadata = NULL){
    if(is.null(new_prof_dt)) new_prof_dt = getTidyProfile(ct2)
    if(is.null(new_obj_history)) new_obj_history = ChIPtsne2.history(ct2)
    if(is.null(new_name_VAR)) new_name_VAR = ct2@name_VAR
    if(is.null(new_position_VAR)) new_position_VAR = ct2@position_VAR
    if(is.null(new_value_VAR)) new_value_VAR = ct2@value_VAR
    if(is.null(new_region_VAR)) new_region_VAR = ct2@region_VAR
    if(is.null(new_fetch_config)) new_fetch_config = ct2@fetch_config
    if(is(ct2, "ChIPtsne2")){
        if(is.null(new_query_gr)) new_query_gr = rowRanges(ct2)
    }else{
        new_query_gr = NULL
    }
    if(is.null(new_sample_metadata)) new_sample_metadata = getSampleMetaData(ct2)

    ChIPtsne2.from_tidy(prof_dt = new_prof_dt,
                        query_gr = new_query_gr,
                        sample_metadata = new_sample_metadata,
                        region_metadata = new_region_metadata,
                        name_VAR = new_name_VAR,
                        position_VAR = new_position_VAR,
                        value_VAR = new_value_VAR,
                        region_VAR = new_region_VAR,
                        obj_history = new_obj_history,
                        fetch_config = new_fetch_config,
                        init = init)
}

#' @export
#' @rdname ct2-clone
setGeneric("cloneChIPtsne2_fromTidy", function(ct2, new_prof_dt = NULL, new_obj_history = NULL, new_name_VAR = NULL,
                                               new_position_VAR = NULL, new_value_VAR = NULL, new_region_VAR = NULL,
                                               new_fetch_config = NULL, new_query_gr = NULL, init = FALSE, new_sample_metadata = NULL,
                                               new_region_metadata = NULL) standardGeneric("cloneChIPtsne2_fromTidy"))

#' @export
#' @rdname ct2-clone
setMethod("cloneChIPtsne2_fromTidy", c("ChIPtsne2_no_rowRanges"), .cloneChIPtsne2_fromTidy)
