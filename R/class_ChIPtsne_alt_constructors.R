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
#' @param init If TRUE, initialize history with birthday, session_info, and chiptsne2_version
#'
#' @return ChIPtsne2 object with tsne results
#' @export
#'
#' @examples
#' query_gr = CTCF_in_10a_overlaps_gr
#' prof_dt = CTCF_in_10a_profiles_dt
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
        init_history = list(birthday = date(), session_info = sessionInfo(), chiptsne2_version = utils::packageDescription("chiptsne2")$Version)
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
    if(!is.null(query_gr)){
        if(is.null(names(query_gr))){
            stop("names() must be set on query_gr. Maybe call seqsetvis::prepare_fetch_GRanges_names on query_gr before fetching prof_dt?")
        }
        if(!all(unique(as.character(prof_dt[[region_VAR]])) %in% names(query_gr))){
            stop("region_VAR: ", region_VAR, " in prof_dt is not consistent with names of query_gr")
        }
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
    if(!is.null(region_metadata)){
        if(!region_VAR %in% colnames(region_metadata)){
            stop("region_VAR: ", region_VAR, " not found in region_metadata")
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
    if(!is.null(query_gr)){
        if(!setequal(names(query_gr), rn)){
            stop("names(query_gr) is not consistent with region_VAR: ", region_VAR, " in prof_dt")
        }
    }

    #create wide profile matrix
    tmp_wide = tidyr::pivot_wider(
        prof_dt,
        names_from = all_of(c(name_VAR, position_VAR)),
        values_from = all_of(c(value_VAR)),
        id_cols = all_of(c(region_VAR))
    )
    prof_mat = as.matrix(tmp_wide[, -1])
    rownames(prof_mat) = tmp_wide[[region_VAR]]
    prof_mat = prof_mat[rn,]

    tmp = c("value_VAR")
    names(tmp) = value_VAR
    #create max assay
    #unsure how to programmatically use summarize the way I want
    prof_max = prof_dt %>%
        dplyr::group_by(.data[[region_VAR]], .data[[name_VAR]]) %>%
        dplyr::summarise(value_VAR = max(.data[[value_VAR]])) %>%
        dplyr::rename(all_of(tmp)) %>%
        tidyr::pivot_wider(
            names_from = all_of(c(name_VAR)),
            id_cols = all_of(c(region_VAR)),
            values_from = all_of(c(value_VAR))
        )
    prof_max_mat = as.matrix(prof_max[, -1])
    rownames(prof_max_mat) = prof_max[[region_VAR]]
    prof_max_mat = prof_max_mat[rn,]

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
                dplyr::select(!any_of(c(drop_vars))) %>%
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

    if(!is.null(query_gr)){
        if(!is.null(region_metadata)){
            #merge region_metadata into query_gr
            #region_metadata is not used after
            query_gr = .add_region_metadata(query_gr, region_metadata, region_VAR, overwrite = TRUE)
            query_gr = query_gr[rn]
        }
    }

    map_dt = prof_dt %>%
        dplyr::select(all_of(c(name_VAR, position_VAR))) %>%
        unique
    map_dt = dplyr::mutate(map_dt, cn = paste(get(name_VAR), get(position_VAR), sep = "_"))
    stopifnot(map_dt$cn == colnames(prof_mat))
    map_dt = dplyr::mutate(map_dt, nr = seq(nrow(map_dt)))
    map_list = split(map_dt$cn, map_dt[[name_VAR]])

    #impose cn and rn
    prof_max_mat = prof_max_mat[rn, cn, drop = FALSE]
    prof_mat = prof_mat[rn, , drop = FALSE]
    map_list = map_list[cn]
    sample_metadata = sample_metadata[cn, , drop = FALSE]

    if(is.null(query_gr)){
        ChIPtsne2_no_rowRanges(
            assay = list(max = prof_max_mat),
            rowToRowMat = prof_mat,
            colToRowMatCols = map_list,
            colData = sample_metadata,
            name_VAR = name_VAR,
            position_VAR = position_VAR,
            value_VAR = value_VAR,
            region_VAR = region_VAR,
            fetch_config = fetch_config,
            metadata = obj_history)
    }else{
        ChIPtsne2(
            assay = list(max = prof_max_mat),
            rowRanges = query_gr,
            rowToRowMat = prof_mat,
            colToRowMatCols = map_list,
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
#' bam_cfg_f = system.file("extdata/bam_config.csv", package = "chiptsne2", mustWork = TRUE)
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
    # name_lev = NULL
    # if(is.factor(sample_metadata[[name_VAR]])){
    #     name_lev = levels(sample_metadata[[name_VAR]])
    # }else if(is.character(sample_metadata[[name_VAR]])){
    #     name_lev = unique(sample_metadata[[name_VAR]])
    # }
    # if(!is.null(name_lev)){
    #     prof_dt[[name_VAR]] = factor(prof_dt[[name_VAR]], levels = name_lev)
    # }
    #
    # prof_dt[[name_VAR]] = droplevels(prof_dt[[name_VAR]])
    # sample_metadata[[name_VAR]] = droplevels(sample_metadata[[name_VAR]])

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
#' @param ct2
#'
#' @return list of history items
#' @export
ChIPtsne2.history = function(ct2){
    ct2@metadata
}

.cloneChIPtsne2 = function(ct2, .rowToRowMat = NULL, .colToRowMatCols = NULL, .name_VAR = NULL, .position_VAR = NULL, .value_VAR = NULL, .region_VAR = NULL, .fetch_config = NULL, .rowRanges = NULL, .colData = NULL, .assays = NULL, .metadata = NULL){
    if(is.null(.rowToRowMat)) .rowToRowMat = rowToRowMat(ct2)
    if(is.null(.colToRowMatCols)) .colToRowMatCols = colToRowMatCols(ct2)
    if(is.null(.name_VAR)) .name_VAR = ct2@name_VAR
    if(is.null(.position_VAR)) .position_VAR = ct2@position_VAR
    if(is.null(.value_VAR)) .value_VAR = ct2@value_VAR
    if(is.null(.region_VAR)) .region_VAR = ct2@region_VAR
    if(is.null(.fetch_config)) .fetch_config = ct2@fetch_config
    if(is.null(.rowRanges)) .rowRanges = rowRanges(ct2)
    if(is.null(.colData)) .colData = colData(ct2)
    if(is.null(.assays)) .assays = as.list(ct2@assays@data)
    if(is.null(.metadata)) .metadata = ct2@metadata


    ChIPtsne2(rowToRowMat = .rowToRowMat,
              colToRowMatCols = .colToRowMatCols,
              name_VAR = .name_VAR,
              position_VAR = .position_VAR,
              value_VAR = .value_VAR,
              region_VAR = .region_VAR,
              fetch_config = .fetch_config,
              rowRanges = .rowRanges,
              colData = .colData,
              assays = .assays,
              metadata = .metadata)
}

#' @export
setGeneric("cloneChIPtsne2", function(ct2, .rowToRowMat = NULL, .colToRowMatCols = NULL, .name_VAR = NULL, .position_VAR = NULL, .value_VAR = NULL, .region_VAR = NULL, .fetch_config = NULL, .rowRanges = NULL, .colData = NULL, .assays = NULL, .metadata = NULL) standardGeneric("cloneChIPtsne2"))

#' @export
setMethod("cloneChIPtsne2", c("ChIPtsne2"), .cloneChIPtsne2)

.cloneChIPtsne2_fromTidy = function(ct2, prof_dt = NULL, obj_history = NULL, name_VAR = NULL, position_VAR = NULL, value_VAR = NULL, region_VAR = NULL, fetch_config = NULL, query_gr = NULL, init = FALSE, sample_metadata = NULL, region_metadata = NULL){
    if(is.null(prof_dt)) prof_dt = getTidyProfile(ct2)
    if(is.null(obj_history)) obj_history = ChIPtsne2.history(ct2)
    if(is.null(name_VAR)) name_VAR = ct2@name_VAR
    if(is.null(position_VAR)) position_VAR = ct2@position_VAR
    if(is.null(value_VAR)) value_VAR = ct2@value_VAR
    if(is.null(region_VAR)) region_VAR = ct2@region_VAR
    if(is.null(fetch_config)) fetch_config = ct2@fetch_config
    if(is.null(query_gr)) query_gr = rowRanges(ct2)
    if(is.null(sample_metadata)) sample_metadata = getSampleMetaData(ct2)

    ChIPtsne2.from_tidy(prof_dt = prof_dt,
                        query_gr = query_gr,
                        sample_metadata = sample_metadata,
                        region_metadata = region_metadata,
                        name_VAR = name_VAR,
                        position_VAR = position_VAR,
                        value_VAR = value_VAR,
                        region_VAR = region_VAR,
                        obj_history = obj_history,
                        fetch_config = fetch_config,
                        init = init)
}

#' @export
setGeneric("cloneChIPtsne2_fromTidy", function(ct2, prof_dt = NULL, obj_history = NULL, name_VAR = NULL, position_VAR = NULL, value_VAR = NULL, region_VAR = NULL, fetch_config = NULL, query_gr = NULL, init = FALSE, sample_metadata = NULL, region_metadata = NULL) standardGeneric("cloneChIPtsne2_fromTidy"))

#' @export
setMethod("cloneChIPtsne2_fromTidy", c("ChIPtsne2"), .cloneChIPtsne2_fromTidy)
