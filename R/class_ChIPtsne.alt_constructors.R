ChIPtsne2.from_config = function(){

}

#' ChIPtsne2.from_tidy
#'
#' @param prof_dt Profile data.table, as returned from seqsetvis::ssvFetch* functions
#' @param query_gr The query GRanges object used to fetch prof_dt.
#' @param meta_dt Additional information for each name_VAR entry in prof_dt.
#' @param name_VAR Variable name that contains sample ids/names. Links prof_dt to meta_dt. Default is "sample".
#' @param position_VAR Variable name that contains positional information in prof_dt. Default is "x".
#' @param value_VAR Variable name that contains signal value information in prof_dt. Default is "y".
#' @param region_VAR Variable name that contains region ID information in prof_dt. Default is "id".
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
                               auto_sample_metadata = TRUE){
    #basic VAR checks
    if(!all(c(name_VAR, position_VAR, value_VAR, region_VAR) %in% colnames(prof_dt))){
        missed = !c(name_VAR, position_VAR, value_VAR, region_VAR) %in% colnames(prof_dt)
        missed_msg = paste(c("name_VAR", "position_VAR", "value_VAR", "region_VAR")[missed],
                           c(name_VAR, position_VAR, value_VAR, region_VAR)[missed], sep = ": ")
        stop(paste(c("Missing required VAR in prof_dt:", missed_msg, collapse = "\n")))
    }
    if(is.null(names(query_gr))){
        stop("names() must be set on query_gr. Maybe call seqsetvis::prepare_fetch_GRanges_names on query_gr before fetching prof_dt?")
    }
    if(!all(unique(as.character(prof_dt[[region_VAR]])) %in% names(query_gr))){
        stop("region_VAR: ", region_VAR, " in prof_dt is not consistent with names of query_gr")
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
            cn = levels(sample_metadata[[name_VAR]])
        }else if(is.character(sample_metadata[[name_VAR]])){
            cn = unique(sample_metadata[[name_VAR]])
        }
    }else{
        if(is.factor(prof_dt[[name_VAR]])){
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
            rn = levels(prof_dt[[region_VAR]])
        }else if(is.character(prof_dt[[region_VAR]])){
            rn = unique(prof_dt[[region_VAR]])
        }
    }
    if(is.null(rn)){
        stop("Could not determine row order from prof_dt from region_VAR: ", region_VAR, "\nIs ", region_VAR, " present and a character or factor?")
    }
    if(!setequal(names(query_gr), rn)){
        stop("names(query_gr) is not consistent with region_VAR: ", region_VAR, " in prof_dt")
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
    prof_mat = prof_mat[names(query_gr),]

    #create max assay
    prof_max = prof_dt %>%
        dplyr::group_by(id, sample) %>%
        dplyr::summarise(y = max(y)) %>%
        tidyr::pivot_wider(
            names_from = all_of(c(name_VAR)),
            id_cols = all_of(c(region_VAR)),
            values_from = all_of(c(value_VAR))
        )
    prof_max_mat = as.matrix(prof_max[, -1])
    rownames(prof_max_mat) = prof_max[[region_VAR]]
    prof_max_mat = prof_max_mat[names(query_gr),]

    xy_dt = tsne_from_profile_mat(prof_mat)

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
                          "cluster_id")
            sample_metadata = prof_dt %>%
                # dplyr::select(all_of(c(name_VAR))) %>%
                dplyr::select(!any_of(c(drop_vars))) %>%
                unique
            if(nrow(sample_metadata) != length(cn)){
                stop("Something has gone wrong attempting to automatically derive sample_metadata. Either supply explicitly or disable with auto_sample_metadata = FALSE")
            }
        }else{
            sample_metadata = data.frame(V1 = cn)
            colnames(sample_metadata) = name_VAR
        }
    }
    sample_metadata = as.data.frame(sample_metadata)
    rownames(sample_metadata) = sample_metadata[[name_VAR]]
    sample_metadata[[name_VAR]] = NULL

    if(!is.null(region_metadata)){
        #merge region_metadata into query_gr
        #region_metadata is not used after
        region_metadata
        new_mcols = cbind(
            mcols(query_gr[region_metadata[[region_VAR]]]),
            as.data.frame(region_metadata %>% dplyr::select(!dplyr::all_of(c(region_VAR))))
        )
        mcols(query_gr) = NULL
        mcols(query_gr) = new_mcols
    }


    map_dt = prof_dt %>%
        dplyr::select(all_of(c(name_VAR, position_VAR))) %>%
        unique
    map_dt = dplyr::mutate(map_dt, cn = paste(get(name_VAR), get(position_VAR), sep = "_"))
    stopifnot(map_dt$cn == colnames(prof_mat))
    map_dt = dplyr::mutate(map_dt, nr = seq(nrow(map_dt)))
    map_list = split(map_dt$n, map_dt[[name_VAR]])

    #impose cn and rn
    prof_max_mat = prof_max_mat[rn, cn]
    query_gr = query_gr[rn]
    prof_mat = prof_mat[rn, ]
    map_list = map_list[cn]

    ChIPtsne2(assay = list(max = prof_max_mat),
              rowRanges = query_gr,
              rowToRowMat = prof_mat,
              colToRowMatCols = map_list,
              colData = sample_metadata,
              metadata = list(time = date()))
}

tsne_from_profile_mat = function(prof_mat){
    tsne_res = Rtsne::Rtsne(prof_mat)
    xy_dt = as.data.frame(tsne_res$Y)
    colnames(xy_dt) =  c("tx", "ty")
    rownames(xy_dt) = rownames(prof_mat)
    xy_dt
}
