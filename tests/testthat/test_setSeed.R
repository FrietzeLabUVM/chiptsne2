testthat::context("setSeed_and_rerun_history")
library(chiptsne2)
library(testthat)

# some functions depend on random seed
# the same seed should result in identical results from these functions
# just as important, a different seed should result in different results from these functions

if(FALSE){
    #since this toy data is small, it's actually not hard to get identical clustering with differnt seeds
    #lets find a seed that works
    prof_dt = getTidyProfile(ct2)
    set.seed(0)
    clust_0 = seqsetvis::ssvSignalClustering(prof_dt, nclust = 3)
    set.seed(555)
    clust_1 = seqsetvis::ssvSignalClustering(prof_dt, nclust = 3)
    all(clust_0$cluster_id == clust_1$cluster_id)
    clust_0[cluster_id != clust_1$cluster_id]

    all(factor(clust_0$id) == factor(clust_1$id))
}

dimRed_to_test = c(dimReduceTSNE, dimReduceUMAP, dimReducePCA)

for(dimRed in dimRed_to_test){
    ct2 = exampleChIPtsne2.with_meta()
    ct2_1 = ct2 %>% setSeed(0) %>% groupRegionsBySignalCluster() %>% dimReduceUMAP() %>% groupRegionsByDimReduceCluster()

    ct2_re = exampleChIPtsne2.with_meta()
    ct2_re = rerun_history(ct2_re, ct2_1)

    altered_history = ChIPtsne2.history(ct2_1)

    altered_history$setSeed$ARG$seed = 555

    ct2_alt = exampleChIPtsne2.with_meta()
    ct2_alt = rerun_history(ct2_alt, altered_history)



    meta_1 = getRegionMetaData(ct2_1)
    meta_re = getRegionMetaData(ct2_re)
    meta_alt = getRegionMetaData(ct2_alt)


    expect_equal(colnames(meta_1), colnames(meta_re))
    expect_equal(meta_1, meta_re)

    expect_equal(colnames(meta_1), colnames(meta_alt))
    should_match = seq(1, 4)
    expect_equal(meta_1[, should_match], meta_alt[, should_match])
    should_not_match = seq(5, 8)

    test_res = logical()
    for(i in should_not_match){
        message(i)
        is_equal = all(as.character(meta_1[, i]) == as.character(meta_alt[, i]))
        test_res[[colnames(meta_1)[i]]] = is_equal
    }
    if(any(test_res)){
        message("the following should not match but do: ", paste(names(test_res)[test_res], collapse = ", "))
    }

    expect_true(all(!test_res))

}
