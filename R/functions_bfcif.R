#' bfcif
#'
#' Conditionally runs a function if it's results don't exist in the cache.  Cache is controlled by BiocFileCache.
#'
#' @param FUN A function that takes zero arguments.  This function can be a wrapper around other functions; ie. FUN = function(){mean(x)}.
#' @param rname The unique identifier for the results in the cache. The recommendation is to use either a unique and meaningful description or digest::digest() on a list containing FUN and it's parameters.
#' @param bfc A BiocFileCache object, typically from BiocFileCache::BiocFileCache().
#' @param version A version indicator string to further distinguish cache entries.  Typically, you want to iterate this value when code outside of FUN or input parameters have changed and you want to force FUN to reevaluate. Default is SQC_CACHE_VERSION option or v3 if option is not set.
#' @param force_overwrite If TRUE, FUN will be rerun regardless of cache state. If it exists, current cache contents will be overwritten. Default is SQC_FORCE_CACHE_OVERWRITE option or FALSE if option is not set.
#' @param return_path_only If TRUE, FUN will not be run and instead the cache path is returned. Default is FALSE.
#' @param verbose If TRUE, status is reported via messages. Default is value of SQC_CACHE_VERBOSE option or FALSE if option not set.
#'
#' @return Result of FUN, from cache if available.
#' @export
#' @importFrom BiocFileCache bfcquery bfcnew bfcrpath BiocFileCache
#'
#' @examples
#' bfc = BiocFileCache::BiocFileCache()
#' test_fun = function(x)mean(seq(10))
#' bfcif(test_fun, "test1", bfc, verbose = TRUE)
#' bfcif(test_fun, "test1", bfc, verbose = TRUE,
#'   return_path_only = TRUE)
#' bfcif(test_fun, "test1", bfc, verbose = TRUE,
#'   force_overwrite = TRUE)
bfcif = function(FUN, rname, bfc = NULL, version = NULL, force_overwrite = FALSE, return_path_only = FALSE, verbose = TRUE){
    if(is.null(bfc)){
        bfc = BiocFileCache::BiocFileCache()
    }
    if(!is.null(version)){
        rname = paste0(rname, "_", version)
    }
    # is rname in cache?
    if(nrow(BiocFileCache::bfcquery(bfc, query = rname, field = "rname", exact = TRUE)) == 0){
        if(verbose) message("results not in cache. ", appendLF = FALSE)
        cache_path = BiocFileCache::bfcnew(bfc, rname = rname)

    }else{
        if(verbose) message("previous cache results found. ", appendLF = FALSE)
        cache_path = BiocFileCache::bfcrpath(bfc, rname)
    }
    if(return_path_only){
        if(verbose) message("returning cache path.")
        return(cache_path)
    }
    # does cached file exist?
    if(file.exists(cache_path) && !force_overwrite){
        if(verbose) message("loading previous cache results...")
        load(BiocFileCache::bfcrpath(bfc, rname))
    }else{
        if(verbose) message("running function...", appendLF = FALSE)
        res = FUN()
        if(verbose) message("caching results...")
        save(res, file = cache_path)
    }
    # return either new results or cached results
    res
}

#' get_args
#'
#' @param env environment to capture variables from. Default of `parent.frame()` retrieves the call function environment.
#' @param to_ignore character vector of variable names to omit. Default of `character()` is nn variables.
#' @param ... named extra variabes to include.
#'
#' @return named list or variables and their values in `env`
#' @export
#'
#' @examples
#' test_fun = function(x = 10, ...){
#'     args = get_args(...)
#'     y = mean(seq(x))
#'     args
#' }
#' test_fun(5)
#' # example with ...
#' test_fun(x = 6, y = 7, z = 8)
#' # example to_ignore
#' test_fun(x = 6, y = 7, z = 8, to_ignore = "x")
get_args = function(env = parent.frame(), to_ignore = character(), ...){
    if(is.function(env)) env = env()
    stopifnot(is.environment(env))
    args = c(as.list(env), list(...))
    args = args[!names(args) %in% to_ignore]
    args[order(names(args))]
}

#' digest_args
#'
#' returns digest::digest results of named list of arguments to calling function
#'
#' @inheritParams get_args
#'
#' @return A md5 digest of all variabes in calling environment. For a function
#'   is the same as calling arguments.
#' @importFrom digest digest
#' @export
#'
#'
#' @examples
#' #The most common usage is to simply collect all local variables in a function
#' test_fun = function(x = 1, y = 2){
#'   digest_args()
#' }
#' test_fun()
#'
#' #Specified variables may be ignored
#' test_fun2 = function(x = 1, y = 2){
#'   digest_args(to_ignore = "x")
#' }
#' test_fun2()
#'
#' #Additional variables can also be added from higher environments
#' global_z = 3
#' test_fun3 = function(x = 1, y = 2){
#'   digest_args(env = parent.frame, to_ignore = character(), z = global_z)
#' }
#' test_fun3()
digest_args = function(env = parent.frame(), to_ignore = character(), ...){
    digest::digest(get_args(env, to_ignore, ...))
}
