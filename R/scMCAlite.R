
#' Mouse Cell Atlas mapping
#'
#' simplized version of origianl scMCA function.
#'
#' @param scdata matrix-like data input, log-normalized data recommended.
#' @param num the number of the most similar cell types returned, default 3.
#'
#' @return a list containing cors_matrix, top_cors, scMCA and scMCA_similarity.
#'
#' @importFrom stats cor
#' @importFrom reshape2 melt
#' @import dplyr
#' @import Matrix
#' @export
scMCAlite <- function(scdata, num = 3){
  
  tst.expr <- scdata
  used.gene <- intersect(rownames(ref.expr), rownames(tst.expr))
  tst.expr <- tst.expr[used.gene, ]
  
  message(length(used.gene), " genes used, calculating correlation coefficient...")
  cors <- cor(log2(ref.expr[used.gene, ] + 1), as.matrix(tst.expr), method = "pearson")
  
  cors.index <- apply(cors, 2, function(x){order(x, decreasing = T)[1:num]})
  cors.index <- as.integer(cors.index) %>% unique() %>% sort()
  
  type.res <- apply(cors, 2, function(x){rownames(cors)[which.max(x)]})
  
  message(length(cors.index), " types expected, extracting top ", num, " results...")
  cors_in = cors[cors.index, ]
  cors_out = melt(cors_in)[c(2, 1, 3)]
  colnames(cors_out)<- c("Cell", "Type", "Score")
  cors_out <- as.data.frame(cors_out %>% group_by(Cell) %>% top_n(n = num, wt = Score))
  
  result <- list()
  result[["cors_matrix"]] <- cors
  result[['top_cors']] <- num
  result[['scMCA']] <- type.res
  result[['scMCA_similarity']] <- cors_out
  return(result)
  
}
