random_model <- function(x, Co) {
  
  # x=vector of s and connectance
  if (length(x) > 1) {
    S <- x[1]
    C <- x[2]
  } else {
    S <- x
    C <- Co
  }
  
  if (S > 1 & C > 0 & is.finite(S) & is.finite(C)) {
    n.i <- sort(runif(S), decreasing = FALSE)
    r.i <- suppressWarnings(rbeta(S, 1, ((1 / (2 * C)) - 1)) * n.i)
    c.i <- suppressWarnings(runif(S, r.i / 2, n.i))
    
    if (any(is.na(c.i)) == FALSE & any(is.na(r.i)) == FALSE & any(is.na(n.i)) == FALSE) {
      a <- matrix(0, nrow = S, ncol = S)
      for (i in 2:S) {
        for (j in 1:S) {
          if (n.i[j] > (c.i[i] - (0.5 * r.i[i])) & n.i[j] <
              (c.i[i] + 0.5 * r.i[i])) {
            a[j, i] <- 1
          }
        }
      }
      
      rownames(a) <- 1:nrow(a)
      colnames(a) <- 1:ncol(a)
    } else {
      a <- NULL
    }
    
    return(a)
  } else {
    return(NULL)
  }
}