split.bseq <- function(object){
  listname <- "cond"
  reStruc <- object$modelInfo$reStruc[[paste0(listname, "ReStruc")]] ## random-effects structure
  nc <- vapply(reStruc, function(x) x$blockSize, numeric(1)) ## number of RE params per block
  nb <- vapply(reStruc, function(x) x$blockReps, numeric(1)) ## number of blocks per RE (may != nlevs in some cases)
  ### splitting the b's into their respective random effects
  nbseq <- rep.int(seq_along(nb), nb * nc)       ## splitting vector
  return(nbseq)
}

extract.rr.fl.b <- function(object){
  listname <- "cond"
  cnms <- object$modelInfo$reTrms[[listname]]$cnms   ## list of (named) terms and X columns
  flist <- object$modelInfo$reTrms[[listname]]$flist ## list of grouping variables
  levs <- lapply(flist, levels)
  reStruc <- object$modelInfo$reStruc[[paste0(listname, "ReStruc")]] ## random-effects structure
  nc <- vapply(reStruc, function(x) x$blockSize, numeric(1)) ## number of RE params per block

  pl <- object$obj$env$parList(object$fit$par, object$fit$parfull)
  nbseq <- split.bseq(object)      ## splitting vector
  ml.b <- split(pl$b, nbseq)
  ml <- ml.b

  for (i in seq_along(ml.b)) {
    ml[[i]] <- matrix(ml.b[[i]], ncol = nc[i], byrow = TRUE,
                      dimnames = list(levs[[i]], cnms[[i]]))
  }

  get_rank <- function(x){
    if(x[["blockCode"]]==9){
      p <- x$blockSize
      nt <- x$blockNumTheta
      rank <- (2*p + 1 - sqrt((2*p+1)^2 - 8*nt))/2
    } else
      rank <- 0
    return(rank)
  }

  rank <- vapply(object$modelInfo$reStruc$condReStruc,
                 get_rank,
                 FUN.VALUE=numeric(1))
  nlv <- rank[rank > 0]
  rrName <- names(nlv)
  rrBlock <- which(rank > 0)
  b = ml[[rrBlock]][,1:nlv]
  fact_load <- object$obj$env$report(object$fit$parfull)$fact_load[[rrBlock]]

  return(list(fl = fact_load, b = b))
}

#For PIRLS
load.file.type <- function(file.type = "acg.*"){
  files <- paste0("Data/PIRLS/", list.files("Data/PIRLS/", pattern = file.type))
  data <- lapply(files, function(x) read.spss(x, use.value.label = TRUE, to.data.frame = TRUE))
}



theme_mine <- function(base_size = 12) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      strip.background = element_blank(),
      strip.text.x = element_text(size = 12),
      strip.text.y = element_text(size = 12),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12,hjust=1),
      axis.ticks =  element_line(colour = "black"),
      axis.title.x = element_text(size=12),
      axis.title.y = element_text(size=12, angle=90),
      panel.background = element_blank(),
      panel.border =element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1.0, "lines"),
      plot.background = element_blank(),
      # plot.spacing = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
      axis.line.x = element_line(color="black", size = 1),
      axis.line.y = element_line(color="black", size = 1)
    )
}


