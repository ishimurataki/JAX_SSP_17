#########################################################################
# ASSOCIATION PLOTTING WITH DIFFERENT SNP FUNCTIONS AS DIFFERENT COLORS #
#########################################################################

#### SETTING UP ALL THE FUNCTIONS ######
align_scan1_map <-
  function(scan1_output, map)
  {
    if(!is.list(map)) stop("map should be a list")
    
    scan1_names <- rownames(scan1_output)
    map_names <- map2markernames(map)
    map_chr <- map2chr(map)
    
    # perfectly fine
    if(length(scan1_names) == length(map_names) &&
       all(scan1_names == map_names))
      return(list(scan1_output=scan1_output, map=map))
    
    if(!any(scan1_names %in% map_names))
      stop("scan1 output and map have no markers in common")
    
    # subset scan1_output to markers in the map
    #    and use order as in map
    scan1_attr <- attributes(scan1_output)
    names_ordered <- map_names[map_names %in% scan1_names]
    scan1_output <- scan1_output[names_ordered,,drop=FALSE]
    for(a in c("SE", "sample_size", "hsq", "class"))
      attr(scan1_output, a) <- scan1_attr[[a]]
    scan1_names <- rownames(scan1_output)
    
    # subset map to markers in scan1_output
    if(any(!(map_names %in% scan1_names))) {
      map_attr <- attributes(map)
      keep <- map_names %in% scan1_names
      uchr <- unique(map_chr[keep])
      map <- map[uchr]
      if("is_x_chr" %in% names(map_attr))
        map_attr$is_x_chr <- map_attr$is_x_chr[uchr]
      for(i in seq_along(map))
        map[[i]] <- map[[i]][names(map[[i]]) %in% scan1_names]
      for(a in c("is_x_chr", "class"))
        attr(map, a) <- map_attr[[a]]
    }
    
    list(scan1_output=scan1_output, map=map)
  }

# grab marker names as a vector
map2markernames <-
  function(map)
  {
    nam <- unlist(lapply(map, names))
    names(nam) <- NULL
    nam
  }

# grab chromosome IDs as a vector
map2chr <-
  function(map)
  {
    chr <- rep(names(map), vapply(map, length, 0))
    names(chr) <- map2markernames(map)
    chr
  }

plot_scan1_experimental <-
  function(x, map, lodcolumn=1, chr=NULL, add=FALSE, gap=25,
           bgcolor="gray90", altbgcolor="gray85", ...)
  {
    if(is.null(map)) stop("map is NULL")
    
    if(!is.list(map)) map <- list(" "=map) # if a vector, treat it as a list with no names
    
    # subset chromosomes
    if(!is.null(chr)) {
      chri <- match(chr, names(map))
      if(any(is.na(chri)))
        stop("Chromosomes ", paste(chr[is.na(chri)], collapse=", "), " not found")
      x <- qtl2scan::subset_scan1(x, map, chr)
      map <- map[chri]
    }
    
    # align scan1 output and map
    tmp <- align_scan1_map(x, map)
    x <- tmp$scan1
    map <- tmp$map
    
    if(nrow(x) != length(unlist(map)))
      stop("nrow(x) [", nrow(x), "] != number of positions in map [",
           length(unlist(map)), "]")
    
    # pull out lod scores
    if(length(lodcolumn) > 1) { # If length > 1, take first value
      warning("lodcolumn should have length 1; one first element used.")
      lodcolumn <- lodcolumn[1]
    }
    if(is.character(lodcolumn)) { # turn column name into integer
      tmp <- match(lodcolumn, colnames(x))
      if(is.na(tmp))
        stop('lodcolumn "', lodcolumn, '" not found')
      lodcolumn <- tmp
    }
    if(lodcolumn < 1 || lodcolumn > ncol(x))
      stop("lodcolumn [", lodcolumn, "] out of range (should be in 1, ..., ", ncol(x), ")")
    lod <<- unclass(x)[,lodcolumn]
    
    # internal function; trick to be able to pull things out of "..."
    #    but still have some defaults for them
    plot_scan1_internal <-
      function(map, lod, add=FALSE, gap,
               bgcolor, altbgcolor,
               lwd=2, col="green", xlab=NULL, ylab="LOD score",
               xlim=NULL, ylim=NULL, xaxs="i", yaxs="i",
               main="", mgp.x=c(2.6, 0.5, 0), mgp.y=c(2.6, 0.5, 0),
               mgp=NULL, las=1,
               hlines=NULL, hlines.col="white", hlines.lwd=1, hlines.lty=1,
               vlines=NULL, vlines.col="white", vlines.lwd=1, vlines.lty=1,
               ...)
      {
        dots <- list(...)
        onechr <- (length(map)==1) # single chromosome
        
        xpos <- map_to_xpos(map, gap)
        chrbound <- map_to_boundaries(map, gap)
        
        if(!add) { # new plot
          if(is.null(ylim))
            ylim <- c(0, max(lod, na.rm=TRUE)*1.02)
          
          if(is.null(xlim)) {
            xlim <- range(xpos, na.rm=TRUE)
            if(!onechr) xlim <- xlim + c(-gap/2, gap/2)
          }
          
          if(is.null(xlab)) {
            if(onechr) {
              if(names(map) == " ") xlab <- "Position"
              else xlab <- paste("Chr", names(map), "position")
            }
            else xlab <- "Chromosome"
          }
          
          # margin parameters
          if(!is.null(mgp)) mgp.x <- mgp.y <- mgp
          
          # make basic plot
          plot(xpos, lod, xlab="", ylab="", xlim=xlim, ylim=ylim,
               xaxs=xaxs, yaxs=yaxs, xaxt="n", yaxt="n", type="n",
               main=main)
          
          # add background rectangles
          u <- par("usr")
          if(!is.null(bgcolor))
            rect(u[1], u[3], u[2], u[4], col=bgcolor, border=NA)
          if(!is.null(altbgcolor) && !onechr) {
            for(i in seq(2, ncol(chrbound), by=2))
              rect(chrbound[1,i], u[3], chrbound[2,i], u[4], col=altbgcolor, border=NA)
          }
          
          # include axis labels?
          if(is.null(dots$xaxt)) dots$xaxt <- par("xaxt")
          if(is.null(dots$yaxt)) dots$yaxt <- par("yaxt")
          
          # add x axis unless par(xaxt="n")
          if(dots$xaxt != "n") {
            if(onechr) {
              axis(side=1, at=pretty(xlim), mgp=mgp.x, las=las, tick=FALSE)
            }
            else {
              loc <- colMeans(chrbound)
              odd <- seq(1, length(map), by=2)
              even <- seq(2, length(map), by=2)
              axis(side=1, at=loc[odd], names(map)[odd],
                   mgp=mgp.x, las=las, tick=FALSE)
              axis(side=1, at=loc[even], names(map)[even],
                   mgp=mgp.x, las=las, tick=FALSE)
            }
          }
          
          # add y axis unless par(yaxt="n")
          if(dots$yaxt != "n") {
            axis(side=2, at=pretty(ylim), mgp=mgp.y, las=las, tick=FALSE)
          }
          dots$xaxt <- dots$yaxt <- NULL # delete those
          
          # x and y axis labels
          title(xlab=xlab, mgp=mgp.x)
          title(ylab=ylab, mgp=mgp.y)
          
          # grid lines
          if(onechr && !(length(vlines)==1 && is.na(vlines))) { # if vlines==NA (or mult chr), skip lines
            if(is.null(vlines)) vlines <- pretty(xlim)
            abline(v=vlines, col=vlines.col, lwd=vlines.lwd, lty=vlines.lty)
          }
          if(!(length(hlines)==1 && is.na(hlines))) { # if hlines==NA, skip lines
            if(is.null(hlines)) hlines <- pretty(ylim)
            abline(h=hlines, col=hlines.col, lwd=hlines.lwd, lty=hlines.lty)
          }
        }
        
        # plot each chromosome
        indexes <- map_to_index(map)
        for(i in seq(along=indexes))
          lines(xpos[indexes[[i]]], lod[indexes[[i]]],
                lwd=lwd, col=col, ...)
        
        # add box just in case
        box()
      }
    
    # make the plot
    plot_scan1_internal(map=map, lod=lod, add=add, gap=gap,
                        bgcolor=bgcolor, altbgcolor=altbgcolor,
                        ...)
  }


# convert map to x-axis positions for plot_scan1
map_to_xpos <-
  function(map, gap)
  {
    if(length(map)==1) return(map[[1]])
    
    chr_range <- vapply(map, range, c(0,1), na.rm=TRUE)
    
    result <- map[[1]]-chr_range[1,1] + gap/2
    for(i in 2:length(map)) {
      result <- c(result,
                  map[[i]] - chr_range[1,i] + gap + max(result, na.rm=TRUE))
    }
    result
  }

# boundaries of chromosomes in plot_scan1
# first row: left edges
# second row: right edges
map_to_boundaries <-
  function(map, gap)
  {
    if(length(map)==1)
      return(cbind(range(map[[1]], na.rm=TRUE)))
    
    # range of each chromosome
    chr_range <- lapply(map, range, na.rm=TRUE)
    
    # corresponding xpos, as matrix with two rows
    startend <- matrix(map_to_xpos(chr_range, gap), nrow=2)
    
    startend[1,] <- startend[1,] - gap/2
    startend[2,] <- startend[2,] + gap/2
    
    startend
  }

# convert map to list of indexes to LOD vector
map_to_index <-
  function(map)
  {
    if(length(map)==1) {
      map[[1]] <- seq(along=map[[1]])
      return(map)
    }
    
    lengths <- vapply(map, length, 0)
    split(1:sum(lengths), rep(seq(along=map), lengths))
  }


plot_snpasso_experimental <-
  function(scan1output, snpinfo, show_all_snps=TRUE, drop.hilit=NA,
           col.hilit="green", col="pink",
           pch=16, cex=0.75, ylim=NULL, add=FALSE, gap=25,
           bgcolor="gray90", altbgcolor="gray85", ...)
  {
    uindex <- unique(snpinfo$index)
    if(length(uindex) != nrow(scan1output))
      stop("Something is wrong with snpinfo$index.\n",
           "      length(unique(snpinfo$index)) [",
           length(unique(snpinfo$index)), "] != nrow(scan1output) [",
           nrow(scan1output), "].")
    
    if(any(snpinfo$index[uindex] != uindex))
      stop("Something is wrong with snpinfo$index.\n",
           "      snpinfo$index[u] should == u for values in snpinfo$index")
    
    map <- snpinfo_to_map(snpinfo)
    
    if(show_all_snps) {
      tmp <- expand_snp_results(scan1output, map, snpinfo)
      scan1output <- tmp$lod
      map <- tmp$map
    }
    
    # maximum LOD
    maxlod <- max(unclass(scan1output)[,1], na.rm=TRUE)
    
    if(is.null(ylim))
      ylim <- c(0, maxlod*1.02)
    
    if(!is.na(drop.hilit) && !is.null(drop.hilit)){
      colors <- vector("character")
      for(i in 1:length(scan1output)){
        snpid <- rownames(scan1output)[i]
        snpfunction <- snpinfo$csq[snpinfo$snp == snpid]
        if(snpfunction == "3_prime_UTR_variant"){
          col = "orange"
        }
        if(snpfunction == "3_prime_UTR_variant,NMD_transcript_variant"){
          col = "orange"
        }
        if(snpfunction == "5_prime_UTR_variant"){
          col = "purple"
        }
        if(snpfunction == "downstream_gene_variant"){
          col = "green"
        }
        if(snpfunction == "intergenic_variant"){
          col = "darkslateblue"
        }
        if(snpfunction == "intron_variant"){
          col = "yellow"
        }
        if(snpfunction == "intron_variant,NMD_transcript_variant"){
          col = "yellow"
        }
        if(snpfunction == "intron_variant,non_coding_transcript_variant"){
          col = "yellow"
        }
        if(snpfunction == "missense_variant"){
          col = "red"
        }
        if(snpfunction == "non_coding_transcript_exon_variant,non_coding_transcript_variant"){
          col = "brown"
        }
        if(snpfunction == "synonymous_variant"){
          col = "orange"
        }
        if(snpfunction == "upstream_gene_variant"){
          col = "black"
        }
        colors <- c(colors, col)
      }
    }
    col <- colors
    plot_scan1_experimental(scan1output, map, lodcolumn=1, bgcolor=bgcolor, altbgcolor=altbgcolor, ylim=ylim,
               gap=gap, add=add, col = col, type="p", cex=cex, pch=pch, ...)
  }


# expand snp association results according to snpinfo
expand_snp_results <-
  function(snp_results, map, snpinfo)
  {
    snpinfo <- split(snpinfo, snpinfo$chr)
    
    if(length(map) != length(snpinfo))
      stop("length(map) [", length(map), "] != length(snpinfo) [",
           length(snpinfo), "]")
    
    if(nrow(snp_results) != length(unlist(map)))
      stop("nrow(snp_results) [", nrow(snp_results), "] != length(unlist(map)) [",
           length(unlist(map)), "]")
    
    lodindex <- split(1:nrow(snp_results), rep(names(map), vapply(map, length, 0)))
    
    result <- NULL
    for(i in seq(along=map)) {
      revindex <- rev_snp_index(snpinfo[[i]])
      
      map[[i]] <- snpinfo[[i]]$pos
      names(map[[i]]) <- snpinfo[[i]]$snp
      result <- rbind(result, unclass(snp_results)[lodindex[[i]],,drop=FALSE][revindex,,drop=FALSE])
      rownames(result) <- snpinfo[[i]]$snp
    }
    
    list(lod=result,
         map=map)
  }

# snpinfo to map
snpinfo_to_map <-
  function(snpinfo)
  {
    uindex <- sort(unique(snpinfo$index))
    if(any(snpinfo$index < 1 | snpinfo$index > nrow(snpinfo)))
      stop("snpinfo$index values outside of range [1, ",
           nrow(snpinfo), "]")
    
    uchr <- unique(snpinfo$chr)
    chr <- factor(snpinfo$chr, levels=uchr)
    
    map <- split(snpinfo$pos, chr)
    snp <- split(snpinfo$snp, chr)
    index <- split(snpinfo$index, chr)
    for(i in seq(along=map)) {
      u <- unique(index[[i]])
      map[[i]] <- map[[i]][u]
      names(map[[i]]) <- snp[[i]][u]
    }
    
    names(map) <- uchr
    
    map
  }

# reverse index
rev_snp_index <-
  function(snpinfo)
  {
    index_spl <- split(1:nrow(snpinfo), snpinfo$index)
    revindex <- rep(seq(along=index_spl), vapply(index_spl, length, 1))
    revindex[unlist(index_spl)] <- revindex
    
    revindex
  }


######### PERFORMING ASSOCIATION MAPPING ##############

# assoc.3 <- assoc_mapping(probs = probs, pheno = ins_test_clinphenos, idx = 3,
#                          addcovar = addcovar, k = K, markers = markers,
#                          chr = 11, start = 82, end = 85, ncores = 4)
# genes = data.frame(chr   = annot.mrna$chr,
#                    start = annot.mrna$start,
#                    stop  = annot.mrna$end,
#                    strand = annot.mrna$strand,
#                    Name   = annot.mrna$symbol, stringsAsFactors = F)
# quartz()
# layout(matrix(1:2, 2, 1))
# par(plt = c(0.1, 0.99, 0.1, 0.9))
# plot_snpasso_experimental(scan1output = assoc.3[[1]], snpinfo = assoc.3[[2]],
#                           drop.hilit = 1, xlim = c(82,85))
# 
# plot_genes(genes[genes$chr == 11 & genes$start > 82e6 & genes$stop < 85e6,],
#            xlim = c(82, 85))