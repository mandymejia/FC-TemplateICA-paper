###############################################################
## FUNCTION FOR PLOTTING FC
###############################################################

#x = FC matrix
#zlim = color scale limits
plot_FC <- function(x, zlim=c(-1,1), title=NULL, cols=c('darkblue','turquoise','white','pink','red'), cols_rev=FALSE, break_by=0.5, digits_legend=1, cor=TRUE, lines = seq(size), labels = NULL, col_lines = 'black', lwd_lines = 1, cex = 0.8){

  require(grDevices)
  old_par <- par(no.readonly = TRUE)

  if(cols_rev) cols <- rev(cols)

  #set color scale and breaks
  breaks <- seq(zlim[1], zlim[2], length.out=100)
  levs <- format(round(seq(zlim[1], zlim[2], break_by), digits_legend), scientific=FALSE, digits = digits_legend)
  palfun <- grDevices::colorRampPalette(cols)
  pal <- palfun(100-1)

  #make plot
  layout(matrix(c(1,2,0,3), nrow=2, ncol=2), widths=c(5,1.2), heights=c(1.2,5))

  #title in top-left cell
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.3, title, cex = 2, col = "black")

  # matrix in bottom-left cell
  if(cor) diag(x) <- NA
  size <- ncol(x)
  x[x >= zlim[2]] <- (zlim[2]-0.0001) #truncate above
  x[x <= zlim[1]] <- (zlim[1]+0.0001) #trunate below
  par(mar=c(1,2,0,1))
  image(seq(size), seq(size), t(x[size:1,]), col=pal, breaks=breaks-1e-8, xaxt="n", yaxt="n", ylab="", xlab="")
  abline(h=size-lines+0.5, v=lines+0.5, col = col_lines, lwd=lwd_lines)
  if(!is.null(labels)){
    if(length(labels) == length(lines) + 1){
      nK <- length(labels)
      for(k in 1:nK){
        if(k == 1) start <- 0 else start <- lines[k-1]
        if(k == nK) end <- size else end <- lines[k]
        at_k <- start + (end-start)/2 + 0.5
        mtext(labels[k], side = 3, at = at_k, line=0, font=2, cex=cex)
        mtext(labels[k], side = 2, at = size - at_k + 1, line=0, font=2, cex=cex)
      }
    }
  }

  # color scale in bottom-right cell
  par(mar=c(1, 0.7, 0, 4))
  image.scale(mat, col=pal, breaks=breaks-1e-8, axis.pos=4, add.axis=FALSE)
  axis(4,at=levs, las=2)
  abline(h=levs)

  par(old_par)
}

#This function creates a color scale for use with the image()
#function. Input parameters should be consistent with those
#used in the corresponding image plot. The "axis.pos" argument
#defines the side of the axis. The "add.axis" argument defines
#whether the axis is added (default: TRUE)or not (FALSE).
image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, axis.pos=1, add.axis=TRUE, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  if(axis.pos %in% c(1,3)){ylim<-c(0,1); xlim<-range(breaks)}
  if(axis.pos %in% c(2,4)){ylim<-range(breaks); xlim<-c(0,1)}
  plot(1,1,t="n",ylim=ylim, xlim=xlim, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", ...)
  for(i in seq(poly)){
    if(axis.pos %in% c(1,3)){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(axis.pos %in% c(2,4)){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
  box()
  if(add.axis) {axis(axis.pos)}
}


