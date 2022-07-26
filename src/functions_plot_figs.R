my_custom_cor_gene <- function(data, mapping, thres.tmp, color = 'black', sizeRange = c(1, 5), ...) {

  # get the x and y data to use the other code
  x <-GGally::eval_data_col(data, mapping$x)
  y <-GGally::eval_data_col(data, mapping$y)

  ct <- cor.test(x,y)
  sig <- symnum(
    ct$p.value, corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("***", "**", "*", ".", " ")
  )
  ct.symbol <- c("***", "**", "*", ".", " ")[(ct$p.value  < c(0.001, 0.01, 0.05, 0.1, 1) )][1]
  r <- unname(ct$estimate)
  rt <- format(r, digits=2)[1]

  ## odds ratio
  exp.x <- 10^(-x)
  qval.x <- p.adjust(exp.x, method = 'fdr')
  exp.y <- 10^(-y)
  qval.y <- p.adjust(exp.y, method = 'fdr')

  ###
  m.tmp <- cbind(qval.x, qval.y)
 # thres.tmp <- 0.2
  a <- as.numeric( sum( rowSums(m.tmp < thres.tmp) == 2 ) ) +1
  b <- as.numeric( sum( ((qval.x<thres.tmp)*1 - (qval.y<thres.tmp)*1) == 1  ) ) +1
  c <- as.numeric( sum( ((qval.x<thres.tmp)*1 - (qval.y<thres.tmp)*1) == -1  ) ) +1
  d <- as.numeric( sum( rowSums(m.tmp < thres.tmp) == 0 ) ) +1

  #or.res <- oddsratio(as.integer(a),as.integer(b),as.integer(c),as.integer(d) )
  or.est <- format( d/b/c*a, digits=2)
  or.pval <- fisher.test(x=matrix(c(a,c,b,d),2,2))$p.value
  print(or.pval)
  or.symbol <- c("***", "**", "*", ".", " ")[( or.pval < c(0.001, 0.01, 0.05, 0.1, 1) )][1]
  # since we can't print it to get the strsize, just use the max size range
  cex <- max(sizeRange)

  # plot the cor value
  ggally_text(
    label = paste0( 'Corr: ',as.character(rt), ct.symbol,'\n', 'OR: ', as.character(or.est), or.symbol),
    mapping = aes(),
    xP = 0.5, yP = 0.5,
    size = 5,
    color = color,
    ...
  ) +
    # add the sig stars
    geom_text(
      aes_string(
        x = 0.8,
        y = 0.8
      ),
      label = " ",
      size = 1,
      color = color,
      ...
    ) +
    # remove all the background stuff and wrap it with a dashed line
    theme_classic() +
    theme(
      panel.background = element_rect(
        color = color,
        linetype = "longdash"
      ),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank()
    )
}

my_custom_cor_meth <- function(data, mapping, thres.tmp, color = 'black', sizeRange = c(1, 5), ...) {

  # get the x and y data to use the other code
  x <-GGally::eval_data_col(data, mapping$x)
  y <-GGally::eval_data_col(data, mapping$y)

  ct <- cor.test(x,y)
  sig <- symnum(
    ct$p.value, corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("***", "**", "*", ".", " ")
  )
  ct.symbol <- c("***", "**", "*", ".", " ")[(ct$p.value  < c(0.001, 0.01, 0.05, 0.1, 1) )][1]
  r <- unname(ct$estimate)
  rt <- format(r, digits=2)[1]

  ## odds ratio
  exp.x <- 10^(-x)
  p0.x <- pi0.est(exp.x)$p0
  qval.x <- qvalue.cal(exp.x, p0.x)
  exp.y <- 10^(-y)
  p0.y <- pi0.est(exp.y)$p0
  qval.y <- qvalue.cal(exp.y, p0.y)

  ###
  m.tmp <- cbind(qval.x, qval.y)
  #thres.tmp <- 0.2
  a <- as.numeric( sum( rowSums(m.tmp < thres.tmp) == 2 ) ) +1
  b <- as.numeric( sum( ((qval.x<thres.tmp)*1 - (qval.y<thres.tmp)*1) == 1  ) ) +1
  c <- as.numeric( sum( ((qval.x<thres.tmp)*1 - (qval.y<thres.tmp)*1) == -1  ) ) +1
  d <- as.numeric( sum( rowSums(m.tmp < thres.tmp) == 0 ) ) +1

  #or.res <- oddsratio(as.integer(a),as.integer(b),as.integer(c),as.integer(d) )
  or.est <- format( d/b/c*a, digits=2)
  or.pval <- fisher.test(x=matrix(c(a,c,b,d),2,2))$p.value
  print(or.pval)
  or.symbol <- c("***", "**", "*", ".", " ")[( or.pval < c(0.001, 0.01, 0.05, 0.1, 1) )][1]
  # since we can't print it to get the strsize, just use the max size range
  cex <- max(sizeRange)

  # plot the cor value
  ggally_text(
    label = paste0( 'Corr: ',as.character(rt), ct.symbol,'\n', 'OR: ', as.character(or.est), or.symbol),
    mapping = aes(),
    xP = 0.5, yP = 0.5,
    size = 5,
    color = color,
    ...
  ) +
    # add the sig stars
    geom_text(
      aes_string(
        x = 0.8,
        y = 0.8
      ),
      label = " ",
      size = 1,
      color = color,
      ...
    ) +
    # remove all the background stuff and wrap it with a dashed line
    theme_classic() +
    theme(
      panel.background = element_rect(
        color = color,
        linetype = "longdash"
      ),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank()
    )
}

de_pattern_show <- function(dp.ix){
  
  if(grepl('indep_all', dp.ix)){
    
    plot(NULL, type = 'l', lty = 1,xlim=c(0,4),ylim = c(0,2), lwd = 6,
         frame.plot=F, yaxt = 'n', xaxt = 'n', ylab='',xlab='', ann = FALSE)
    
    segments(x0=c(1,2,3,4), x1=c(1,2,3,4), y0=0.5, y1=1.5, lwd = 6)
    
    text(x=0.1, y = 1, labels= 'a', cex =4)
    
  }else if(grepl('corr_single', dp.ix)){
    
    plot(NULL, type = 'l', lty = 1,xlim=c(0,4),ylim = c(0,2), lwd = 6,
         frame.plot=F, yaxt = 'n', xaxt = 'n', ylab='',xlab='', ann = FALSE)
    
    segments(x0=c(1),x1=c(4),y0=1.5,y1=1.5, lwd = 6)
    segments(x0=c(1,2,3,4),x1=c(1,2,3,4),y0=0.5,y1=1.5, lwd = 6)
    
    text(x=0.1, y = 1, labels= 'b', cex =4)
    
  }else if(grepl('indep_first', dp.ix)){
    
    plot(NULL, type = 'l', lty = 1,xlim=c(0,4),ylim = c(0,2), lwd = 6,
         frame.plot=F, yaxt = 'n', xaxt = 'n', ylab='',xlab='', ann = FALSE)
    
    segments(x0=c(3),x1=c(4),y0=1.5,y1=1.5, lwd = 6)
    segments(x0=c(1,2,3,4),x1=c(1,2,3,4),y0=0.5,y1=1.5, lwd = 6)
    
    text(x=0.1, y = 1, labels= 'c', cex =4)
    
  }else if(grepl('indep_second',dp.ix)){
    
    plot(NULL, type = 'l', lty = 1,xlim=c(0,4),ylim = c(0,2), lwd = 6,
         frame.plot=F, yaxt = 'n', xaxt = 'n', ylab='',xlab='', ann = FALSE)
    
    segments(x0=c(1),x1=c(2),y0=1.5,y1=1.5, lwd = 6)
    segments(x0=c(1,2,3,4),x1=c(1,2,3,4),y0=0.5,y1=1.5, lwd = 6)
    
    text(x=0.1, y = 1, labels= 'd', cex =4)
    
  }else if(grepl('indep_group',dp.ix)){
    
    plot(NULL, type = 'l', lty = 1,xlim=c(0,4),ylim = c(0,2), lwd = 6,
         frame.plot=F, yaxt = 'n', xaxt = 'n', ylab='',xlab='', ann = FALSE)
    
    segments(x0=c(1,3),x1=c(2,4),y0=1.5,y1=1.5, lwd = 6)
    segments(x0=c(1,2,3,4),x1=c(1,2,3,4),y0=0.5,y1=1.5, lwd = 6)
    
    text(x=0.1, y = 1, labels= 'e', cex =4)
    
  }else if(grepl('corr_tree',dp.ix)){
    
    plot(NULL, type = 'l', lty = 1,xlim=c(0,4),ylim = c(0,2), lwd = 6,
         frame.plot=F, yaxt = 'n', xaxt = 'n', ylab='',xlab='', ann = FALSE)
    
    segments(x0=c(1,3),x1=c(2,4),y0=1,,y1=1, lwd = 6)
    segments(x0=c(1,2,3,4),x1=c(1,2,3,4),y0=0.5,y1=1, lwd = 6)
    segments(x0=c(1.5),x1=c(3.5),y0=1.5,,y1=1.5, lwd = 6)
    segments(x0=c(1.5,3.5),x1=c(1.5,3.5),y0=1,y1=1.5, lwd = 6)
    
    text(x=0.1, y = 1, labels= 'f', cex =4)
    
  }
}








