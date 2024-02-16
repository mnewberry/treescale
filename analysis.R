source("dplfit.R")
source("plot-common.R")
library(data.table)

# Steps to reproduce the analysis of
# Fractal scaling and the aesthetics of trees
# (c) Jingyi Gao and Mitchell Newberry 2024 
# licensed under GPL version 3.0, see LICENSE

set.seed(2)

read_sample = function(nlsv, work, rep)
  data.frame(x=read.delim(nlsv,header=FALSE)$V1, work=work, rep=rep)

remax = function(df,n=500) {
  xs = sort(df$x,decreasing=TRUE)[1:min(n,nrow(df))] #take top NNN vals
  xs = xs / min(xs)
  return(data.frame(x=xs,work=df$work[1],rep=df$rep[1])) }

# Renormalize measurements to make the max value equal truemax
normmax = function(truemax, df) {
  df$x = df$x*truemax/max(df$x)
  return(df) }

# Renormalize measurements to make the min value equal truemin
normmin = function(truemin, df) {
  df$x = df$x/truemin
  return(df) }

# Assemble all data, renormalizing pixel values into physical units, or
# physical units into multiples of xmin in the case of real tree measurements
# from Bentley et al. 2013.
samples = rbind(
  transform(read.delim("samples/grove5.tsv"),
    x = diameter/0.5, work="Simulated trees", rep=alpha)[,c("x","work","rep")],
  normmin(0.017, subset(
    read_sample("samples/pinon.nlsv", "Real trees", "Piñon"), x > 0.0075)),
  normmin(0.012, subset(
    read_sample("samples/balsa.nlsv", "Real trees", "Balsa"), x > 0.006)),
  normmin(0.13, subset(
    read_sample("samples/ponderosa.nlsv", "Real trees", "Ponderosa"),x>0.065)),
  normmax(24.5, read_sample("samples/gkm.nlsv",  "Klimt", "M")),
  normmax(24.6, read_sample("samples/gkj.nlsv",  "Klimt", "J")),
  normmax(25.73, read_sample("samples/mgm.nlsv",  "Goshun", "M")),
  normmax(25.73, read_sample("samples/mgj.nlsv",  "Goshun", "J")),
  normmax(10.1, read_sample("samples/ss1aj.nlsv","Sidi Saiyyed", "1aJ")),
  normmax(13.9, read_sample("samples/ss1am.nlsv","Sidi Saiyyed", "1aM")),
  normmax(14.7, read_sample("samples/ss2aj.nlsv","Sidi Saiyyed", "2aJ")),
  normmax(16.4, read_sample("samples/ss2am.nlsv","Sidi Saiyyed", "2aM")),
  normmax(9.42, read_sample("samples/gtm.nlsv",  "Mondrian", "M")),
  normmax(7.39, read_sample("samples/gtj.nlsv",  "Mondrian", "J")),
  normmax(9.3, read_sample("samples/gta.nlsv",  "Mondrian", "A")))

# General function to make tail distribution pots and conduct power law fits.
plplot <- function (xms, lambdas, data, title) {
  fit = function(xm, lambda) { return(function (data) {
    lengths = data$x
    fitres = dplfit_bin_se(lengths, lambda, xm)
    alpha = fitres$alpha
    se = fitres$se
    label=sprintf("alpha == %0.1f %%+-%% %0.1f",alpha,se)
    xmoffset = log(1-length(lengths[lengths < xm])/length(lengths))
    yoffset = length(lengths[lengths >= xm])/length(lengths)
    basemult = log(xm,10)/log(xm) 
    return(data.frame(est = "Discrete MLE",
      xm=xm, lambda=lambda,alpha=alpha, label=label,
      intercept=log10(yoffset)-(-alpha*log10(xm)),
      slope=-alpha,rep=data$rep[1])) }) }
  ecdfs = rbindlist(by(data, data$rep,
    function (df) cbind(make_ecdf_df(df$x),
       work=df$work[1],rep=df$rep[1])))
  fits=data.frame()
  for(i in 1:length(lambdas)) {
    xm = xms[i] ; lambda = lambdas[i] ;
    fits = rbind(fits, rbindlist(by(data, data$rep, fit(xm, lambda)))) }

  return(ggplot(fits) +
    geom_vline(xintercept=xm,size=1*pt,linetype="22") +
    geom_abline(data=fits,aes(intercept=intercept,slope=slope,color=rep),
      size=1*pt) +
    geom_line(data=ecdfs,aes(x=x,y=y,group=rep,color=rep),size=1*pt) +
    percentlog10y + 
    scale_x_log10() + 
    labs(y="percent of branches exeeding d", x="branch diameter, d") +
    ggtitle(title) +
    default_theme + theme(legend.pos = "none"))
}

# Artwork samples
label=function(x) sprintf("alpha == %0.1f",x)
ylim=c(1.4/500,1)
plplot(1.8,2,
    subset(samples, work=="Sidi Saiyyed"),"Sidi Saiyyed") +
  geom_text(aes(group=rep,color=rep,
      vjust=1.4*as.numeric(factor(rep))-0.8,label=label),
    x=log10(max(subset(samples, work=="Sidi Saiyyed")$x)),
    y=0,hjust=1,size=7*pt,parse=TRUE) +
  xlab(expression(phantom("()"))) +
  coord_fixed(ylim=ylim)
plplot(0.3,2,
    subset(samples, work=="Goshun"),"Goshun") +
  geom_text(aes(group=rep,color=rep,
      vjust=1.4*as.numeric(factor(rep))-0.8,label=label),
    x=log10(max(subset(samples, work=="Goshun")$x)),
    y=0,hjust=1,size=7*pt,parse=TRUE) +
  xlab(expression(paste("branch diameter, ", d, " (cm)"))) +
  coord_fixed(ylim=ylim)
plplot(2.4,2,
    subset(samples, work=="Klimt"),"Klimt") +
  geom_text(aes(group=rep,color=rep,
      vjust=1.4*as.numeric(factor(rep))-0.8,label=label),
    x=log10(max(subset(samples, work=="Klimt")$x)),
    y=0,hjust=1,size=7*pt,parse=TRUE) +
  xlab(expression(phantom("()"))) +
  coord_fixed(ylim=ylim)

# Real and simulated trees
ylim=c(0.0015,1)
plplot(1,2,
    subset(samples, work=="Simulated trees"),"Simulated") +
  geom_text(aes(group=rep,color=rep,
      vjust=1.4*as.numeric(factor(rep))-1.1,label=label),
    x=log10(max(subset(samples, work=="Simulated trees")$x)),
    y=0,hjust=1,size=7*pt,parse=TRUE) +
  xlab(expression(phantom("()"))) +
  coord_fixed(ylim=ylim)
plplot(1,2,
    subset(samples, work=="Real trees"),"Measured") +
  geom_text(aes(group=rep,color=rep,
      vjust=1.4*as.numeric(factor(rep))-0.8,label=label),
    x=log10(max(subset(samples, work=="Real trees")$x)),
    y=0,hjust=1,size=7*pt,parse=TRUE) +
  xlab(expression(paste("branch diameter, ", d, " (normalized)"))) +
  coord_fixed(ylim=ylim) 

# Mondrian Grey tree
ylim=c(0.007,1)
plplot(1.6,2,
    subset(samples, work=="Mondrian"),"Mondrian") +
  geom_text(aes(group=rep,color=rep,
      vjust=1.4*as.numeric(factor(rep))-0.8,label=label),
    x=log10(max(subset(samples, work=="Mondrian")$x)),
    y=0,hjust=1,size=7*pt,parse=TRUE) +
  xlab(expression(paste("branch diameter, ", d, " (cm)"))) +
  coord_fixed(ylim=ylim)
