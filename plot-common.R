library(ggplot2)

#
# Helpful functions
#
mmperpt = 1/72 * 25.4
pt = mmperpt

percentlog10y = scale_y_log10(
    breaks=c(1,0.1,0.01),
    labels=c(100,10,1))

margin_none = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")

default_theme = theme(
  text = element_text(size=9,margin=margin_none),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(size=0.75*pt,color="#808080",linetype="13"),
  axis.ticks = element_line(size=1*pt),
  axis.ticks.length = unit(1.5,"pt"),
  axis.title.y.right=element_text(margin=margin_none),
  axis.title.y.left=element_text(margin=margin_none),
  axis.title.x.top=element_text(margin=margin_none),
  axis.title.x.bottom=element_text(margin=margin_none),
  title=element_text(size=9,margin=margin_none),
  #plot.margin=unit(c(0,0.1,0,0),c("cm","cm","cm","cm")),
  plot.margin=margin_none,
  panel.spacing=unit(1,"mm"),
  plot.title=element_text(margin=margin(t=0,b=0.5,r=0,l=0,unit="mm")),
  axis.text = element_text(color="#000000",size=7,margin=margin_none),
  strip.background=element_rect(color=NA,fill=NA),
  strip.text.x=element_text(size=9,margin=margin_none),
    #margin=margin(t=unit(3,"mm"),r=0,b=unit(1,"mm"),l=0)),
  panel.background = element_rect(
    fill="transparent", colour = "#000000",size=1*pt),
  plot.background = element_rect(
    fill="transparent", colour = NA))

# make_discrete_cdf - create data frame with points for the top and bottom of
# discontinuities in a discrete cumulative distribution function
make_discrete_cdf = function(xs, ys) {
  return(data.frame(
    x=rep(xs,each=2)[1:(2*length(xs)-1)],
    y=c(1,rep(ys[1:(length(ys)-1)],each=2)))) }

# make_ecdf_df - create a dataframe in the same form as make_discrete_cdf for
# the cdf of an empirical sample
make_ecdf_df = function(vs) {
  ecd = ecdf(vs)
  xs = environment(ecd)$x
  ys = environment(ecd)$y
  return(make_discrete_cdf(xs, 1-ys)) }
