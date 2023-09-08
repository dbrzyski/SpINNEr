## Clear the workspace
rm(list = ls())

## Libraries
library(ggplot2)

## Download data
# datapath <- "C:\\Users\\brzyski\\Dropbox\\NuclearNorm\\NuclearLasso"
# savepath <- "C:\\Users\\brzyski\\Dropbox\\NuclearNorm\\NuclearLasso\\JournalVersion\\figures"
# setwd(datapath)
load("G:\\Dropbox\\NuclearNorm\\NuclearLasso\\simulationCode\\Results\\SimulationResults.Rdata")

## Plot settings
plot_labels <- c("SpINNEr", "ElasNet", "Nuclear", "LASSO", "Ridge", "CPM")
plot_colors <- c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "#F09438", "#666666")
plot_shape <- c(25, 22, 23, 24, 21, 4)
## X-axis Labels
x_label <- rep(c(expression(2^-3), expression(2^-2), expression(2^-1), expression(2^0), expression(2^1), 
                 expression(2^2), expression(2^3), expression(2^4), expression(2^5)), 6)

## PLOT: scenario 1
df    <- results_setting1
xbeg  <-  -2.95
xend  <-  4.8
ybeg  <-  0.05
yend  <-  2.05

p <- ggplot(df, aes(k, err_mean, colour = type, shape = type, fill = type,  group = type))
p <- p + geom_line(mapping = aes(x = k, y = err_mean, color = type), size = 1, linetype = 1)
p <- p + scale_color_manual(values = plot_colors, labels = plot_labels)
p <- p + scale_shape_manual(values = plot_shape, labels = plot_labels)
p <- p + scale_fill_manual(name = "Method", values = plot_colors, labels = plot_labels)
p <- p + geom_errorbar(aes(ymin = err_lower, ymax = err_upper), width = 0.25) 
p <- p + geom_point(df, mapping = aes(x = k, y = err_mean,colour = type, shape = type), size = 6)
p <- p + coord_cartesian( xlim = c(xbeg,  xend),  ylim = c(ybeg, yend) )
p <- p + labs(title = "", x = "s", y = 'Estimated MSEr', color = "Method", shape = "Method")
p <- p + theme_bw(base_size = 25)
# p <- p + theme( legend.position = c(0.242, 0.434), legend.justification = c("left", "top"))
p <- p + theme( legend.position = "none")
p <- p + theme(legend.text = element_text(size = 18))
p <- p + theme(legend.background = element_rect(colour = 'black'))
p <- p + scale_x_continuous("s", labels = x_label, breaks = df$k)
p <- p + theme(panel.grid.major = element_line(size = 1), panel.grid.minor = element_line(size = 0.7))
p <- p + theme(axis.text.x = element_text(size = 20))
p <- p + theme(axis.text.y = element_text(size = 20, angle = 0))
p <- p + theme(axis.title.y = element_text(vjust = -3))
p <- p + theme(axis.title.x = element_text(vjust = 7))
p <- p + theme(axis.ticks.length = unit(-0.25, "cm"), axis.text.x = element_text(margin=unit(c(0.35,0.35,0.35,0.35), "cm")), axis.text.y = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm")) )
p <- p + theme(plot.margin=unit(c(-1.15, 0.15, -0.45, -.2),"cm"))
p <- p + theme(legend.key.width = unit(1.05,"cm"))
p <- p + theme(legend.key.height = unit(0.55,"cm"))
# setwd(savepath)
pdf("scenario1.pdf", width = 7, height = 6, pointsize = 10)
plot(p)
dev.off()

#-------------------------------------------------------------------------------------------
## PLOT: scenario 2
df    <- results_setting2
xbeg  <-  100
xend  <-  300
ybeg  <-  0.00
yend  <-  2.75

p <- ggplot(df, aes(n, err_mean, colour= type, shape= type, fill = type,  group=type))
p <- p + geom_line(mapping = aes(x = n, y = err_mean, color = type), size = 1, linetype = 1)
p <- p + scale_color_manual(values=plot_colors, labels = plot_labels)
p <- p + scale_shape_manual(values = plot_shape, labels = plot_labels)
p <- p + scale_fill_manual(name="Method", values=plot_colors, labels = plot_labels)
p <- p + geom_errorbar(aes(ymin = err_lower, ymax = err_upper), width = 7) 
p <- p + geom_point(df, mapping=aes(x = n, y = err_mean,colour= type, shape = type), size = 3)
p <- p + coord_cartesian( xlim = c(xbeg,  xend),  ylim = c(ybeg, yend) )
p <- p + labs(title="", x = "n", y = 'Estimated MSEr', color = "Method", shape = "Method")
p <- p + theme_bw(base_size = 25)
# p <- p + theme( legend.position = c(0.055, 0.43), legend.justification = c("left", "top"))
p <- p + theme( legend.position = "none")
p <- p + theme(legend.text=element_text(size=18))
p <- p + theme(legend.background = element_rect(colour = 'black'))
p <- p + scale_x_continuous('n', labels = as.character(df$n), breaks = df$n)
p <- p + theme(panel.grid.major = element_line(size = 1), panel.grid.minor = element_line(size = 0.7))
p <- p + theme(axis.text.x=element_text(size=20))
p <- p + theme(axis.text.y=element_text(size=20, angle=0))
p <- p + theme(axis.title.y = element_text(vjust=-2.5))
p <- p + theme(axis.title.x = element_text(vjust=7))
p <- p + theme(axis.ticks.length=unit(-0.25, "cm"), axis.text.x = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm")), axis.text.y = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm")) )
p <- p + theme(plot.margin=unit(c(-1.15, 0.15, -0.45, -.2),"cm"))
p <- p + theme(legend.key.width = unit(1.05,"cm"))
p <- p + theme(legend.key.height = unit(0.55,"cm"))

# setwd(savepath)
pdf("scenario2.pdf", width=7, height=6, pointsize=10)
plot(p)
dev.off()

#-------------------------------------------------------------------------------------------
## PLOT: scenario 3: original scale
df    <- results_setting3
xbeg  <-  -2.95
xend  <-  4.8
ybeg  <-  0.05
yend  <-  367

p <- ggplot(df, aes(k, err_mean, colour= type, shape= type, fill = type,  group=type))
p <- p + geom_line(mapping = aes(x = k, y = err_mean, color = type), size = 1, linetype = 1) 
p <- p + scale_color_manual(values=plot_colors, labels = plot_labels)
p <- p + scale_shape_manual(values = plot_shape, labels = plot_labels)
p <- p + scale_fill_manual(name="Method", values=plot_colors, labels = plot_labels)
p <- p + geom_errorbar(aes(ymin = err_lower, ymax = err_upper), width = 0.25) 
p <- p + geom_point(df, mapping=aes(x = k, y = err_mean,colour= type, shape = type), size = 3)
p <- p + coord_cartesian( xlim = c(xbeg,  xend),  ylim = c(ybeg, yend) )
p <- p + labs(title="", x = "s", y = 'Estimated MSEr', color = "Method", shape = "Method")
p <- p + theme_bw(base_size = 25)
p <- p + theme( legend.position = c(0.635, 0.82), legend.justification = c("left", "top"))
p <- p + theme(legend.text=element_text(size=18))
p <- p + theme(legend.background = element_rect(colour = 'black'))
p <- p + scale_x_continuous("s", labels = x_label, breaks = df$k)
p <- p + theme(panel.grid.major = element_line(size = 1), panel.grid.minor = element_line(size = 0.7))
p <- p + theme(axis.text.x=element_text(size=20))
p <- p + theme(axis.text.y=element_text(size=20, angle=0))
p <- p + theme(axis.title.y = element_text(vjust=-2.5))
p <- p + theme(axis.title.x = element_text(vjust=7))
p <- p + theme(axis.ticks.length=unit(-0.25, "cm"), axis.text.x = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm")), axis.text.y = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm")) )
p <- p + theme(plot.margin=unit(c(-1.15, 0.15, -0.45, -.2),"cm"))
p <- p + theme(legend.key.width = unit(1.05,"cm"))
p <- p + theme(legend.key.height = unit(0.55,"cm"))

# setwd(savepath)
pdf("scenario3_original.pdf", width=7, height=6, pointsize=10)
plot(p)
dev.off()

#-------------------------------------------------------------------------------------------
## PLOT: scenario 3: log scale 
p <- p + scale_y_log10()

# setwd(savepath)
pdf("scenario3_logscale.pdf", width=7, height=6, pointsize=10)
plot(p)
dev.off()


