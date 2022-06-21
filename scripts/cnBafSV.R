library(ggplot2)
library(scales)
library(gtable)
library(grid)
library(DNAcopy)

args = commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=T)

# Read-depth
p = ggplot(data=x, aes(x=start, y=x[,6]))
p = p + geom_point(pch=21, size=0.5)
p = p + xlab("chr2")
p = p + ylab("Copy-number")
p = p + scale_x_continuous(labels=comma)
p = p + scale_y_continuous(labels=comma, breaks=0:12, limits=c(0,12))
p = p + theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave("cov.png", p, width=10, height=7)
if (length(args) == 1) { quit(); }

# Read-depth + SVs
sv = read.table(args[2], header=F)
colnames(sv) = c("chr","start","end","type","id")
p2 = p + geom_curve(data=sv, aes(x=start, xend=end, col=type), y=8, yend=8, curvature=-0.5)
p2 = p2 + labs(colour="SV type")
seg=segments.summary(segment(smooth.CNA(CNA(log(x[,6]), rep("chr2", nrow(x)), x$start, data.type="logratio", sampleid="tumor")), undo.splits="sdundo", undo.SD=2))
p3=p2 + geom_segment(data=seg, aes(x=loc.start, y=exp(seg.median), xend=loc.end, yend=exp(seg.median)), colour="darkorange")
p3=p3 + theme(legend.position="bottom")
ggsave("cov.png", p3, width=10, height=7)
if (length(args) == 2) { quit(); }

# Read-depth + SVs + BAF
baf = read.table(args[3], header=F)
colnames(baf) = c("pos", "baf")
q = ggplot(data=baf, aes(x=pos, y=baf))
q = q + geom_point(fill="black", colour="black", size=0.1, shape=21)
q = q + ylab("B-Allele Frequency") + xlab("chr2")
q = q + scale_x_continuous(labels=comma)
q = q + ylim(0,1)
q = q + theme(axis.text.x = element_text(angle=45, hjust=1))
p3 = p3 + xlim(min(baf$pos), max(baf$pos))
g1 = ggplotGrob(p3)
g2 = ggplotGrob(q)
g = rbind(g1, g2, size="first")
g$widths = unit.pmax(g1$widths, g2$widths)
ggsave(g, file="cov.png", width=10, height=7)
print(warnings())
