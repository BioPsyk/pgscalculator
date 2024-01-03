args <- commandArgs(trailingOnly=TRUE)

library(data.table)
poX <- args[1]
clX <- args[2]
clhX <- args[3]
po <- fread(poX)
cl <- fread(clX)

me <- merge(po, cl[,c("RSID","B")], by.x="Name", by.y="RSID")
setorder(me, Chrom, Position)

# support variables
me$pos <- 1:length(me$Position)
me$blownUp <- 1*(abs(me$A1Effect) > 10*abs(me$B) & abs(me$A1Effect) > quantile(abs(me$B), 0.95))

# write plot table
write.table(me, file = stdout(), row.names = FALSE, quote = FALSE)

# plot
plot(me$Position, me$A1Effect, col=c("black", "grey")[1 +me$Chrom %% 2])
points(me$Position[me$blownUp == 1], me$A1Effect[me$blownUp == 1], col=2, pch=16)
abline(h=quantile(abs(me$B), 0.95),col=2, lty=2)
abline(h=-quantile(abs(me$B), 0.95),col=2, lty=2)


