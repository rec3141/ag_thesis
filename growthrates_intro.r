### Calculate bacterial growth rates from plate reader data

library(growthrates) # calculates growth rates
library(reshape2) # melt/cast
library(readxl) # read XLS files

# setwd("~/Downloads/20200124-Research2-master")

# read in data
dat.in = read_excel("total_raw_data_2.xlsx", sheet = "Data")
dat.in$`Start Time` = rep(dat.in[[1,3]],nrow(dat.in))
dat.in$`TD` = as.numeric(difftime(dat.in$`Time Taken`,dat.in$`Start Time`,units = "hours"))
dat.in$`Time Taken` = NULL
dat.in$`Start Time` = NULL

# cut out bad data
dat.in = dat.in[!(dat.in$Temperature == 4 & dat.in$`Plate Number`== 3 & dat.in$TD > 40),]
dat.in = dat.in[!(dat.in$Temperature == 11 & dat.in$`Plate Number` == 8 & dat.in$TD > 40),]
dat.in = dat.in[!(dat.in$Temperature == 17 & dat.in$TD > 40),]
dat.in$`Plate Number` = NULL

# melt data
dat.all = reshape2::melt(dat.in,id = c("TD","Temperature","Time Point","Strain"))
colnames(dat.all) = c("time","temp","timepoint","strain","replicate","value")
dat.all = dat.all[,c("strain","replicate","temp","time","value")]

#plot all data
pdf(file="xyplot_strain_temps.pdf", width=24, height=24)
xyplot(value ~ time|strain+as.factor(temp), data = dat.all, groups = replicate, pch = 16, cex = 0.5)
dev.off()

## also: growthrates_demo.r

# fit splines to all data
many_spline_fits = all_splines(value ~ time | replicate + temp + strain,
                                data = dat.all, spar = 0.5)

many_spline_res = results(many_spline_fits)

# plot all fits
pdf(file="many_spline_fits.pdf",width=16,height=8)

par(mfrow = c(4, 8))
par(mar = c(2.5, 4, 2, 1))

plot(many_spline_fits)

xyplot(mumax ~ temp|strain, data = many_spline_res, layout = c(7, 4))

dev.off()

#output data
grd = many_spline_res[,c("strain","replicate","temp","mumax","r2")]
rownames(grd) = NULL

write.table(grd, file="growthrates.tsv",sep="\t",quote=F)

# run ratk.r to predict OGT using Ratkowsky model

