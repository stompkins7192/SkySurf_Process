invisible(library(NISTunits, quietly = TRUE))
invisible(library(pracma, quietly = TRUE))
invisible(library(Rcpp, quietly = TRUE))
invisible(library(RANN, quietly = TRUE))
invisible(library(celestial, quietly = TRUE))
invisible(library(foreach, quietly = TRUE))
invisible(library(magicaxis, quietly = TRUE))
invisible(library('data.table', quietly = TRUE))
invisible(library('Rfits', quietly = TRUE))
invisible(library(ProFound, quietly = TRUE))

invisible(library(Rfits, quietly = TRUE))
invisible(library(Cairo, quietly = TRUE))
invisible(library(doParallel, quietly = TRUE))
library(Rwcs, quietly = TRUE)


inputargs = commandArgs(TRUE)
filter = as.character(inputargs[1])
visit_start = as.numeric(inputargs[2]) #visit number to start at
visit_end = as.numeric(inputargs[3]) #visit number to end at
slope = as.numeric(inputargs[4]) # Slope of Mag vs log10(R50) space divider
yint = as.numeric(inputargs[5]) #y intercept of the relationship
bc = as.numeric(inputargs[6])   # brightness cut, all objects brighter than bd will be flagged.
Ncpu = as.numeric(inputargs[7]) # Number of CPU's to use. Choose wisely.


#visit = 1416

registerDoParallel(cl <- makeCluster(Ncpu, outfile = ""))
j = 1





cat(paste("Total # of Visits to Process =", (visit_end - visit_start + 1)), "\n",
    "On", Ncpu, "Cores ", "\n")

foreach(j = visit_start:visit_end, .verbose = FALSE,
        .packages = c("pracma" ,"NISTunits", "RANN", "Rcpp","Rfits", "Cairo", "stringr",
                      "magicaxis", "data.table",
                      "Rwcs", "MASS", "ProFound", "celestial", "plotrix")) %dopar% {
                        
                        visit = j
start_time = Sys.time()
filter = filter

path = paste0("~/../../../media/scott/SSD/dds/", filter, "/", filter, "_drizzle_outputs/run_2/Visit_", visit,"/drizzle_products/")
wht = Rfits_point(filename = paste0(path, "Visit_", visit, "_drz_wht.fits"))

p2 = strsplit(path, paste0(filter, "_drizzle_outputs"))
p2 = unlist(p2)
p2 = p2[1]




#power law defining mask radius
maskrad = function(x){
  m = 3.845 - 0.1289*x
  m = 10^m
  return(m*2)
}

#return objects above line to catalog
f1 = function(x){
  y = slope*x + yint
  return(y)
}

#
pd = function(v1,v2){
  s = (v1 - v2) / ((v1 + v2)/2)
  return(100*s)
}


stars = read.csv(file = paste0(p2, "/profound_outputs/Gaia_CSVs/", filter, "_", visit, "_stars.csv"))
objects = read.csv(paste0(p2, "profound_outputs/Image_Info/segstats_",
                          filter, "_", visit, ".csv"))
#draw line

b = which(log10(stars$R50) < f1(stars$ABmag))

stars = stars[b,]
r = maskrad(stars$ABmag)
bright = which(objects$mag < bc)
r2 = maskrad(objects$mag[which(objects$mag < bc)])



test = profoundApplyMask(wht[,]$imDat, mask = "disc",
                         xcen = c(stars$xcen, objects$xcen[bright]), ycen = c(stars$ycen, objects$ycen[bright]),
                         xsize = c(maskrad(stars$ABmag), maskrad(objects$mag[bright])),
                         ysize =c(maskrad(stars$ABmag), maskrad(objects$mag[bright])))



test$image[which(is.na(test$image))] = 0
a = length(which(test$image != 0))*pixarea(wht) / 3600^2
{
  #CairoPNG(filename = paste0("/home/scott/Pictures/08-23/wfc3ir_history/testdir/",
   #                          visit, ".png"),
    #       width=dim(test$image)[1], height=dim(test$image)[2],units="px")
  #magimage(wht[,]$imDat, qdiff = TRUE)
  #a = length(which(wht[,]$imDat != 0))*pixarea(wht) / 3600^2
  #legend("topright", legend = c(paste("area =", signif(a, 5), "deg^2")),
   #      text.col = "cyan", bty = "n", border = "n",
    #    cex = 6)
  
  
  
  #dev.off() 
}


{
#  CairoPNG(filename = paste0("/home/scott/Pictures/08-23/wfc3ir_history/testdir/",
 #                            visit, "_sm.png"),
  #         width=dim(test$image)[1], height=dim(test$image)[2],units="px")
  #magimage(test$image, qdiff = TRUE)

  
  #legend("topright", legend = c(paste("area =", signif(a, 5), "deg^2")),
   #      text.col = "cyan", bty = "n", border = "n",
    #     cex = 6)
  
  

  


  below = which(log10(objects$R50) < f1(objects$mag))
  above = which(log10(objects$R50) > f1(objects$mag))
  z = which((objects$mag) < bc)#flag objects brighter than bc with "3"
  objects$secondary_flag = 0
  objects$starflag[below] = 1
  objects$starflag[above] = 0
 
  
  for(i in 1:length(objects$segID)){
    
    if(test$mask[objects$xcen[i], objects$ycen[i]] !=0){
      
      objects$secondary_flag[i] = 2
      
    }
    
  }
  flags = which(objects$secondary_flag == 2)
  objects$secondary_flag[z] = 3
  #points(objects$xcen[flags], objects$ycen[flags], col = "cyan", pch = "*",
   #      cex = 5)
  
  
}
#cat(paste("original area = ", length(which(wht[,]$imDat != 0))*pixarea(wht) / 3600^2, "\n"))
#cat(paste("new area = ",a, "\n"))
#overwrite segstats with new column added
write.csv(objects, file = paste0(p2, "profound_outputs/Image_Info/segstats_",
                filter, "_", visit, ".csv"), row.names = FALSE)


orig_area = length(which(wht[,]$imDat != 0))*pixarea(wht) / 3600^2
new_area = a
objects$area = a

redo = cbind(filter, as.numeric(visit), 
             as.numeric(orig_area), as.numeric(new_area),
            length(unique(stars$segstat_ID)), length(flags))



redo = data.frame(redo)

names = c("filter", "visit", "orig_area", "new_area", "N_stars", "N_flag")
colnames(redo) = names


filename =  paste0(strsplit(path, filter)[[1]][1], filter, "/profound_outputs/Image_Info/new_areas_",
                   filter, ".csv")

if(visit == visit_start){
  
  write.table(redo, file = filename, append = FALSE, row.names = FALSE, col.names = TRUE,
              sep = ",")
  
} else{
  
  write.table(redo, file = filename, append = TRUE, row.names = FALSE, col.names = FALSE,
              sep = ",")
  
}
end_time = Sys.time()
end_time = (difftime(end_time, start_time, units="secs")[[1]])
cat(paste("Execution Time For Visit", visit, "=", signif(end_time, 4), "seconds", "\n"))

}




#end doParallel
invisible(stopCluster(cl))

