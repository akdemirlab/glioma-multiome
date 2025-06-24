
load({'{output_path_from_ascat.sc.R}/result_object_final.Rda')

output<-"**Out_path**"

outdir<-paste0(output,'/qc/')
dir.create(outdir)

project_name<-"final_p29"


nrec<-16
dist<-0.03
probs<-0.1

#this is a modified version of the default ascat.sc filtering function.

getFilters_update <- function (res, probs = probs, outdir = outdir, projectname = project_name, thresholdNrec = nrec, distance_filter = dist) 
{
    filters <- NULL
    try({
        allT <- res$allTracks.processed
        allS <- res$allSolutions
        getloess <- function(qu, nr) {
            nms <- paste0("n", 1:length(qu))
            names(qu) <- nms
            quo <- qu[order(nr, decreasing = F)]
            fitted <- stats::runmed(quo, k = 31, endrule = "keep")
            names(fitted) <- names(quo)
            list(fitted = fitted[nms], residuals = quo[nms] - 
                fitted[nms])
        }
        getQuality.SD <- function(allT) {
            sapply(allT, function(x) {
                median(abs(diff(unlist(lapply(x$lCTS, function(y) y$smoothed)))))
            })
        }   

##########################################      FILTERING      ########################################################


       nrecords  <-  sapply(allT,function(x) sum(unlist(lapply(x$lCTS,function(y) y$records))))
        ambiguous <- sapply(allS,function(x) x$ambiguous)
        doublet <- sapply(allS,function(x) if(!is.null(x$bestfit)) !x$bestfit$ambiguous else F)
        qualities <- getQuality.SD(allT)
        thresholdQual <- quantile(qualities,probs=1-probs)## removing 10% cells with lowest read counts
        keep <- qualities<=thresholdQual & log2(nrecords)>=thresholdNrec & !ambiguous & !doublet
        keep2 <- !(qualities<thresholdQual & log2(nrecords)<thresholdNrec)
        keep2 <- keep2 & !ambiguous
        ll <- getloess(qualities[keep2], log2(nrecords)[keep2])
        ord <- order(log2(nrecords)[keep2],decreasing=F)
       filters <- (1:length(nrecords)) %in% (which(keep2)[ll$residuals <= 
            distance_filter]) & keep

##########################################      PLOTTING      ########################################################
        pdf(paste0(outdir,projectname, "_qc.pdf"))
        plot(qualities,
             log2(nrecords),
             xlab="Noise logr",
             ylab="Total number of reads",
             pch=ifelse(doublet,15,19),
             cex=ifelse(doublet,.3,.1),
             col=ifelse(ambiguous,rgb(1,0,.5,.5),rgb(0,0,0,.5)))
        points(ll$fitted[ord],log2(nrecords)[keep2][ord],type="l",col=rgb(1,0,0,.5),lwd=1.5)
        abline(v=thresholdQual,h=thresholdNrec)
        points(qualities[keep2],
               log2(nrecords)[keep2],
               col=ifelse(abs(ll$residuals)>distance_filter, rgb(1,0,0,1), rgb(0,0,0,0)),cex=.15)
        try({plot(density(log2(nrecords[!doublet])))
            polygon(density(log2(nrecords[doublet])),border="red")
        },silent=T)
        plot(qualities,
             log2(nrecords),
             xlab="Noise logr",
             ylab="Total number of reads",
             pch=ifelse(filters,15,19),
             cex=ifelse(filters,.3,.1),
             col=ifelse(filters,rgb(0, 1, 0, 1),rgb(0,0,0,.5)))
        
            abline(v=thresholdQual,h=thresholdNrec)

        dev.off()


        
    })
res$filters <- filters
return(res)
}





res_update<-getFilters_update(res, probs = probs, outdir=outdir, projectname = project_name, thresholdNrec = nrec, distance_filter = dist)
length(res_update$filters)
#####
passing_cells<-names(res_update$filters[res_update$filters==TRUE])


outdir<-paste0(outdir,'/output')
dir.create(outdir)

printResults_all <- function(res,
                             ismedian=FALSE,
                             outdir=outdir,
                             projectname=project_name,
                             svinput=NULL,
                             lSVinput=NULL,
                             rainbowChr=TRUE)
{
    GAMMA <- 1
    if(any(names(res)=="gamma")) GAMMA <- res$gamma
    createDir <- paste0("mkdir ",outdir,"/profiles_",projectname)
    system(createDir)
    for(i in 1:length(res$allTracks.processed))
    {
        png(paste0(outdir,"/profiles_",projectname,"/",names(res$allTracks)[i],".png"), width = 5500, height = 2496, res=300)
        try({
            suppressWarnings(plotSolution(res$allTracks.processed[[i]],
                                          purity=res$allSolutions[[i]]$purity,
                                          ploidy=res$allSolutions[[i]]$ploidy,
                                          gamma=GAMMA,
                                          ismedian=ismedian,
                                          allchr=res$chr,
                                          ismale=if(!is.null(res$sex)) res$sex[i]=="male" else "female",
                                          sol=res$allSolutions[[i]],
                                          svinput=if(!is.null(lSVinput)) lSVinput[[i]] else svinput,
                                          rainbowChr=rainbowChr))
            title(names(res$allTracks)[i])
        })
        dev.off()
        try({
            writeProfile(prof=res$allProfiles[[i]],
                         samplename=paste0(names(res$allTracks)[i],"_",projectname),
                         outdir=outdir)
        })
    }

    createDirRefit <- paste0("mkdir ",outdir,"/profiles_",projectname,"_refitted")
    system(createDirRefit)
    if(any(grepl("refitted",names(res))))
    {
        for(i in 1:length(res$allTracks.processed))
        {
            png(paste0(outdir,"/profiles_",projectname,"_refitted/",names(res$allTracks)[i],".png"), width = 5500, height = 2496, res=300)
            try({
                suppressWarnings(plotSolution(res$allTracks.processed[[i]],
                                              purity=res$allSolutions.refitted.auto[[i]]$purity,
                                              ploidy=res$allSolutions.refitted.auto[[i]]$ploidy,
                                              ismale=if(!is.null(res$sex)) res$sex[i]=="male" else "female",
                                              allchr=res$chr,
                                              gamma=GAMMA,
                                              ismedian=ismedian,
                                              svinput=if(!is.null(lSVinput)) lSVinput[[i]] else svinput,
                                              sol=res$allSolutions[[i]],
                                              rainbowChr=rainbowChr))
                title(paste0(names(res$allTracks)[i],"-refitted"))
            })
            dev.off()
            try({
                writeProfile(prof=res$allProfiles.refitted.auto[[i]],
                             samplename=paste0(names(res$allTracks)[i],"_",projectname,"_refitted"),
                             outdir=outdir)
            })
        }
     
    }
    .mytry <- function(x,retVal=NA,...)
    {
        out <- try(x,silent=T,...)
        if(inherits(out,"try-error")) return(retVal)
        out
    }
    getploidy <- function(tt)
    {
        tt <- data.frame(chromosome=as.character(tt[,"chromosome"]),
                         start=as.numeric(tt[,"start"]),
                         end=as.numeric(tt[,"end"]),
                         total_copy_number=as.numeric(tt[,"total_copy_number"]))
        sizes <- (tt$end-tt$start)/1000000
        isna <- is.na(sizes) | is.na(tt$total_copy_number)
        sum(tt$total_copy_number[!isna]*sizes[!isna],na.rm=T)/sum(sizes[!isna],na.rm=T)
    }
    try({res <- append(res,
                       list(summary=list(allSols=data.frame(samplename=names(res$allTracks),
                                                            purity=sapply(res$allSolutions,function(x) .mytry(x$purity)),
                                                            ploidy=sapply(res$allSolutions,function(x) .mytry(x$ploidy)),
                                                            ploidy.tumour=sapply(res$allProfiles,function(x) .mytry(getploidy(x)))),
                                         allSols.refitted=if(!any(grepl("refitted",names(res)))) NULL
                                                          else
                                                              data.frame(samplename=names(res$allTracks),
                                                                         purity=sapply(res$allSolutions.refitted.auto,function(x) .mytry(x$purity)),
                                                                         ploidy=sapply(res$allSolutions.refitted.auto,function(x) .mytry(x$ploidy)),
                                                                         ploidy.tumour=sapply(res$allProfiles.refitted.auto,function(x) .mytry(getploidy(x)))))))})
    outdir <- gsub("/$","",outdir)
    try(write.table(res$summary$allSols,
                    file=paste0(outdir,
                                "/summary_",
                                projectname,".txt"),
                    sep="\t",quote=F,
                    col.names=T,row.names=T))
    try(write.table(res$summary$allSols.refitted,
                    file=paste0(outdir,"/summary_",projectname,"_refitted.txt"),
                    sep="\t",quote=F,col.names=T,row.names=T))
    save(res, file=paste0(outdir,"/result_object_",projectname,".Rda"))
}

printResults_all(res_update,outdir=outdir)
