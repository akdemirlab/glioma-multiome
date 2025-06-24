
files_1 <- list.files(path = '***Datapath***', full.names = TRUE)
bams <- files_1[!grepl(".bai", basename(files_1))]
bams <- bams[!grepl("Ctrl",bams)]



#intervals from AmpliconArchitect. I used the largest cycle to define the ecDNA region.
# # Interval        1       chr7    53444680        55474670

svinput=data.frame(chr=c("chr7","chr7"),
                   pos=c(53444680,55474670))


extracted_TIMEPOINT <- gsub(".*SC521_(.*)_W.*", "\\1", bams)

#fit ploidy constraints based on Dapi FACS plots.
extracted_TIMEPOINT <- gsub(".*SC521_(.*)_W.*", "\\1", bams)
ploidy_list<-list()
i=0
for(i in 1:length(extracted_TIMEPOINT)){
        ploidy_vector<-c()
        cell_ploidy<-extracted_TIMEPOINT[i]
        if (cell_ploidy=="1"){
            ploidy_vector <- seq(3.5,6, 0.01)
            ploidy_list[[i]]<-ploidy_vector
        }
        else if (cell_ploidy=="2"){
            ploidy_vector <- seq(3.5,6, 0.01)
            ploidy_list[[i]]<-ploidy_vector
        }           
        i<-i+1
        }

print(i)
ploidy_list

length(extracted_TIMEPOINT)

#purity
purity_list<-rep(list(seq(.999, 1, 0.001)), length(bams))

print('*****************************************************************************************************')
run_sc_sequencing(tumour_bams=bams, ##vector of full paths to your bam files
                          allchr=paste0("chr",c(1:22,"X", "Y")), 
                          sex=rep("male",length(bams)), ##sex for each bam
                          binsize=500000, ##bin size - reduce at higher depths
                          chrstring_bam="chr",
                          purs = purity_list,##value of purity
                          ploidies = ploidy_list, ##value of ploidy
                          maxtumourpsi=6, ## maximum tumour poidy considered
                          build="hg38", ##build, either hg19 or hg38
                          probs_filters=0.1,
                          sc_filters=TRUE,
                          MC.CORES=64,##number of cores
                          projectname=project_name,
                          segmentation_alpha=0.1,
                          predict_refit=TRUE,
                          outdir=output,
                          smooth_sc=FALSE,
                          print_results=TRUE,
                          multipcf=TRUE,
                          svinput=svinput)



