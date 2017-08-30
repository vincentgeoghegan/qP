#qP 0.86
#Vincent Geoghegan, Lancaster University, 2016
#to do
#calculate CV of normalised DNA
#make initialDNA work if there is only one technical replicate
#add ability to take info from other columns for annotation later
#updata calculateratios to take list
#sort out column names
#ability to calculate ratios by sample name or treatment
#normalise data without having to specify sample groups

library("ggplot2")
library("reshape2")
library("boot")

# Multiple plot function
# from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


ReadData <- function(file_name, reference_gene, target_gene, experiment){
  #processes input data, extracts names of all samples into a vector
  #extracts names of dilution series into a vector
  data <- read.csv(file_name, header=TRUE, stringsAsFactors=FALSE)
  #get names of samples
  all_names <- unique(data$Biological.Set.Name[data$Biological.Set.Name !=""])
  #extract names of dilution series
  targets <- as.character(unique(data$Target[data$Target != ""]))
  #generate output file name
  out_fname <- gsub('.csv','_qP.csv',file_name)
  
  #get columns other than target, Sample and cq for reconstructing sample annotation later
  descriptors <- data[, -which(colnames(data) %in% c('Target','Cq', 'Sample'))]
  
  Reference <- reference_gene
  Target <- target_gene
  Experiment <- experiment
  
  meta_data <- data.frame(Reference, Target, Experiment, stringsAsFactors=FALSE)
  
  prepared_data <- list(out_fname=out_fname, data=data, all_names=all_names, targets=targets, 
                        metadata=meta_data, descriptors=descriptors)
  
  return(prepared_data)
  
}

PrimerEff <- function(data_list){
  #calculates primer efficiencies based on dilution series
  #outputs graphs of primer efficiencies
  
  #extract data frame containing all raw data
  raw_data <- data_list$data
  #TO DO: remove dilution series if there are <2 data points
  #loop over dilutions calculating efficiencies
  calculate_efficiencies <- function(x){
    
    search_string <- as.character(x)
    dilution_search <- paste("^", "d", sep="")
    #extract data for the dilution series for this target
    dilution_set <- raw_data[which(grepl(search_string, raw_data$Target)==TRUE & grepl(dilution_search, raw_data$Biological.Set.Name)==TRUE),]
    #have to make sample column numeric because if there were sample names elsewhere will be treated as factor
    dilution_set$log2.copy <- log2(as.numeric(dilution_set$Sample))
    #generate linear model
    dilution_set_lm <- lm(Cq ~ log2.copy, data=dilution_set)
    #get slope of line
    slope <- coef(dilution_set_lm)[[2]]
    #calculate primer efficiency
    primer_efficiency <- 2^(-1/slope)
    #set the position for efficiency label
    label_x <- mean(dilution_set$log2.copy, na.rm=TRUE)
    label_y <- mean(dilution_set$Cq, na.rm=TRUE)
    
    #test how many data points we have for the primer, if there are
    #<2 data points for standard curve we cannot plot anything, so output dummy graph
    if(length(dilution_set$log2.copy) < 2 | length(dilution_set$Cq) < 2){
      dummy_df <- data.frame(x=c(1,2,3,4,5), y=c(1,2,3,4,5))
      eff_plot <- ggplot(data=dummy_df, mapping=aes(x=x, y=y)) +
        geom_blank() +
        ggtitle(paste(dilution_set$Target[1], "Error: insufficient data points"))
      
    }else{
      #sufficient data points so plot dilution curve
      
      eff_plot <- ggplot(dilution_set, aes(x=log2.copy, y=Cq)) +
        geom_point(shape=19, size=3) +    # Use hollow circles
        geom_smooth(method=lm)+
        annotate("text", x = label_x, y = label_y, label = primer_efficiency)+
        ggtitle(paste(dilution_set$Target[1], "efficiency"))
    }
    
    #store primer efficiency and plot in list
    return(list(primer_efficiency, eff_plot, dilution_set_lm))
    
  }
  
  #make sure each list element is named after the sample
  effs_plots <- setNames(lapply(data_list$targets, function(x) calculate_efficiencies(x)),data_list$targets)
  
  
  #extract every second element (the plots)
  plots <- lapply(effs_plots, "[[", 2)
  
  #to do: fix outputting to pdf so that it works for lots of graphs
  
  #write plots
  pdf("Efficiencies.pdf")
  for(i in 1:length(plots)){
    print(plots[i])
  }
  #multiplot(plotlist=plots, cols=2)
  dev.off()
  return(effs_plots)
}


#get all control and treated samples, ie everything except dilutions
#sample_names <- as.character(unique(comparisons_df$sample_names))



InitialDNA <- function(prepared_data, effs_plots){
  #calculates relative initial amount of DNA in reaction

  #extract the original data
  raw_data <- prepared_data$data
  
  metadata <- prepared_data$metadata
  effs_plots <- effs_plots
  
  calculate_initialDNA <- function(x){
    #extract data for the sample
    this_sample <- raw_data[which(raw_data$Biological.Set.Name==x),]
    
    #get target name for sample
    target_name <- as.character(metadata$Target)
    #get reference name for sample
    ref_name <- as.character(metadata$Reference)
    
    t_names <- c(target_name, ref_name)
    print(t_names)
    
    getInitial_targ_ref <- function(x){
      
      #target reactions
      targets_only <- this_sample[which(this_sample$Target==x),]
      
      #extract the efficiencies
      efficiency_list <- lapply(effs_plots, "[[", 1)
      #get the efficiency for the primers
      targets_only$Efficiency <- as.numeric(efficiency_list[x][[1]])
      #get rows with non-NA Cq values
      targets_only_nonNA <- targets_only[which(is.na(targets_only$Cq)==FALSE),]
      #get rows with NA Cq values
      targets_only_NA <- targets_only[which(is.na(targets_only$Cq)==TRUE),]
      #found NA Cq values
      if(nrow(targets_only_NA)!=0){
        #set initial DNA to NA
        targets_only_NA$Initial <- NA
        targets_only_NA$CV <- NA
      }
      
      #calculate initial DNA for each technical replicate so coefficient of variance can be calculated
      #initial amount of DNA is inversely proportional to E^Cq value
      #so divide 1 by E^Cq, initial amount of DNA is directly proportional to this
      #only pass in rows with non-NA Cq values
      targets_only_nonNA$Initial <- 1/(targets_only_nonNA$Efficiency^(as.numeric(targets_only_nonNA$Cq)))
      
      #**********************************
      
      #calculate coefficient of variance for initial DNA
      #calculate standard deviation first
      initialDNA_sd <- sd(targets_only_nonNA$Initial, na.rm=TRUE)
      #calculate mean of initial DNA
      initialDNA_mean <- mean(targets_only_nonNA$Initial, na.rm=TRUE)
      #calculate coefficient of variance
      initialDNA_cv <- (initialDNA_sd/initialDNA_mean)*100
      #if there is only one technical replicate, CV will come out as NA
      #so set it to 1
      
      #make cv column
      targets_only_nonNA$CV <- initialDNA_cv
      
      #**********************
      #combine non-NA and NA rows again
      targets_only <- rbind(targets_only_nonNA, targets_only_NA)
      
      #re-calculate initial DNA, this time using means of Cq values of technical replicates
      #aggregate data on mean of Cq
      target_sub_mean <- aggregate(Cq ~ Biological.Set.Name + Target + Efficiency + CV, 
                                   data=targets_only[,c('Cq','Biological.Set.Name','Target','Efficiency','CV')],
                                   mean, na.rm=TRUE, na.action=na.pass)
      
      #initial amount of DNA is inversely proportional to E^Cq value
      #so divide 1 by E^Cq, initial amount of DNA is directly proportional to this
      target_sub_mean$Initial <- 1/(target_sub_mean$Efficiency^(as.numeric(target_sub_mean$Cq)))
      #if length of first element in efficiency and plot list is <3 then efficiencies have been supplied by user
      #and data is assumed to be relative quantification so no need to calculate copies
      if(length(effs_plots[[1]]) < 3){
        target_sub_mean$Copies <- 1
      #efficiencies have been calculated from dilution series
      }else{
        #extract linear models
        lms <- lapply(effs_plots, "[[", 3)
        #calculate copies, for absolute quantification data, will not have meaning for relative data
        #extract parameters for line
        intercept <- coef(lms[[x]])[1]
        gradient <- coef(lms[[x]])[2]
        print(this_sample$Biological.Set.Name)
        target_sub_mean$Copies <- 2^((target_sub_mean$Cq - intercept)/gradient)
        
      }
      
      return(target_sub_mean)
    }
    id <- setNames(lapply(t_names, function(x) getInitial_targ_ref(x)), t_names)
    print(id)
    return(id)
    
  }
  all_names <- prepared_data$all_names
  #get all sample names, ie everything that isn't the dilution series which begins with 'd'
  sample_names <- all_names[!grepl("^d", all_names)]
  initialDNA <- setNames(lapply(sample_names, function(x) calculate_initialDNA(x)), sample_names)
  #add file name to list
  initialDNA$out_fname <- prepared_data$out_fname
  #add metadata
  initialDNA$metadata <- metadata
  #add descriptors
  initialDNA$descriptors <- prepared_data$descriptors
  return(initialDNA)
}


NormaliseDNA <- function(initialDNA, control_name='control', control_samples=NULL, treated_name='treated', treated_samples=NULL){
  #TO DO calculate CV of normalised DNA
  
  #generate vector of control and treated sample names
  control <- strsplit(control_samples,';')[[1]]
  treated <- strsplit(treated_samples,';')[[1]]
  
  metadata <- initialDNA$metadata
  
  #normalise amount of DNA to reference
  normalise <- function(x){
    #get initial amount of DNA for target
    target_initial <- x[[1]]$Initial
    
    #get initial amount of DNA for ref
    ref_initial <- x[[2]]$Initial
    #normalise
    x[[1]]$Normalised <- target_initial/ref_initial
    
    #get treatment name
    if(is.element(x[[1]]$Biological.Set.Name,control)==TRUE){
      treatment_name <- control_name
    }else{
      treatment_name <- treated_name
    }
    
    #calculate target per reference, only makes sense with absolute data
    target_copies <- x[[1]]$Copies
    target_name <- names(x[1])
    ref_copies <- x[[2]]$Copies
    ref_name <- names(x[2])
    #use target and ref name to name the absolute quantification column
    absolute_name <- paste(target_name, ref_name, sep='.per.')
    x[[1]][,paste(absolute_name)] <- target_copies/ref_copies
    
    #make second data frame same size
    x[[2]]$Normalised <- "NA"
    x[[2]][,paste(absolute_name)] <- "NA"
    
    #put in treatment name
    x[[1]]$Treatment <- treatment_name
    x[[2]]$Treatment <- treatment_name
    
    return(x)    
  }
  #output file name
  out_file_name <- initialDNA$out_fname
  len_id <- length(initialDNA)
  #pass in initial DNA list without elements corresponding to file name, metadata and descriptors
  normalisedDNA <- lapply(initialDNA[-which(names(initialDNA) %in% c('out_fname', 'metadata', 'descriptors'))], function(x) normalise(x))
  
  #unlist the raw data to write results as a table
  unlist_results <- function(x){
    #unlist the first target
    element1 <- unlist(x[1])
    #unlist the second target
    element2 <- unlist(x[2])
    #remove first two rows which just describe the sample name and target name
    #this information will be present as column names and row names anyway
    element <- c(element1[c(-1,-2)], element2[c(-1,-2)])
    return(element)
  }
  #unlist raw data, returning a matrix
  unlisted <- sapply(normalisedDNA, function(x) unlist_results(x))
  
  
  #unlist the raw data to write results as a table in 'long' format
  unlist_results_long <- function(x){
    #unlist the first target
    element1 <- unlist(x[1])
    #unlist the second target
    element2 <- unlist(x[2])
    #remove first two rows which just describe the sample name and target name
    #this information will be present as column names and row names anyway
    #also remove last row which is treatment name
    element <- c(element1[c(-1,-2,-length(element1))], element2[c(-1,-2,-length(element2))])
    return(element)
  }
  unlisted_long <- sapply(normalisedDNA, function(x) unlist_results_long(x))
  #write out list of raw results to file
  out_file <- file(out_file_name, open="w")  #creates a file in append mode
  #write out data, col.names=NA will include an empty cell for row names header
  write.table(unlisted, file=out_file, sep=",", dec=".", quote=FALSE, col.names=NA)
  
  close(out_file)  #close connection to file.csv
  #write.csv(output, "qPCR results summary.csv", row.names=TRUE)
  #make data frame
  nd <- data.frame(unlisted_long)
  nd$Parameters <- rownames(nd)
  #melt data frame to long format
  nd_long <- melt(nd, id.vars='Parameters')
  #get the treatment name for each sample
  annotate_nd_long <- function(x){
    #get treatment name
    if(is.element(x,control)==TRUE){
      treatment_name <- control_name
    }else{
      treatment_name <- treated_name
    }
    return(treatment_name)
  }
  #get name of target only
  #note '.' means match any character in regex so make it a character class
  #so it is interpreted literally
  targets_col <- sapply(strsplit(nd_long[,1],'[.]'),'[[',1)
  #get name of parameter
  #param_col <- sapply(strsplit(nd_long[,1],'[.]'),'[[',2)
  param_col <- sapply(nd_long[,1], 
                      function(x){ifelse(grepl('[.]per[.]', x),'Absolute',strsplit(x, '[.]')[[1]][2])})
  #get name of experiment
  experiment <- metadata$Experiment
  #generate column describing the experiment
  exp_col <- rep(experiment, length(param_col))
  #get treatment name for each sample
  treatments <- sapply(nd_long[,2], function(x) annotate_nd_long(x))
  #bind all new data
  nd_long <- cbind(nd_long,targets_col)
  nd_long <- cbind(nd_long,param_col)
  nd_long <- cbind(nd_long, experiment)
  ann_nd_long <- cbind(nd_long, treatments)
  #rename columns
  colnames(ann_nd_long) <- c('Description','Sample','Value','Target','Parameter','Experiment','Treatment')
  ann_nd_long_fname <- gsub('_qP.csv','_qP_l.csv',out_file_name)
  write.table(ann_nd_long, file=ann_nd_long_fname, sep=",", dec=".", quote=FALSE, row.names=FALSE)
  #return the long data frame and data frame containing all information to label data
  normalised_data <- list(normalised=ann_nd_long, descriptors=initialDNA$descriptors, outfile=ann_nd_long_fname)
  return(normalised_data)
  
}

#annotates data with sample descriptors if they have been supplied, e.g. timepoint, treatment, concentration
AnnotateData <- function(normalised_data){
  normalised_df <- normalised_data$normalised
  descriptors <- normalised_data$descriptors
  #get the corresponding descriptors for each row
  get_descriptors <- function(x){
    #get rid of X in front of sample names
    s_name <- substr(as.character(x),2,nchar(as.character(x)))
    desc <- descriptors[which(descriptors$Biological.Set.Name == s_name),]
    return(desc[1,])
  }
  desc <- lapply(normalised_df$Sample, function(x) get_descriptors(x))
  #turn list into columns of data frame
  desc_df <- as.data.frame(t(sapply(desc, function(x) return(x))))
  #add columns of data frame to normalised data frame
  annotated_data <- cbind(normalised_df,desc_df)
  ad_fname <- gsub('_qP_l.csv','_qP_la.csv', normalised_data$outfile)
  print(ad_fname)
  write.table(as.matrix(annotated_data), file=ad_fname, sep=",", dec=".", quote=FALSE, row.names=FALSE)
  return(annotated_data)
  
}

CalculateRatios <- function(normalised, target, control, treated, type='Relative', CI_type='Boot', p_value='Wilcox',ByTreatment=TRUE){
  #calculates ratios between normalised DNA in control and treated samples
  #to specify individual samples to be compared set ByTreatment to FALSE#
  #and list samples as 's1;s2;s3' and 's4;s5;s6' for example
  if(ByTreatment==TRUE){
    if(type=='Relative'){
      #control data
      data_normalised_c <- normalised[which(normalised$Target==target & normalised$Parameter=='Normalised' 
                                            & normalised$Treatment==control),]
      #treated data
      data_normalised_t <- normalised[which(normalised$Target==target & normalised$Parameter=='Normalised' 
                                            & normalised$Treatment==treated),]
      
      #ratio <- mean(as.numeric(data_normalised_t$Value))/mean(as.numeric(data_normalised_c$Value))
      #data is absolute quantification
    }else{
      #control data
      data_normalised_c <- normalised[which(normalised$Target==target & normalised$Parameter=='Absolute' 
                                            & normalised$Treatment==control),]
      #treated data
      data_normalised_t <- normalised[which(normalised$Target==target & normalised$Parameter=='Absolute' 
                                            & normalised$Treatment==treated),]
      
      #ratio <- mean(as.numeric(data_normalised_t$Value))/mean(as.numeric(data_normalised_c$Value))
    }
    #individual samples have been specified 
  }else{
    control_samples <- strsplit(control,';')[[1]]
    treated_samples <- strsplit(treated,';')[[1]]
    if(type=='Relative'){
      data_normalised_c <- normalised[which(is.element(normalised$Sample, control_samples)==TRUE & normalised$Parameter=='Normalised' 
                                            & normalised$Target==target),]
      data_normalised_t <- normalised[which(is.element(normalised$Sample, treated_samples)==TRUE & normalised$Parameter=='Normalised' 
                                            & normalised$Target==target),]
      #ratio <- mean(as.numeric(data_normalised_t$Value))/mean(as.numeric(data_normalised_c$Value))
      
      #data is absolute quantification
    }else{
      data_normalised_c <- normalised[which(is.element(normalised$Sample, control_samples)==TRUE & normalised$Parameter=='Absolute' 
                                            & normalised$Target==target),]
      data_normalised_t <- normalised[which(is.element(normalised$Sample, treated_samples)==TRUE & normalised$Parameter=='Absolute' 
                                            & normalised$Target==target),]
      #ratio <- mean(as.numeric(data_normalised_t$Value))/mean(as.numeric(data_normalised_c$Value))
    }
  }
  ratio_matrix <- as.vector(outer(as.numeric(data_normalised_t$Value), as.numeric(data_normalised_c$Value), "/"))
  mean_t <- mean(as.numeric(data_normalised_t$Value), na.rm=TRUE)
  mean_c <- mean(as.numeric(data_normalised_c$Value), na.rm=TRUE)
  
  ratio_mean <- mean_t/mean_c
  #calculate confidence intervals
  #bootstrapping
  boot_CI <- function(treated, control){
    b_t <- boot(treated, function(u,i) mean(u[i]), R = 1000)
    boot_result <- boot.ci(b_t, type = c("perc"))
    upper_lower <- tail(as.vector(boot_result$percent),2)
    CI_lower_t <- upper_lower[1]
    CI_upper_t <- upper_lower[2]
    
    b_c <- boot(control, function(u,i) mean(u[i]), R = 1000)
    boot_result <- boot.ci(b_c, type = c("perc"))
    upper_lower <- tail(as.vector(boot_result$percent),2)
    CI_lower_c <- upper_lower[1]
    CI_upper_c <- upper_lower[2]
    return(c(CI_lower_t/CI_lower_c, CI_upper_t/CI_upper_c))
    
  }
  if(CI_type=='Boot'){
    CI <- boot_CI(as.numeric(data_normalised_t$Value), as.numeric(data_normalised_c$Value))
  }
  if(p_value=='Wilcox'){
    w_test <- wilcox.test(as.numeric(data_normalised_t$Value), as.numeric(data_normalised_c$Value))
    p_value <- w_test['p.value']
  }
  return(c(ratio_mean, CI, p_value))
  
}





