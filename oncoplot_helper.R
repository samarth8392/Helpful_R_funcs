mutation_colors <- c(
  "Nonsense mutation" = "#ad7aff",
  "Nonsense_Mutation" = "#ad7aff", 
  "Missense mutation" = "#377EB8",
  "Missense_Mutation" = "#377EB8",
  "Frameshift deletion" = "#4DAF4A",
  "Frame_Shift_Del" = "#4DAF4A",
  "In Frame Ins" = "#ff008c",
  "In_Frame_Ins" = "#ff008c",
  "Splice Site" = "#FF7F00", 
  "Splice_Site" = "#FF7F00",
  "Multi Hit" = "#FFFF33",
  "Multi_Hit" = "#FFFF33",
  "Frameshift insertion" = "#A65628",
  "Frame_Shift_Ins" = "#A65628",
  "In Frame Del" = "#f781bf",
  "In_Frame_Del" = "#f781bf",
  "Translation Start site" = "#400085",
  "Translation_Start_Site" = "#400085",
  "Nonstop Mutation" = "#b68dfc",
  "Nonstop_Mutation" = "#b68dfc",
  "Amplification" = "red",
  "Deletion" = "blue",
  "Splice Region" = "red",
  "Splice_Region" = "red",
  "Silent" = "blue",
  "Germline Missense Mutation" = "#377EB8",
  "Germline_Missense_Mutation" = "#377EB8",
  "Germline Frame Shift Del" = "#4DAF4A",
  "Germline_Frame_Shift_Del" = "#4DAF4A",
  "Germline Frameshift Del" = "#4DAF4A",
  "Germline_Frameshift_Del" = "#4DAF4A",
  "Germline Frameshift deletion" = "#4DAF4A",
  "Germline_Frameshift_Deletion" = "#4DAF4A",
  "Germline Nonsense Mutation" = "#ad7aff",
  "Germline_Nonsense_Mutation" = "#ad7aff",
  "Pathogenic" = "black",
  "COSMIC" = "grey25",
  "VUS" = "grey50"
)

### List defining functions for color and shape of cells in oncoplot
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # "0" = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
  #             gp = gpar(fill = "#CCCCCC", col = NA))
  # },
  "Nonsense Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = mutation_colors["Nonsense mutation"], col = NA))
  },
  "Nonsense_Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = mutation_colors["Nonsense mutation"], col = NA))
  },
  "Missense Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = mutation_colors["Missense mutation"], col = NA))
  },
  "Missense_Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = mutation_colors["Missense mutation"], col = NA))
  },
  "Frame_Shift_Del" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = mutation_colors["Frameshift deletion"], col = NA))
  },
  "Frame Shift Del" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = mutation_colors["Frameshift deletion"], col = NA))
  },
  "In Frame Ins" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = mutation_colors["In frame insertion"], col = NA))
  },
  "Splice_Site" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Splice site"], col = NA))
  },
  "Splice Site" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Splice site"], col = NA))
  },
  "Multi Hit" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Multi hit mutation"], col = NA))
  },
  "Multi_Hit" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Multi hit mutation"], col = NA))
  },
  "Frame Shift Ins" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Frameshift insertion"], col = NA))
  },
  "In Frame Del" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = mutation_colors["In frame deletion"], col = NA))
  },
  "Nonstop Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Nonstop mutation"], col = NA))
  },
  "Nonstop_Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Nonstop mutation"], col = NA))
  },
  "Translation Start Site" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Translation Start site"], col = NA))
  },
  "Amp" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.25, "mm"), h-unit(0.25, "mm"),
              gp = gpar(fill = mutation_colors["Amp"], col = NA))
  },
  "Del" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.25, "mm"), h-unit(0.25, "mm"),
              gp = gpar(fill = mutation_colors["Del"], col = NA))
  },
  ##Silent mutations
  "Silent" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Silent"], col = NA))
  },
  "Splice Region" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Splice region"], col = NA))
  },
  
  "no variants" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              # gp = gpar(fill = "#e0e0e0", col = NA))
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  "Pathogenic" = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.8,
              gp = gpar(col = mutation_colors["Pathogenic"], fill = NA, lwd=5))
  },
  "COSMIC" = function(x, y, w, h) {
   # grid.rect(x, y, w*0.8, h*0.8,
  #            gp = gpar(col = mutation_colors["COSMIC"], fill = NA, lwd=5))
    grid.points(x, y,  pch = 16, size=w*.5, 
                gp = gpar(fill = mutation_colors["COSMIC"],col=mutation_colors["COSMIC"]))
  },
  "VUS" = function(x, y, w, h) {
    # grid.points(x, y, pch = 3, size=w,gp=gpar(col=col["VUS"], lwd=3))
    # grid.rect(x, y, w*0.2, h-unit(0.5, "mm"),
    #           gp = gpar(fill = col["VUS"], col = NA))
    grid.rect(x, y, w*0.6, h*0.8,
              gp = gpar(col = mutation_colors["VUS"], fill = NA, lwd=5))
  },
  "Germline Missense Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w*0.6, h*0.3,
              gp = gpar(fill = mutation_colors["Germline missense mutation"]))
  },
  "Germline Frame Shift Del" = function(x, y, w, h) {
    grid.rect(x, y, w*0.6, h*0.3,
              gp = gpar(fill = mutation_colors["Germline frameshift deletion"]))
  },
  "Germline Nonsense Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w*0.6, h*0.3,
              gp = gpar(fill = mutation_colors["Germline nonsense mutation"]))
  }
)




createMatCHM=function(m, g = NULL, chatty = TRUE, add_missing = FALSE,COSMIC=FALSE,CLINVAR=FALSE,includeSilent=F){
  #browser()
  if(is.null(g)){
    stop("Please provde atleast two genes!")
  }
  
  subMaf = subsetMaf(maf = m, genes = g, includeSyn = FALSE, mafObj = FALSE,dropLevels = F)
  
  if(includeSilent){
    subMaf=m@maf.silent %>% filter(Hugo_Symbol %in% g) %>% bind_rows(subMaf,.)
    
  }
  t1=subMaf %>% dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification) %>%
    group_by(Hugo_Symbol, Tumor_Sample_Barcode) %>%
    mutate(oncomat_label1=paste0(ifelse(length(unique(as.character(Variant_Classification)))>1, "Multi_Hit",unique(as.character(Variant_Classification))))) %>%
    distinct(Hugo_Symbol,Tumor_Sample_Barcode,oncomat_label1,.keep_all = T) #%>%
  #mutate(oncomat_label1=gsub("_"," ",oncomat_label1),Variant_Classification=gsub("_"," ",Variant_Classification))
  
  if(CLINVAR){
    if("CLIN_SIG" %in%names(subMaf)){
      t2=subMaf %>% 
        dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode, CLIN_SIG) %>%
        group_by(Hugo_Symbol, Tumor_Sample_Barcode) %>%
        mutate(oncomat_label2=ifelse(grepl("pathogenic",CLIN_SIG)|grepl("Pathogenic",CLIN_SIG)|grepl("conflicting",CLIN_SIG),"Pathogenic",
                                     ifelse(grepl("uncertain",CLIN_SIG)|grepl("Uncertain",CLIN_SIG), "VUS",''))) %>%
        filter(oncomat_label2!="") %>%
        distinct(Hugo_Symbol,Tumor_Sample_Barcode,oncomat_label2,.keep_all = T) %>%
        filter(!(n()>1 & oncomat_label2 %in% "VUS"))
      
    }else{
      t2=subMaf %>% dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode, CLIN_SIG=CLNSIG) %>%
        group_by(Hugo_Symbol, Tumor_Sample_Barcode) %>%
        mutate(oncomat_label2=ifelse(grepl("pathogenic",CLIN_SIG)|grepl("Pathogenic",CLIN_SIG)|grepl("conflicting",CLIN_SIG),"Pathogenic",
                                     ifelse(grepl("uncertain",CLIN_SIG)|grepl("Uncertain",CLIN_SIG), "VUS",''))) %>%
        filter(oncomat_label2!="") %>%
        distinct(Hugo_Symbol,Tumor_Sample_Barcode,oncomat_label2,.keep_all = T) %>%
        arrange(Tumor_Sample_Barcode,oncomat_label2) %>%
        filter(!(n()>1 & oncomat_label2 %in% "VUS"))
    }
    
  }
  if(COSMIC){
    t3=subMaf %>% dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode, cosmic92_noncoding,cosmic92_coding) %>%
      mutate(COSMIC=ifelse(grepl("ID=",cosmic92_coding)|grepl("ID=",cosmic92_noncoding),"COSMIC",NA)) %>% 
      dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode, COSMIC) %>%
      group_by(Hugo_Symbol, Tumor_Sample_Barcode)  %>% 
      mutate(oncomat_label3=ifelse(grepl("COSMIC",COSMIC),"COSMIC","")) %>%
      filter(oncomat_label3!="") %>%
      distinct(Hugo_Symbol,Tumor_Sample_Barcode,oncomat_label3,.keep_all = T)
    
  }
  ##Try to create the final ONCOMAT!
  if(exists('t3')){
    oncomat1=full_join(t1,t2,by=c("Hugo_Symbol",'Tumor_Sample_Barcode')) %>% 
      full_join(.,t3,by=c("Hugo_Symbol",'Tumor_Sample_Barcode')) %>% 
      unite(oncomat_label,starts_with("oncomat_label"),sep=";",na.rm=T) %>%
      dplyr::select(Hugo_Symbol,Tumor_Sample_Barcode,oncomat_label) %>%
      ungroup() %>%
      pivot_wider(names_from=Tumor_Sample_Barcode,values_from=oncomat_label,values_fill = "")
  }else{
    oncomat1=full_join(t1,t2,by=c("Hugo_Symbol",'Tumor_Sample_Barcode')) %>% 
      #full_join(.,t3,by=c("Hugo_Symbol",'Tumor_Sample_Barcode')) %>% 
      unite(oncomat_label,starts_with("oncomat_label"),sep=";",na.rm=T) %>%
      dplyr::select(Hugo_Symbol,Tumor_Sample_Barcode,oncomat_label) %>%
      ungroup() %>%
      pivot_wider(names_from=Tumor_Sample_Barcode,values_from=oncomat_label,values_fill = "")
  }
  
  
  if(!all(unique(m@clinical.data$Tumor_Sample_Barcode) %in% colnames(oncomat1)[-1])){
    mis1=m@clinical.data$Tumor_Sample_Barcode[!m@clinical.data$Tumor_Sample_Barcode %in% colnames(oncomat1)[-1]]
    mis1=paste(mis1,sep=",")
    oncomat1[,mis1]=""
  }
  
  
  #convert to matrix
  #data.table::setDF(oncomat)
  oncomat=data.frame(oncomat1)
  rownames(oncomat) = oncomat$Hugo_Symbol
  oncomat = as.matrix(oncomat[,-1, drop = FALSE])
  colnames(oncomat)=colnames(oncomat1)[-1]
  
  variant.classes = as.character(unique(subMaf[,Variant_Classification]))
  variant.classes = c('',variant.classes, 'Multi_Hit')
  names(variant.classes) = 0:(length(variant.classes)-1)
  
  #Complex variant classes will be assigned a single integer.
  vc.onc = unique(unlist(apply(oncomat, 2, unique)))
  vc.onc = vc.onc[!vc.onc %in% names(variant.classes)]
  names(vc.onc) = rep(as.character(as.numeric(names(variant.classes)[length(variant.classes)])+1), length(vc.onc))
  variant.classes2 = c(variant.classes, vc.onc)
  
  oncomat.copy <- oncomat
  #Make a numeric coded matrix
  for(i in 1:length(variant.classes2)){
    oncomat[oncomat == variant.classes2[i]] = names(variant.classes2)[i]
  }
  
  #If maf has only one gene
  if(nrow(oncomat) == 1){
    mdf  = t(matrix(as.numeric(oncomat)))
    rownames(mdf) = rownames(oncomat)
    colnames(mdf) = colnames(oncomat)
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }
  
  #convert from character to numeric
  mdf = as.matrix(apply(oncomat, 2, function(x) as.numeric(as.character(x))))
  rownames(mdf) = rownames(oncomat.copy)
  
  
  #If MAF file contains a single sample, simple sorting is enuf.
  if(ncol(mdf) == 1){
    sampleId = colnames(mdf)
    mdf = as.matrix(mdf[order(mdf, decreasing = TRUE),])
    colnames(mdf) = sampleId
    
    oncomat.copy = as.matrix(oncomat.copy[rownames(mdf),])
    colnames(oncomat.copy) = sampleId
    
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  } else{
    #Sort by rows as well columns if >1 samples present in MAF
    #Add total variants per gene
    mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
      length(x[x != "0"])
    }))
    #Sort by total variants
    mdf = mdf[order(mdf[, ncol(mdf)], decreasing = TRUE), ]
    #colnames(mdf) = gsub(pattern = "^X", replacement = "", colnames(mdf))
    nMut = mdf[, ncol(mdf)]
    
    mdf = mdf[, -ncol(mdf)]
    
    mdf.temp.copy = mdf #temp copy of original unsorted numeric coded matrix
    
    mdf[mdf != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
    tmdf = t(mdf) #transposematrix
    mdf = t(tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), decreasing = TRUE)), ]) #sort
    
    mdf.temp.copy = mdf.temp.copy[rownames(mdf),] #organise original matrix into sorted matrix
    mdf.temp.copy = mdf.temp.copy[,colnames(mdf)]
    mdf = mdf.temp.copy
    
    #organise original character matrix into sorted matrix
    oncomat.copy <- oncomat.copy[,colnames(mdf)]
    oncomat.copy <- oncomat.copy[rownames(mdf),]
    
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }
}



renamefxn=function(x){
  
  x1=ifelse(grepl("frame shift",x,ignore.case = T),gsub("frame shift","Frameshift",x,ignore.case = T),x)
  
  
  x1=ifelse(grepl("\\bDel\\b",x1,ignore.case = T),gsub("\\bDel\\b","deletion",x1,ignore.case = T),x1)
 # x1=ifelse(grepl("\\bDel",x1,ignore.case = T),sub("\\bDel","deletion",x1,ignore.case = T),x1)
  
  x1=ifelse(grepl("\\bIns\\b",x1,ignore.case = T),gsub("\\bIns\\b","insertion",x1,ignore.case = T),x1)
#  x1=ifelse(grepl("\\bIns",x1,ignore.case = T),sub("\\bIns","insertion",x1,ignore.case = T),x1)
  
  x1=ifelse(grepl("In Frame",x1,ignore.case = T),gsub("In Frame","In frame",x1,ignore.case = T),x1)
  
  x1=ifelse(grepl("Multi hit",x1,ignore.case = T),gsub("Multi hit","Multi hit mutation",x1,ignore.case = T),x1)
  x1=ifelse(grepl("Mutation",x1,ignore.case = T),gsub("Mutation","mutation",x1,ignore.case = T),x1)
  x1=ifelse(grepl("Site",x1,ignore.case = T),gsub("Site","site",x1,ignore.case = T),x1)
  x1=ifelse(grepl("Vus",x1,ignore.case = T),gsub("Vus","VUS",x1,ignore.case = T),x1)
    #x1=ifelse(x1 %in% c("VUS","COSMIC"),x1,str_to_sentence(x1,))
  
  
  return(x1)
}




renameoncomat=function(dt){
  require(stringi)
  if(is.matrix(dt)){
    tempname=colnames(dt)
    tempdt=data.frame(dt)
    tempdt1=mutate_all(tempdt,renamefxn)
    tempdt1=as.matrix(tempdt1)
    colnames(tempdt1)=tempname
    tempdt1
    }else if(is.data.frame(dt)){
    tempdt1=mutate_all(dt,~renamefxn(.x))
    tempdt1
    }else{
      renamefxn(dt)
  }
  
  
}


createoncomat=function(maf,genes,CLINVAR=F,COSMIC=F,add_missing=F){
  oncomat <- createMatCHM(maf, g=genes, add_missing = T,COSMIC,CLINVAR = T)$oncoMatrix
  genes1=intersect(genes,rownames(oncomat))
  if(add_missing){
    missg=genes[!(genes %in% genes1)]
    missgmat=matrix(ncol = ncol(oncomat),nrow = length(missg),data = "")
    rownames(missgmat)=missg
    missgmat=rbind(oncomat,missgmat)
    oncomat=missgmat[match(genes,rownames(missgmat)), , drop=F] 
  }else{
  oncomat <- oncomat[match(genes1,rownames(oncomat)), , drop=F]
  }
  onco_genes <- rownames(oncomat)
  
  oncomat=gsub("_"," ",oncomat)
  oncomat=renameoncomat(oncomat)
  
}

names(mutation_colors) <- gsub("_"," ",names(mutation_colors))
names(mutation_colors)=renameoncomat(names(mutation_colors))
names(alter_fun)=renameoncomat(names(alter_fun))
