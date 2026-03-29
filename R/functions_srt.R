#' Testing each row against the others
#' @param x Matrix with a function per row
#' @param method Testing procedure
#' @param B Number of resamples for bootstrap
#' 
oneRow = function(x,method=c("all","bootstrap"),B=NULL){
    
    if(method == "all"){
        sig = NULL
        for(k in 1:nrow(x)){
            cat("Row ",k,"\n")
            env = apply(x[-k,],2,range,na.rm=TRUE)
            if(any(x[k,] < env[1,] | x[k,] > env[2,],na.rm=TRUE))
                sig = c(sig,k)
        }
    }

    if(method == "bootstrap"){
        sig = NULL
        rows = 1:nrow(x)
        for(k in 1:nrow(x)){
            cat("Row ",k,"\n")
            sel = sample(rows[-k],B,replace=TRUE)
            env = apply(x[sel,],2,range,na.rm=TRUE)
            if(any(x[k,] < env[1,] | x[k,] > env[2,],na.rm=TRUE))
                sig = c(sig,k)
        }
    }
    sig
}

#' Testing the direction of the test
#' @param x Matrix with a function per row
#' @param sig Significant associations previously calculated
#'  
oneRowSide = function(x,sig){
    side = matrix(0,nrow=length(sig),ncol=ncol(x))
        for(k in 1:length(sig)){
            env = apply(x[-sig[k],],2,range,na.rm=TRUE)
            side[k,which(x[sig[k],] < env[1,])] = -1
            side[k,which(x[sig[k],] > env[2,])] = 1
        }
    side
}




#' A multitype point point is made from datasets
#' @param dirbase Directory for the files
#' @param case Case to be analyzed
#' @param ff Fraction of the points to be used
#' @param type Celular type
#' @param min.n Minimum number of points
#' @param num.genes Number of genes to be analyzed (100 or 423)
#' @export
do.ppp = function(dirbase,case,ff=0.02,type=NULL,min.n = NULL,num.genes=100){
    namebase = paste0("Case",case,"_")
    cases = paste0(c(paste0("00",0:9),paste0("0",10:99)),".csv")
    if(num.genes >100)
        cases = paste0(c(paste0("00",0:9),paste0("0",10:99),100:(num.genes-1)),".csv")
    
    casesn = paste0(dirbase,namebase,cases)
    
    x = vector("list",num.genes)
    for(i in 1:num.genes){
        names(x)[i] = paste0("x",i)
        x[[i]] = NA
        if(file.exists(casesn[i])){
            u = read.csv(casesn[i],header=TRUE,sep=",")
            if(!is.null(type)) u = u[u[,5] == type,]
            if(ff < 1){
                n = round(nrow(u)*ff)
                sel = sample(nrow(u),n)
            }
            if(ff > .99){
                n = nrow(u)
                sel = 1:nrow(u)
            }
            xlimits = range(u[sel,1])
            ylimits = range(u[sel,2])
            if(nrow(u) > min.n)
                x[[i]] =  ppp(u[sel,1], u[sel,2],xlimits,ylimits)
            else
                x[[i]] = NA
        }
    }
    x
}

#' Graph showing important associations
#' @description
#' Graph showing important associations
#' @param file_x File name with the multitype point pattern
#' @param file_sig File names with the vector with the important pairs
#' @param file_sig_side File names with the matrix where each row shows the
#' @param file_annotation File name with the annotation file
#'  
do_graph = function(file_x,file_sig,file_sig_side,file_annotation){
    load(file_annotation)
    load(file_x)
    sel = which(!is.na(x))
    pairs1 = t(combn(sel,2))
    load(file_sig)
    pairs1 = pairs1[sig,]
    sel1 = sort(unique(c(pairs1[,1],pairs1[,2])))
    n = matrix(0,nrow=length(sel1),ncol=length(sel1))
    rownames(n) = sel1
    colnames(n) = sel1
    for(k in 1:length(sig))
        n[rownames(n) == pairs1[k,1],colnames(n) == pairs1[k,2]] = 1

    load(file_sig_side)
    type = NULL
    for(i in c(-1,1))
        type = cbind(type,apply(side,1,function(x) sum(x == i)))

    out = rep(NA,nrow(type))
    for(k in 1:nrow(type)){
        if(type[k,1] >0 & type[k,2] >0) out[k] = 3 # up and low
        if(type[k,1] >0 & type[k,2] ==0)  out[k] = 1 ## low
        if(type[k,1] == 0 & type[k,2] > 0)  out[k] = 2 ## up
    }
    out_length = 10*apply(type,1,sum)/ncol(side)
     
    g = graph_from_adjacency_matrix(
        adjmatrix = n,
        mode = "max",
        weighted = NULL,
        diag = TRUE,
        )

    pairs1_hgnc_symbol = NULL
    for(k in 1:nrow(pairs1))
        pairs1_hgnc_symbol = rbind(pairs1_hgnc_symbol,
            c(hgnc_symbol[pairs1[k,1]],hgnc_symbol[pairs1[k,2]]))

    
    V(g)$label =  hgnc_symbol[sel1]
    V(g)$label.cex = 1
    E(g)$label1 = pairs1_hgnc_symbol[,1]
    E(g)$label2 = pairs1_hgnc_symbol[,2]
    E(g)$width  =  out_length
    E(g)$out = out
    E(g)$color[E(g)$out == 1] = "red"
    E(g)$color[E(g)$out == 2] = "black"
    E(g)$color[E(g)$out == 3] = "blue"

    
   
    g
}


#' Envelopes
#' @description
#' Given three matrices  and a pair of genes the envelope of the rest of pairs
#' is evaluated
#' @param i Case
#' @param j Window
#' @param ff Sampling fraction
#' @param desc Descriptive used 
#' @param genes Pair of genes 
#' @param dir_figures Directory where to save the figures
#' @param file_annotation Annotation file
#' @export 

do_env= function(xd,i,j,ff,desc,genes,dir_figures,file_annotation){
    file_x = paste0("../data/local/x_l_case_",i,"_area_",j,"_ff_",1,".rda")
    load(file_x)
    sel = which(!is.na(x)) 
    pairs1 = t(combn(sel,2))
    x1 = x[sel]
    pairs11 =  t(combn(length(sel),2))
    file_sig = paste0("results/local/sig_",desc,"_l_",i,"_ff_",ff,"_area_",j,".rda")
    load(file_sig)
    file_desc = paste0("results/local/",desc,"_l_",i,"_ff_",ff,"_area_",j,".rda")
    load(file_desc)
    switch(desc,
           "kc" = {
               xd = kc
           },
           "Lc" = {
               xd = Lc
           },
          "pcf" = {
               xd = pcf
           })
    load(file_annotation)
    pair0 = range(which(hgnc_symbol == genes[1]),which(hgnc_symbol == genes[2]))
    row0 = which(pairs1[,1] == pair0[1] & pairs1[,2] == pair0[2] )
    df = data.frame(x0 = 1:ncol(xd),y0=xd[row0,],low0 = apply(xd[-row0,],2,min),
                    up0 = apply(xd[-row0,],2,max))
    df = na.omit(df)
        p = ggplot(df,aes(x=x0,y=y0)) + geom_line() +
            geom_ribbon(data=df, aes(x=x0, ymin=low0, ymax=up0),alpha=.2) + 
            xlab("distance") + ylab("Descriptive") 
 
    file_env = paste0(dir_figures,"Case_",i,"_Area_",j,"_",genes[1],"_",genes[2],
                      "_",desc,"_ff_",ff,"_env.png")
    ggsave(file_env,p)
    p
}
 
