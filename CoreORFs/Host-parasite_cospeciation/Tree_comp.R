# install.packages('dendextend')
# install.packages("colorspace")
# install.packages("ape")
# install.packages("Rcpp")
# install.packages("phytools")
# install.packages("RRphylo")
# install.packages("manipulate")
# install.packages("doParallel")
# BiocManager::install("DECIPHER")
# install.packages("optparse")
# install.packages("ade4")
# install.packages("adephylo")


library(dendextend)
#library(colorspace)
library(ape)
library(Rphylip)
library(Rcpp)
library(phytools)
library(RRphylo)
#library(manipulate)
#library(DECIPHER)
#library(gtools)
library(doParallel)
library("optparse")
library("treespace")
library("ade4")
library("adephylo")

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-i", "--input_dir"), type="character", default=NA, 
              help="Path to the Input Trees Directory ", metavar="character"),
  make_option(c("-p", "--tree_file_pattern"), type="character", default='*.nwk$', 
              help="Pattern to find trees", metavar="character"),
  make_option(c("-o","--output_dir"), type="character", default=NA, 
              help="Path to the Output Directory", metavar="character"),
  make_option(c("-D","--lable2drop"), type="character", default=NA, 
              help="Tip to remove from all the trees", metavar="character"),
  make_option(c("-G","--host_outgroup"), type="character", default=NA, 
              help="Host OutGroup", metavar="character"),
  make_option(c("-B","--bac2host_file"), type="character", default=NA, 
              help="File connecting Host to Bacteria [tab delimited with 2 columns 'Bac' and 'Host' ]", metavar="character"),
  make_option(c("-g","--tree_outgroup"), type="character", default=NA, 
              help="Host OutGroup", metavar="character"),
  make_option(c("-c","--cpu"), type="numeric", default = 1, 
              help="Number of CPUs to Use ", metavar="character"),
  make_option(c("-n","--nboot"), type="numeric", default = 1000, 
              help="Number of Bootstraps to Use for P-Value Calculation ", metavar="character"),
  make_option(c("-r","--rooted_trees"), action="store_true",default=FALSE, 
              help="The Trees are Rooded!", metavar="character"),
  make_option(c("-m","--method"), type="character", default='correlation', 
              help="Method for calculating dist [default: correlation]", metavar="character"),
  make_option(c("-H", "--host_tree"), type="character", default=NA, 
              help="Path to the Host Tree File", metavar="character")
); 


#Get user information
opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", 
                                    option_list=option_list,
                                    epilogue="\n\nAuthor:Liron Levin");
opt = optparse::parse_args(opt_parser);


dir           = opt$input_dir
out_dir       = opt$output_dir
cpu           = opt$cpu
Nboots        = opt$nboot
root          = opt$rooted_trees
host_tree     = opt$host_tree
method        = opt$method
host_outgroup = opt$host_outgroup
tree_outgroup = opt$tree_outgroup
bac2host_file = opt$bac2host_file

original_tree_outgroup = tree_outgroup

bac2host = read.delim(bac2host_file)

reorder_tip = function(tree){
  is_tip <- tree$edge[,2] <= length(tree$tip.label)
  ordered_tips <- tree$tip.label[tree$edge[is_tip, 2]]
  return(ordered_tips)
}

untangle = function(tree1,tree2){
  tree1_roteted = tree1
  tree2_roteted = tree2
  original_dist = 0
  test_dist     = -1
  while (test_dist<original_dist){
    original_dist = sum(abs(sapply(1:length(tree1_roteted$tip.label),  function(x) x - which(reorder_tip(tree1_roteted)[x] == reorder_tip(tree2_roteted)))))
    temp_tree = tree1_roteted
    for (node in unique(temp_tree$edge[,1])){
      temp_temp_tree = ape::rotate(phy = temp_tree ,node)
      orig    = sum(abs(sapply(1:length(temp_tree$tip.label),  function(x) x - which(reorder_tip(temp_tree)[x] == reorder_tip(tree2)))))
      roteted = sum(abs(sapply(1:length(temp_temp_tree$tip.label),  function(x) x - which(reorder_tip(temp_temp_tree)[x] == reorder_tip(tree2)))))
      if (roteted<orig){
        temp_tree = temp_temp_tree
      }
    }
    tree1_roteted = temp_tree
    
    temp_tree = tree2_roteted
    for (node in unique(temp_tree$edge[,1])){
      temp_temp_tree = ape::rotate(phy = temp_tree ,node)
      orig    = sum(abs(sapply(1:length(temp_tree$tip.label),  function(x) x - which(reorder_tip(temp_tree)[x] == reorder_tip(tree1_roteted)))))
      roteted = sum(abs(sapply(1:length(temp_temp_tree$tip.label),  function(x) x - which(reorder_tip(temp_temp_tree)[x] == reorder_tip(tree1_roteted)))))
      if (roteted<orig){
        temp_tree = temp_temp_tree
      }
    }
    tree2_roteted = temp_tree
    
    test_dist = sum(abs(sapply(1:length(tree1_roteted$tip.label),  function(x) x - which(reorder_tip(tree1_roteted)[x] == reorder_tip(tree2_roteted)))))
  }
  # plotTreeDiff(ape::as.phylo(tree1_roteted), ape::as.phylo(tree2_roteted), use.edge.length=FALSE, 
  #             treesFacing = TRUE,sizeOfDifferences = T,tipMatch = F)
  return(c(tree1_roteted,tree2_roteted))
 
}

culculate_cor <- function(original_Host_tree,tree,bac2host,method,plot_dend=NA){ 
   Host_tree = original_Host_tree
   if ( method == 'dist.topo' ){
        Host_tree = original_Host_tree
        for (host in unique(bac2host[,"Host"])){
              temp =  reorder_tip(tree)[reorder_tip(tree)  %in%  bac2host[bac2host$Host==x,"Bac"]]
              dato = data.frame(bind=temp,reference=rep(host,length(temp)),poly=F)
              Host_tree = RRphylo::tree.merger(backbone=Host_tree,data=dato,plot=F)
              Host_tree = ape::drop.tip(phy = Host_tree,tip = host )
        }
        x <- untangle(Host_tree,tree)
    } else {
        
        for (host in unique(bac2host[,"Host"])){
            temp =  reorder_tip(tree)[reorder_tip(tree)  %in%  bac2host[bac2host$Host==x,"Bac"]]
            ref =host
            # permutations(length(temp), 1, temp, repeats = FALSE)
            # for (tip in sample(temp)){
            for (tip in temp){
              dato = data.frame(bind=tip,reference=ref,poly=F)
              Host_tree = RRphylo::tree.merger(backbone=Host_tree,data=dato,plot=F)
              ref = tip
            }
            Host_tree = ape::drop.tip(phy = Host_tree,tip = host )
        }
        
        #Host_tree$edge.length = Host_tree$edge.length+0.0000001
        
        Host_tree_r = stats::as.dendrogram(Host_tree)
        tree_r      = stats::as.dendrogram(tree)
        
        dl <- dendextend::dendlist(Host_tree_r, tree_r)
        x = dendextend::untangle(dl,method = "step2side")
        
        new_tree = ape::as.phylo(x[[1]])
        
        Host_tree = original_Host_tree
        tip.label = rev(reorder_tip(new_tree))
        for (host in unique(bac2host[,"Host"])){
            temp =  reorder_tip(tree)[reorder_tip(tree)  %in%  bac2host[bac2host$Host==x,"Bac"]]
            ref =host
            for (tip in temp){
              dato = data.frame(bind=tip,reference=ref,poly=F)
              Host_tree = RRphylo::tree.merger(backbone=Host_tree,data=dato,plot=F)
              ref = tip
            }
            Host_tree = ape::drop.tip(phy = Host_tree,tip = host )
        }
        #Host_tree$edge.length = Host_tree$edge.length+0.0000001
        Host_tree_r = stats::as.dendrogram(Host_tree)
        dl <- dendextend::dendlist(Host_tree_r, tree_r)
        x = dendextend::untangle(dl,method = "step2side")
    }

  if (!is.na(plot_dend)){
    # pdf(file = plot_dend)
    # tanglegram(x,center = F,
    # common_subtrees_color_lines = FALSE, highlight_distinct_edges  = FALSE, highlight_branches_lwd = FALSE,
    # margin_inner = 10,rank_branches = T,lab.cex=0.5,columns_width=c(5,1,5))
    # dev.off()
    if (method!='dist.topo') {
      pdf(file = plot_dend)
      tanglegram(x,center = F,
                 common_subtrees_color_lines = FALSE, highlight_distinct_edges  = FALSE, highlight_branches_lwd = FALSE,
                 margin_inner = 10,rank_branches = T,lab.cex=0.5,columns_width=c(5,1,5))
      dev.off()
      ape::write.tree(phy =ape::as.phylo(x[[1]]) ,file = paste(plot_dend,'.HOST.nwk') )
      ape::write.tree(phy =ape::as.phylo(x[[2]]) ,file = paste(plot_dend,'.nwk') )
    } else {
      
      pdf(file = plot_dend,width =20 ,height =20 )
      plotTreeDiff(x[[1]], x[[2]], use.edge.length=FALSE, 
                   treesFacing = TRUE,sizeOfDifferences = T,tipMatch = F)
      dev.off()
      ape::write.tree(phy =ape::as.phylo(x[[1]]) ,file = paste(plot_dend,'.HOST.nwk') )
      ape::write.tree(phy =ape::as.phylo(x[[2]]) ,file = paste(plot_dend,'.nwk') )
    }
  }
  if (method=='correlation'){
    tree_dist = dendextend::cor.dendlist(dl)[2,1]
  }else if (method == 'treeDist' ){
    tree_dist = treeDist(ape::as.phylo(x[[1]]), ape::as.phylo(x[[2]]))
  }else if (method=='gammaindex'){
    tree_dist =  cor_bakers_gamma(ape::as.phylo(x[[1]]), ape::as.phylo(x[[2]]))
  }else if (method=='dist.topo'){
    tree_dist = dist.topo(x[[1]], x[[2]])[1]
  }
  
  return(tree_dist)
}


data = data.frame()

cl <- makeCluster(cpu)
registerDoParallel(cl)
getDoParWorkers()

#tree      = ReadDendrogram(file = 'Trees/Rflav.proteinortho.tsv.OrthoGroup2382.fasta.aln.contree')
original_Host_tree = ape::read.tree(file = host_tree)

for (file in list.files(path =dir ,pattern = opt$tree_file_pattern)){
  print(file.path(dir,file))
  tree          = ape::read.tree(file = file.path(dir,file))
  if (method == "mantel"){
    if (!is.na(opt$lable2drop)){
      tree = ape::drop.tip(phy = tree,tip = opt$lable2drop )
    }
    Host_tree = original_Host_tree
    for (host in unique(bac2host[,"Host"])){
      # temp = reorder_tip(tree)[stringr::str_starts(string = reorder_tip(tree),pattern = host) ]
      temp =  reorder_tip(tree)[reorder_tip(tree)  %in%  bac2host[bac2host$Host==x,"Bac"]]
      # dato = data.frame(bind=temp,reference=rep(Host_lables[host],length(temp)),poly=F)
      dato = data.frame(bind=temp,reference=rep(host,length(temp)),poly=F)
      Host_tree = RRphylo::tree.merger(backbone=Host_tree,data=dato,plot=F)
      Host_tree = ape::drop.tip(phy = Host_tree,tip = host )
    }
    host_dist = adephylo::distTips(x = Host_tree,method = "patristic")
    tree_dist = adephylo::distTips(x = tree,method = "patristic")
    
    tree_dist = as.matrix(tree_dist)
    host_dist = as.matrix(host_dist)
    host_dist <- host_dist[rownames(tree_dist),colnames(tree_dist)]
    results = ade4::mantel.rtest(m1 = as.dist(host_dist),as.dist(tree_dist),nrepet = Nboots)
    data[file,'Original'] = results$obs
    data[file,'P_Value']  = results$pvalue
    write.csv(tree_dist,file =file.path(out_dir,paste(file,'.dist.csv',sep='') ))
    write.csv(host_dist,file =file.path(out_dir,paste("Host",'.dist.csv',sep='') ))
    
  }else{
      tree$edge.length = tree$edge.length+0.0000001
        if (!root){
            tree_outgroup = tree$tip.label[stringr::str_starts(string = tree$tip.label,pattern = original_tree_outgroup) ]
            tree          = chronos(root(tree              ,outgroup = tree_outgroup,resolve.root = T),lambda=0)
            if (!ape::is.rooted(original_Host_tree)){
                Host_tree     = chronos(root(original_Host_tree,outgroup = host_outgroup,resolve.root = T),lambda=0)
            }else{
                Host_tree     = chronos(original_Host_tree)
            }
        }else{
            tree          = chronos(tree)
            Host_tree     = chronos(original_Host_tree)
        }
        
        if (!is.na(opt$lable2drop)){
         tree = ape::drop.tip(phy = tree,tip = opt$lable2drop )
        }
      original_cor = culculate_cor(Host_tree,tree,Host_lables,method,plot_dend = paste(file.path(out_dir,file),'.pdf') )
      
      pack = gtools::loadedPackages(silent = T)
      cor_dist = foreach(x=1:Nboots , .combine=rbind,.packages=pack$Name,.verbose=F) %dopar% {
        temp_tree = tree
        if (!root){
            temp_tree$tip.label = c(temp_tree$tip.label[1] , sample(temp_tree$tip.label[-1]))
        }else{
            temp_tree$tip.label = sample(temp_tree$tip.label)
        }
        cor_dist = culculate_cor(Host_tree,temp_tree,Host_lables,method)
        print(cor_dist)
      }
      #hist(cor_dist)
      data[file,'Original'] = original_cor
      if ((method == "correlation") || (method == "gammaindex") ){
        data[file,'P_Value']       = sum(cor_dist > original_cor)/length(cor_dist)
      }else{
        data[file,'P_Value']       = sum(cor_dist < original_cor)/length(cor_dist)
      }
      data[file,'Rendom_Dist']   = list(list(as.vector(cor_dist)))
      
      
    }
}
save(data,file =file.path(out_dir,paste('Results_',method,sep='') ))
if (method != "mantel"){
    data[,3] = paste(data[,3],collapse = ',',sep = ',')
}
write.csv(data,file =file.path(out_dir,paste('Results_',method,'.csv',sep='') ))
stopCluster(cl)




# entanglement(x)
# 
# 
# 
# 
# 
# tanglegram(x)
# tanglegram(dl)
# cor.dendlist(dl)
# cor_bakers_gamma(Host_tree_r, tree_r)
# 


# host_dist = adephylo::distTips(x = Host_tree,method = "sumDD")
# 
# host_dist =as.matrix(adephylo::distTips(x = x[[1]],method = "sumDD"))
# tree_dist =as.matrix(adephylo::distTips(x = tree,method = "sumDD"))
# 
# ape::mantel.test(host_dist,tree_dist,graph = TRUE)


