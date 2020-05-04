###################################################
###
work.directory <- "D:/network_alignment/"
DATA.FOLD <- paste0(work.directory,'Data/')
NETWORK.FOLD <- paste0(DATA.FOLD,'R_output/')
if(!file.exists(NETWORK.FOLD))
{
dir.create(NETWORK.FOLD)
}
NETWORK.FOLD1 <- paste0(DATA.FOLD,'R_input/')
if(!file.exists(NETWORK.FOLD1))
{
dir.create(NETWORK.FOLD1)
}
NETWORK.FOLD2 <- paste0(DATA.FOLD,'JAR_output/')
if(!file.exists(NETWORK.FOLD2))
{
dir.create(NETWORK.FOLD2)
}
NETWORK.FOLD3 <- paste0(work.directory,'FINAL_result/')
if(!file.exists(NETWORK.FOLD3))
{
dir.create(NETWORK.FOLD3)
}
GOKEGG.FOLD <- paste0(DATA.FOLD, "GOKEGG/")
if(!file.exists(GOKEGG.FOLD)){dir.create(GOKEGG.FOLD)}

###########################################################################################################
###########################################################################################################
#################################Network Division#########################################################
###
if (!requireNamespace("BiocManager", quietly = TRUE))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/CRAN/") 
    install.packages("BiocManager")
### 
BiocManager::install("org.Hs.eg.db")
BiocManager::install("GO.db")
###
library(org.Hs.eg.db)
library(GO.db)
###
GOTERM <- as.list(GOTERM)
###
go2allegs <- as.list(org.Hs.egGO2ALLEGS)
length(go2allegs) 
eg2symbol <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(eg2symbol)
eg2symbol <- as.list(eg2symbol[mapped_genes])

BP.file<- paste0(GOKEGG.FOLD,'hsa.GO.BP.txt')

for(i in 1:length(go2allegs)){
go.id<- names(go2allegs[i])
if(length(GOTERM[[ go.id ]])>0){
go.name<- Term(GOTERM[[ go.id ]])
go.type<- Ontology(GOTERM[[ go.id ]])
go.eg<- as.character(go2allegs[[i]])
go.symbol<- unique(as.character(eg2symbol[go.eg]))

if(go.type=='BP'){ 
sink(BP.file,append = T)
cat(paste(go.id,go.name,sep='_'),'\t',go.type,'\t',sep = '')
for(j in 1:length(go.symbol)){
cat(go.symbol[j],'\t',sep='')
}
cat('\n')
sink()
}
}
}
					
###################################################
###
bpjg=list()
m=0
bp<- file(paste0(GOKEGG.FOLD,'hsa.GO.BP.txt'), "r")
line=readLines(bp,n=1)
while(length(line) != 0 ) {
   m=m+1
   bpjg[[m]]=line  
   line=readLines(bp,n=1)

}
close(bp)
###
fqu=list()
symbollist=list() 
for(i in 1:length(bpjg)){
fqu[[i]]=unlist(strsplit(bpjg[[i]], "\t"))
symbollist[[i]]=fqu[[i]][3:length(fqu[[i]])]
}

###################################################         
###
HPRD.fname <- paste0(DATA.FOLD, "PPI network1.txt")
YEAST.fname <- paste0(DATA.FOLD, "PPI network2.txt")
rtab<- read.delim(HPRD.fname,header=F,stringsAsFactors=F)
ytab<- read.delim(YEAST.fname,header=F,stringsAsFactors=F)
###
rtab<- rtab[rtab[,1]!=rtab[,2],]
ytab<- ytab[ytab[,1]!=ytab[,2],]
###
rjs=vector()
yjs=vector()  
for(i in 1:length(symbollist)){		
interaction<- rtab[(rtab[,1]%in% symbollist[[i]])&(rtab[,2]%in%symbollist[[i]]),]

if(nrow(interaction)>0){
  rjs=rbind(rjs,i)
  HPRD.interaction.pair <- apply(interaction, 1, function(current.row)
	{
	return(paste0(current.row[1], "_", current.row[2]))
	})
				
	HPRD.interaction.pair.unique <- unique(HPRD.interaction.pair)						
	HPRD.unique.ind <- match(HPRD.interaction.pair.unique, HPRD.interaction.pair)				
	HPRD.net <- interaction[HPRD.unique.ind,]
						
	write.table(HPRD.net,file=paste0(NETWORK.FOLD,paste0('HPRD',i,'.txt')), row.names=F, quote=F, sep='\t')
}
interaction<- ytab[(ytab[,1]%in% symbollist[[i]])&(ytab[,2]%in%symbollist[[i]]),]		

if(nrow(interaction)>0){
 yjs=rbind(yjs,i)
 YEAST.interaction.pair <- apply(interaction, 1, function(current.row)
 {
	 return(paste0(current.row[1], "_", current.row[2]))
 })
				
	YEAST.interaction.pair.unique <- unique(YEAST.interaction.pair)						
	YEAST.unique.ind <- match(YEAST.interaction.pair.unique,YEAST.interaction.pair)				
	YEAST.net <- interaction[YEAST.unique.ind,]
						
	write.table(YEAST.net, file=paste0(NETWORK.FOLD,paste0('YEAST',i,'.txt')), row.names=F, quote=F, sep='\t')
 }

}

###################################################			
##
JJI=intersect(rjs,yjs)

###################################################		
###
gs=rep(0,length(JJI))
for(i in 1:length(JJI)){
jbp<-read.table(paste0(NETWORK.FOLD,paste0('HPRD',JJI[i],'.txt')),header=T)
gs[i]=length(unique(unlist(jbp)))
}
###
mwz=which(gs==max(gs))

###
jbp<-read.table(paste0(NETWORK.FOLD,paste0('HPRD',JJI[mwz],'.txt')),header=T)
write.table(jbp,file= paste0(NETWORK.FOLD1,'hsa.GO.RBH_max.txt'), row.names=F, quote=F, sep='\t')					
jbp<-read.table(paste0(NETWORK.FOLD,paste0('YEAST',JJI[mwz],'.txt')),header=T)    
write.table(jbp,file= paste0(NETWORK.FOLD1,'hsa.GO.YBH_max.txt'), row.names=F, quote=F, sep='\t')					

##################################################################################################
###################################################	
###################################################		
###
rtab<- read.delim(HPRD.fname,header=F,stringsAsFactors=F)		
ytab<- read.delim(YEAST.fname,header=F,stringsAsFactors=F)	
###
rtab<- rtab[rtab[,1]!=rtab[,2],]
ytab<- ytab[ytab[,1]!=ytab[,2],]
###
urgene=unique(unlist(rtab))
uygene=unique(unlist(ytab))
###
install.packages("foreach") 
install.packages("doParallel") 
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl)
###
fun <- function(x)
{
fun=which(x==urgene)
return(fun)
}
###
rg=unlist(rtab)
hrg<- foreach(x=1:length(rg),.combine='rbind') %dopar% fun(rg[x])
stopCluster(cl)
###
bhrg=matrix(hrg,nrow(rtab),2)
for(i in 1:nrow(bhrg))
{
if(as.numeric(bhrg[i,1])>as.numeric(bhrg[i,2]))
{
k=bhrg[i,1]
bhrg[i,1]=bhrg[i,2]
bhrg[i,2]=k
}
}
ubhrg=unique(bhrg)
nrow(ubhrg)
rhzd <- apply(bhrg, 1, function(current.row)
{
return(paste0(current.row[1],'_',current.row[2]))
})
urhzd=unique(rhzd)
rind <- match(urhzd,rhzd)
qcrg <- rtab[rind,]
nrow(qcrg)
rphfile<-paste0(NETWORK.FOLD1,'wholehumanquchonginteraction.txt')
sink(rphfile,append = T)
for(i in 1:nrow(qcrg))
{
cat(qcrg[i,1])
cat('\t')
cat(qcrg[i,2]) 
cat('\n') 
}  
sink() 

###################################################	
###
cl <- makeCluster(8)
registerDoParallel(cl)
###
fun <- function(x){
fun=which(x==uygene)
return(fun)
}
yg=unlist(ytab)
hyg<- foreach(x=1:length(yg),.combine='rbind') %dopar% fun(yg[x])
stopCluster(cl)
bhyg=matrix(hyg,nrow(ytab),2)
for(i in 1:nrow(bhyg))
{
if(as.numeric(bhyg[i,1])>as.numeric(bhyg[i,2]))
{
k=bhyg[i,1]
bhyg[i,1]=bhyg[i,2]
bhyg[i,2]=k
}
}
ubhyg=unique(bhyg)
nrow(ubhyg)
yhzd <- apply(bhyg, 1, function(current.row)
{
return(paste0(current.row[1],'_',current.row[2]))
})
uyhzd=unique(yhzd)
yind <- match(uyhzd,yhzd)
qcyg <- ytab[yind,]
nrow(qcyg)
yphfile<-paste0(NETWORK.FOLD1,'wholeyeastquchonginteraction.txt')
sink(yphfile,append = T)
for(i in 1:nrow(qcyg))
{
cat(qcyg[i,1])
cat('\t')
cat(qcyg[i,2]) 
cat('\n') 
}  
sink() 


###################################################	
###################################################		
###    
RHZG<- read.delim(paste0(NETWORK.FOLD1,'hsa.GO.RBH_max.txt'),header=T,stringsAsFactors=F)
urgene=unique(unlist(RHZG))
####################################
###
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl)
###
fun <- function(x)
{
fun=which(x==urgene)
return(fun)
}
###
rg=unlist(RHZG)
hrg<- foreach(x=1:length(rg),.combine='rbind') %dopar% fun(rg[x])
stopCluster(cl)
###
bhrg=matrix(hrg,nrow(RHZG),2)
for(i in 1:nrow(bhrg))
{
if(as.numeric(bhrg[i,1])>as.numeric(bhrg[i,2]))
{
k=bhrg[i,1]
bhrg[i,1]=bhrg[i,2]
bhrg[i,2]=k
}
}
ubhrg=unique(bhrg)
###
rphfile<-paste0(NETWORK.FOLD1,'humanquchonginteraction_max.txt')
sink(rphfile,append = T)
for(i in 1:nrow(ubhrg))
{
cat(ubhrg[i,1],'\t')
cat(ubhrg[i,2],'\n') 
}  
sink() 

###################################################	
YHZG<- read.delim(paste0(NETWORK.FOLD1,'hsa.GO.YBH_max.txt'),header=T,stringsAsFactors=F)
uygene=unique(unlist(YHZG))
####################################
cl <- makeCluster(8)
registerDoParallel(cl)
###
fun <- function(x){
fun=which(x==uygene)
return(fun)
}
###
yg=unlist(YHZG)
hyg<- foreach(x=1:length(yg),.combine='rbind') %dopar% fun(yg[x])
stopCluster(cl)
###
bhyg=matrix(hyg,nrow(YHZG),2)
for(i in 1:nrow(bhyg))
{
if(as.numeric(bhyg[i,1])>as.numeric(bhyg[i,2]))
{
k=bhyg[i,1]
bhyg[i,1]=bhyg[i,2]
bhyg[i,2]=k
}
}
ubhyg=unique(bhyg)
###
yphfile<-paste0(NETWORK.FOLD1,'yeastquchonginteraction_max.txt')
sink(yphfile,append = T)
for(i in 1:nrow(ubhyg))
{
cat(ubhyg[i,1],'\t')
cat(ubhyg[i,2],'\n') 
}  
sink() 




###########################################################################################################
###########################################################################################################
#################################Network Alignment#########################################################
###
tygene<- read.delim(paste0(DATA.FOLD,'mart_export.txt'),header=F,stringsAsFactors=F)
###
RHZG<- read.delim(paste0(NETWORK.FOLD1,'hsa.GO.RBH_max.txt'),header=T,stringsAsFactors=F)
urgene=unique(unlist(RHZG))
YHZG<- read.delim(paste0(NETWORK.FOLD1,'hsa.GO.YBH_max.txt'),header=T,stringsAsFactors=F)
uygene=unique(unlist(YHZG))
###
rphfile<-paste0(NETWORK.FOLD1,'humanquchonginteraction_max.txt')
ubhrg<-read.delim(rphfile,header=F,stringsAsFactors=F)
yphfile<-paste0(NETWORK.FOLD1,'yeastquchonginteraction_max.txt')
ubhyg<-read.delim(yphfile,header=F,stringsAsFactors=F)
###
tyjz=matrix(0,nrow=length(urgene),ncol=length(uygene))
for(i in 1:length(urgene))
{ 
tyjz[i,which(urgene[i]==uygene)]=1
}      
rwz=list()
for(i in 1:length(urgene))
{
rwz[[i]]=which(tygene[,1]==urgene[[i]])
}
ywz=list()
for(j in 1:length(uygene))
{
ywz[[j]]=which(tygene[,2]==uygene[[j]])
}     
for(i in 1:length(urgene))
{
for(j in 1:length(uygene)) 
{  
wz=intersect(rwz[[i]],ywz[[j]])
if(length(wz)==1)
{
tyjz[i,j]=1
} 
}
}
rtab<- read.delim(rphfile,header=F,stringsAsFactors=F)
rljjz=matrix(0,nrow=length(urgene),ncol=length(urgene))
for(i in 1:nrow(rtab))
{
rljjz[rtab[i,1],rtab[i,2]]=1
rljjz[rtab[i,2],rtab[i,1]]=1  
}
ytab<- read.delim(yphfile,header=F,stringsAsFactors=F)
yljjz=matrix(0,nrow=length(uygene),ncol=length(uygene))
for(i in 1:nrow(ytab))
{
yljjz[ytab[i,1],ytab[i,2]]=1
yljjz[ytab[i,2],ytab[i,1]]=1  
}
rljgene=list()
for(i in 1:length(urgene))
{
rljgene[[i]]=which(rljjz[,i]==1)
}
yljgene=list()
for(i in 1:length(uygene))
{
yljgene[[i]]=which(yljjz[,i]==1)
}
xtyjz=matrix(0,nrow=length(urgene),ncol=length(uygene))
for(i in 1:length(urgene))
{
for(j in 1:length(uygene))
{  
tyh=sum(tyjz[rljgene[[i]],yljgene[[j]]])    
xtyjz[i,j]=tyjz[i,j]+tyh/(length(rljgene[[i]])*length(yljgene[[j]]))
}
}
human_command <- paste0('java -jar ',work.directory,'undigraph.jar ',NETWORK.FOLD1,'humanquchonginteraction_max.txt -all ',NETWORK.FOLD2,'humanresult_max.txt')
yeast_command <- paste0('java -jar ',work.directory,'undigraph.jar ',NETWORK.FOLD1,'yeastquchonginteraction_max.txt -all ',NETWORK.FOLD2,'yeastresult_max.txt')
system(human_command)
system(yeast_command)
ROBIT<-read.table(paste0(NETWORK.FOLD2,'humanresult_max.txt'),header=F)
rojz=matrix(nrow=length(urgene),ncol=73)
rojz=ROBIT
YOBIT<-read.table(paste0(NETWORK.FOLD2,'yeastresult_max.txt'),header=F)
yojz=matrix(nrow=length(uygene),ncol=73)
yojz=YOBIT
###
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl)
###
fun <- function(x){
if(x%%length(uygene)>0){
return(1/(1+sqrt(sum((rojz[ceiling(x/length(uygene)),]-yojz[x%%length(uygene),])^2))))
}
if(x%%length(uygene)==0){
return(1/(1+sqrt(sum((rojz[ceiling(x/length(uygene)),]-yojz[length(uygene),])^2))))
}
}
###
ryoxl<- foreach(x=1:(length(urgene)*length(uygene)),.combine='rbind') %dopar% fun(x)
stopCluster(cl)
###
ryojz=matrix(ryoxl,length(urgene),length(uygene),T);
rdn=sum(tyjz);
xs1=rdn/min(length(urgene),length(uygene));
xs2=1-xs1;
zjz=xs1*xtyjz+xs2*ryojz;
write.table(zjz,file= paste0(NETWORK.FOLD3,'lianpeijuzhen_max.txt'), row.names=F, col.names=F, quote=F, sep='\t')
write.table(tyjz,file= paste0(NETWORK.FOLD3,'tongyuanjuzhen_max.txt'), row.names=F, col.names=F, quote=F, sep='\t')

###################################
matlablj='D:/network_alignment/cs.exe'
system(matlablj)

###################################
lpjgjz<-read.table(paste0(NETWORK.FOLD3,'alignment_result_max.txt'),header=F,sep=',')
mpgene=list()
for(i in 1:length(uygene)){
mpwz=which(lpjgjz[,i]==1)
mpgene[[i]]=cbind(urgene[mpwz],uygene[i])
}
###
sink(paste0(NETWORK.FOLD,'new match-pair gene.txt'),append = F)
for(i in 1:length(uygene)){
    cat(mpgene[[i]][1])
    cat('\t') 
    cat(mpgene[[i]][2])
    cat('\n') 
}  
sink()
###
sink(paste0(NETWORK.FOLD,'all match-pair gene.txt'),append = T)
for(i in 1:length(uygene)){
    cat(mpgene[[i]][1])
    cat('\t') 
    cat(mpgene[[i]][2])
    cat('\n') 
}  
sink()
###################################
#####
wqcrg<- read.delim(paste0(NETWORK.FOLD1, "wholehumanquchonginteraction.txt"),header=F,stringsAsFactors=F)		
wqcyg<- read.delim(paste0(NETWORK.FOLD1, "wholeyeastquchonginteraction.txt"),header=F,stringsAsFactors=F)
###
rmtg='a'
ymtg='b'
###
while(length(rmtg)>0&length(ymtg)>0){
###
mpgene<-read.delim(paste0(NETWORK.FOLD,'new match-pair gene.txt'),header=F,stringsAsFactors=F)
ampgene<-read.delim(paste0(NETWORK.FOLD,'all match-pair gene.txt'),header=F,stringsAsFactors=F)
###
###
rmpg=list()
ympg=list()
for(i in 1:nrow(mpgene)){
rmpg[[i]]=mpgene[i,1]
ympg[[i]]=mpgene[i,2]
}
rmpg=unique(unlist(rmpg))
ympg=unique(unlist(ympg))
###all
armpg=list()
aympg=list()
for(i in 1:nrow(ampgene)){
armpg[[i]]=ampgene[i,1]
aympg[[i]]=ampgene[i,2]
}
armpg=unique(unlist(armpg))
aympg=unique(unlist(aympg))
###
rsg=list()
ysg=list()
for(p in 1:nrow(mpgene)){
rljwz1=which(wqcrg[,1]==rmpg[p])
rljwz2=which(wqcrg[,2]==rmpg[p])
rljwz=c(rljwz1,rljwz2)
###
rljhz=wqcrg[rljwz,]
zrgene=unique(unlist(rljhz))
rsg[[p]]=zrgene
###############################
yljwz1=which(wqcyg[,1]==ympg[p])
yljwz2=which(wqcyg[,2]==ympg[p])
yljwz=c(yljwz1,yljwz2)
###
yljhz=wqcyg[yljwz,]
zygene=unique(unlist(yljhz))
ysg[[p]]=zygene
}
rtg=unique(unlist(rsg))
ytg=unique(unlist(ysg))
##################################
###
rmtg=setdiff(rtg,armpg)
ymtg=setdiff(ytg,aympg)
length(rmtg)
length(ymtg)

##################################
###
if(length(rmtg)>0&length(ymtg)>0){
###
ntyjz=matrix(0,nrow=length(rmtg),ncol=length(ymtg))
for(i in 1:length(rmtg))
{ 
ntyjz[i,which(rmtg[i]==ymtg)]=1
}      
rwz=list()
for(i in 1:length(rmtg))
{
rwz[[i]]=which(tygene[,1]==rmtg[[i]])
}
ywz=list()
for(j in 1:length(ymtg))
{
ywz[[j]]=which(tygene[,2]==ymtg[[j]])
}     
for(i in 1:length(rmtg))
{
for(j in 1:length(ymtg)) 
{  
wz=intersect(rwz[[i]],ywz[[j]])
if(length(wz)==1)
{
ntyjz[i,j]=1
} 
}
}
###
write.table(ntyjz,file= paste0(NETWORK.FOLD3,'neighbor-tongyuanjuzhen(0-1).txt'), row.names=F, col.names=F, quote=F, sep='\t')

###
nppjz=matrix(0,nrow=length(rmtg),ncol=length(ymtg))
###
for(p in 1:nrow(mpgene)){
###
rnwz=vector()
for(i in 1:length(rsg[[p]])){
rop=which(rmtg==rsg[[p]][i])
if(length(rop)>0){
rnwz=cbind(rnwz,rop)
}
}
###
ynwz=vector()
for(j in 1:length(ysg[[p]])){
yop=which(ymtg==ysg[[p]][j])
if(length(yop)>0){
ynwz=cbind(ynwz,yop)
}
}
###
nppjz[rnwz,ynwz]=1
}

###
###
for(i in 1:nrow(ntyjz)){
dh=ntyjz[i,]
dh[dh==1]=4.4
dh[dh==0]=-1.6
ntyjz[i,]=dh
}
###
write.table(ntyjz,file= paste0(NETWORK.FOLD3,'neighbor-tongyuanjuzhen.txt'), row.names=F, col.names=F, quote=F, sep='\t')

###
for(i in 1:nrow(nppjz)){
bh=nppjz[i,]
bh[bh==1]=1.6
bh[bh==0]=-0.3
nppjz[i,]=bh
}
###
nlpjz=ntyjz+nppjz
###
write.table(nlpjz,file= paste0(NETWORK.FOLD3,'neighbor-lianpeijuzhen.txt'), row.names=F, col.names=F, quote=F, sep='\t')
###
matlablj='D:/network_alignment/wcs.exe'
system(matlablj)

###
nlpjgjz<-read.table(paste0(NETWORK.FOLD3,'neighbor-alignment_result.txt'),header=F,sep=',')
###
an1=0
nrppd=list()
nyppd=list()
for(i in 1:length(rmtg)){
npwz=which(nlpjgjz[i,]==1)
  if(length(npwz)>0){
          an1=an1+1
   nrppd[[an1]]=rmtg[i]
   nyppd[[an1]]=ymtg[npwz]
  }
}


###
sink(paste0(NETWORK.FOLD,'new match-pair gene.txt'),append = F)
for(i in 1:an1){ 
      cat(nrppd[[i]])
      cat('\t')
      cat(nyppd[[i]]) 
      cat('\n') 
}  
sink()
###
sink(paste0(NETWORK.FOLD,'all match-pair gene.txt'),append = T)
for(i in 1:an1){ 
      cat(nrppd[[i]])
      cat('\t')
      cat(nyppd[[i]]) 
      cat('\n') 
}  
sink()
 }
}












  
