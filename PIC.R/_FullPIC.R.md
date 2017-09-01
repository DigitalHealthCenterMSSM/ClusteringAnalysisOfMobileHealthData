# PSB Cluster
## Unabridged PIC R Code
```markdown
library("RColorBrewer") #brewer.pal
library("caTools") #runmean

################FUNCTIONS DEFINED

my.impute.p=function(r, p, alpha)
{
  temp=alpha*(1-p)+p
  if(r*temp>1) ### make sure alpha*r<=1
    {
     if(p==1)
     {
       alpha=1      
     } else{
        alpha=(1/r-p)/(1-p)
     }           
    temp=alpha*(1-p)+p
    }     

  p1=p/temp 
  r0=r*temp/alpha

  temp1=(1-alpha*r0)*p+alpha*(1-r0)*(1-p)
  p.hat=(1-alpha*r0)*p/temp1
  if(p.hat<0) p.hat=0  #### due to numeric reasons, sometimes p.hat is given a small negative value 
  return(p.hat)  
}




##performs missing data imputation
my.impute=function(x.v, cen=7,alpha,pa=mypa) ### cen: the censoring period; i.e. if the last day is within the last XX days, the users are deemed as an active user on day X
{
   temp=which(!is.na(x.v))
   st=1
   en=length(x.v)
      if(en>(length(x.v)-cen))
      en=length(x.v)
   r=mean(!is.na(x.v[st:en]))#response rate within enrollment period
      ####calculate p.t based on 1 month neighborhood
 	p.t=runmean(x.v,30)
 	p.t[is.na(p.t)]=rep(mean(x.v, na.rm=T),sum(is.na(p.t)))
   	p=mean(x.v, na.rm=T)	
 	x.v.act=x.v[st:en]
   nan=sum(is.na(x.v.act))
   phat=p.t[st:en][is.na(x.v.act)]
   index=c(st:en)[is.na(x.v.act)]
   if (is.na(alpha)) {alpha=sample(1:2,1,prob=c(1-pa,pa),replace=TRUE)}
   if(alpha==2){
   phat=mapply(my.impute.p,rep(r,length(index)),p.t[index],rep(alpha,length(index)))
 	}
     
   x.v.act[is.na(x.v.act)]=unlist(lapply(phat,function(x) sample(0:1,1,prob=c(1-x,x),replace=TRUE)))
   result=rep(NA, length(x.v))
   result[st:en]=x.v.act
   return(result)
}

##performs k-means clustering on each iteration of 
get_clusters=function(pc.concat,c.n){###k-means
temp=apply(pc.concat, 2, my.stand)
clu.try1=kmeans(temp, centers=c.n,iter.max=100)
return(clu.try1$cluster)
}

### filter to select users without at least ts responses

get_sel=function(data.m,ts=50){
	end.day=max(as.numeric(colnames(data.m)))
data.na=apply(is.na(data.m), 1, sum)
data.m[data.na<=(end.day-ts), ]  
}

#########################

my.stand<-function(x.v)
{
   result=(x.v-mean(x.v))/sd(x.v)
}

################infere complete data using imputation + smoothing

get_PIC=function(data.sel,seed=1,c.n=3, sm.n=mysmooth, t0=25, t1=160, K=20,mypa=0)
{
set.seed(seed)
dataI=data.sel

for(i in 1:nrow(data.sel))
{
   dataI[i,]=my.impute(data.sel[i,],alpha=alpha[i],pa=mypa)

}

dataI.sm=t(apply(dataI,1,function(x) runmean(x,sm.n)))
dataI.sm[,t0:t1]
}

##############


get_multi.consensus=function(data.m,B=100,t0=25,t1=160,sm.n=15,ts=50,mypa=0,clus.nv=c(3:7),mystrat=get_PIC){

data.sel=get_sel(data.m,ts=50)
temp.users=row.names(data.sel)
a.n=length(temp.users)

multi.consen=list(NULL)
multi.new=list(NULL)
for(i in 1:length(clus.nv))
{
 cur.cn=clus.nv[i]
 impu.consen=array(0, dim=c(a.n, a.n, B))
 imp.new=array(NA, dim=c(a.n, length(t0:t1), B))

 for(j in 1:B)
 {
    print(c(i,j))
    
    data.imp=mystrat(data.sel,mystrat,seed=j,t0=25,K=20,t1=160,sm.n=15,mypa=mypa)
    try1=get_clusters(data.imp,c.n=cur.cn)
    temp=unique(try1)
    cur.con=matrix(0, nrow=a.n, ncol=a.n)
    for(h in temp)
     {
       cur.con[try1==h, try1==h]=1       
     }    
     impu.consen[,,j]=cur.con
     imp.new[,,j]=data.imp
      }
multi.consen[[i]]=impu.consen #adjacency matrix
multi.new[[i]]=imp.new #imputed data


}

list(multi.consen=multi.consen,multi.new=multi.new)
}


#############this function 

get_cl=function(multi.consen,data.m,data.n,i=2,scheme="Greens",my.ylab,clus.nv=2:3,ts=50,t0=25,t1=160){
	
data.sel=get_sel(data.m,ts=50)
temp.users=row.names(data.sel)
a.n=length(temp.users)

c.n=clus.nv[i]
colfunc <- colorRampPalette(brewer.pal(9,"YlOrRd"))
clu.micc=t(apply(multi.consen[[i]], c(1,2), mean))
row.names(clu.micc)<-temp.users

#alternative: use hierarchical clustering
# dm=dist(clu.micc)
# temp1=hclust(dm, method="complete")# "mcquitty")# ="average")
# clu.try0=cutree(temp1, c.n)

clu.try0=get_clusters(clu.micc,c.n)  #alternative: use k-means clustering

#first order by table rank
###relabel these so that they are numbered  in order of overall mean 

clu.try1=rank(table(clu.try0))[clu.try0]
names(clu.try1)=temp.users
temp=order(clu.try1)

data.n=apply(data.n[[i]],c(1,2),mean) ###use if using the consensus result
mean.curve=apply(data.n,2,function(x) tapply(x,factor(clu.try1[temp]),mean,na.rm=T))

marg.mean=apply(mean.curve,1,mean,na.rm=T)
mycol.mean<-colorRampPalette(brewer.pal(5,scheme))(c.n)[rank(marg.mean)]

colnames(mean.curve)=colnames(data.sel[t0:t1])
mycol<-mycol.mean[clu.try1[temp]]
names(mycol)=names(clu.try1[temp])

list(clu.micc=clu.micc,mycol=mycol,mean.curve=mean.curve)
}


########INPUTS

#DaySym.m : raw NxT data matrix with missing data
#alpha : vector of length N identifying users with non-ignorable missing 
#B number of imputations
#sm.n smoothing factor
#t0 start date
#t1 end date
#ts filter to select users without at least ts responses
#clus.nv consider several cluster 

#########DATA ANALYSIS


emp.pa=length(alpha[alpha%in%2])/length(alpha[alpha%in%c(1,2)])

#perform consensus clustering based on PIC startegy
mc.PIC_day=get_multi.consensus(DaySym.m,t0=25,t1=160,mystrat=get_PIC,
B=10,sm.n=15,ts=50,mypa=emp.pa,clus.nv=2:3)

#######collect data ouputs
DaySym.n=mc.PIC_day$multi.new #list of m=length(clus.nv) items, where each items is an NxTxB predicted complete data matrux
mc.PIC_day=mc.PIC_day$multi.consen #list of m=length(clus.nv) items, where each item is an NxNXB adjacency matrix


##get consensus cluster assignment 
clust_day=get_cl(mc.PIC_day,DaySym.m,DaySym.n,i=2,"Blues",my.ylab="Day symptoms (PIC)",clus.nv=2:3,ts=50)


#######OUTPUTS


clu.micc=clust_day$clu.micc #consensus adjacency matrix
mycol=clust_day$mycol #cluster assignments color-coded based on the overall mean 
mean.curve=clust_day$mean.curve #predicted mean curve



matplot(t(mean.curve),type="l",col=names(table(mycol)),lwd=3,lty=1)
```
