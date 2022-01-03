#copyright 2021, Auton Systems LLC

library(plyr)
library(ggplot2)
library(data.table)

cv = function(formula,model,predict,data,k=10,blocks=NULL){
  folds = rep(1:k,ceiling(nrow(data)/k))[order(runif(ceiling(nrow(data)/k)*k))][1:nrow(data)]
  if(!is.null(blocks)){
    print("using blocking")
    # assign each block to a fold
    blk = unique(blocks)
    blk_folds = rep(1:k,ceiling(length(blk)/k))[order(runif(ceiling(length(blk)/k)*k))][1:length(blk)]
    # broadcast out to observations
    folds = blk_folds[match(blocks,blk)]
  }
  sdf = rbindlist(lapply(1:k,function(fold){
    print(paste0("fold: ",fold))
    data.frame(i=which(folds==fold),score=predict(model(formula,data[folds!=fold,]), data[folds==fold,]))
  }))
  sdf = sdf[order(sdf$i),]
  return(sdf$score)
}

read_params = function(path,keep = c("initialI","initialNs","initialPrR","beta","nu","lambda","mutCost","epsilon","lambdaAntigenic")){
  f = scan(paste0(path,"/parameters_load.yml"),what=character(),sep='\n',quote=NULL)
  f = strsplit(f,':')
  names = unlist(lapply(f,function(s) s[1]))
  values = unlist(lapply(f,function(s) gsub('\\[|\\]','',s[2])))
  idx = names%in%keep
  values = as.numeric(values[idx])
  names(values) = names[idx]
  return(values)
}

#sims = list.dirs(path = "./sims", recursive = TRUE)
sims = list.dirs(path = "./sims_fixedparam", recursive = TRUE)
#sims = sims[grep("\\./sims/job[0-9]+",sims)]
sims = sims[grep("\\./sims_fixedparam/job[0-9]+",sims)]
sims_finished = sims[unlist(lapply(sims,function(d) "out.antigenicSamples"%in%list.files(d)))]

rfpredict = function(fit,data){
  predict(fit,data,type="prob")[,"TRUE"]
}
glmfit = function(formula,data){
  glm(formula,family=binomial(logit),data=data)
}
glmpredict = function(fit,data){
  predict(fit,data)
}


# for each finished sim, count number of antigenic variants
prop = 0.05
dur = 3*28

counts = rbindlist(lapply(sims_finished,function(d){
  params = read_params(d)
  p = read.csv(paste0(d,"/out.infectionsByPhenotype"))
  p = ddply(p,"type",summarize,day=day,I=cumsum(infections-recoveries))
  p = join(p,ddply(p,"day",summarize,totI=sum(I)))
  
  prevalence = ddply(p,"type",summarize,mx=max(I/totI),cnt=sum(I/totI > prop),first=min(day))
  prevalence = prevalence[prevalence$first< 1095-365,] # don't consider tail end of simulation
  transient = prevalence[prevalence$cnt>0 & prevalence$cnt<dur,]
  persistant = prevalence[prevalence$cnt>=dur,]
  retdf = data.frame(as.list(params))
  retdf$sim = d
  retdf$transient = nrow(transient)
  retdf$persistant = nrow(persistant)
  retdf$mx = NA
  if(length(unique(prevalence$type))>1) retdf$mx = max(prevalence$mx[prevalence$type!=0])
  return( retdf )
}))

##########################################
library(ROCR)
library(randomForest)
df = counts[,c("initialNs","initialI","initialPrR","beta","nu","lambda","mutCost","epsilon","lambdaAntigenic")]
df$y = as.factor( counts$transient>0 | counts$persistant>1 )
scores = cv(y~.,randomForest,rfpredict,df)
perf <- performance(prediction(scores, df$y=="TRUE"),"tpr","fpr")
plot(perf)
abline(0,1)
############################################


fitness_features = rbindlist(lapply(sims_finished,function(d){
  p = read.csv(paste0(d,"/out.infectionsByPhenotype"))
  p = ddply(p,"type",summarize,day=day,I=cumsum(infections-recoveries))
  p = join(p,ddply(p,"day",summarize,totI=sum(I)))

  vF = read.delim(paste0(d,"/out.viralFitnessSeries"))
  vF$date = floor((vF$date)*365)
  vF2 = vF[!is.na(vF$type),]
  vF_ = vF[is.na(vF$type),]
  
  prevalence = ddply(p,"type",summarize,cnt=sum(I/totI > prop),birth=min(day),first=min(c(max(day),day[I/totI > prop])),last=max(c(-1,day[I/totI > prop])))
  prevalence = prevalence[prevalence$first< 1095-365,] # don't consider tail end of simulation
  transient = prevalence[prevalence$cnt>0 & prevalence$cnt<dur,]
  persistant = prevalence[prevalence$cnt>=dur,]
  
  types = rbind(transient,persistant)
  types = types[types$type!=0,]
  
  vF_$y = FALSE
  if(nrow(types)>0){
  for(i in 1:nrow(types)){
    vF_$y[vF_$date>=types$birth[i] & vF_$date<=types$first[i]] = TRUE 
  }
    vF_ = vF_[vF_$date<=max(types$first),]
    }
  vF_ = vF_[vF_$date< 1095-365,]
  vF_$type=NULL
  vF_$date=NULL
  vF_$block = d
  
  return( vF_ )
}))
fitness_features$y = as.factor(fitness_features$y)

#df2 = fitness_features
#blocks = df2$block
#df2$block = NULL
#scores = cv(y~.,randomForest,rfpredict,df2,blocks = blocks)
#perf <- performance(prediction(scores, df2$y=="TRUE"),"tpr","fpr")
#plot(perf)
#abline(0,1)

df2 = rbind(fitness_features[fitness_features$y=="TRUE",],
            fitness_features[sample(which(fitness_features$y=="FALSE"),sum(fitness_features$y=="TRUE")),])
blocks = df2$block
df2$block = NULL
scores = cv(y~.,randomForest,rfpredict,df2,blocks=blocks)
perf <- performance(prediction(scores, df2$y=="TRUE"),"tpr","fpr")
plot(perf)
abline(0,1)

df2 = rbind(fitness_features[fitness_features$y=="TRUE",],
            fitness_features[sample(which(fitness_features$y=="FALSE" & fitness_features$varSigma>0),sum(fitness_features$y=="TRUE")),])
blocks = df2$block
df2$block = NULL
scores = cv(y~.,glmfit,glmpredict,df2,blocks=blocks)
slot(performance(prediction(scores, df2$y=="TRUE"),"auc"),'y.values')[[1]]
perf <- performance(prediction(scores, df2$y=="TRUE"),"tpr","fpr")
plot(perf)
abline(0,1)
summary(glm(y~.,family=binomial(logit),data=df2))

ggplot(fitness_features)+geom_density(aes(x=varSigma,color=y))


############################################


cumulative_features = rbindlist(lapply(sims_finished,function(d){
  params = read_params(d)
  p = read.csv(paste0(d,"/out.infectionsByPhenotype"))
  p = ddply(p,"type",summarize,day=day,I=cumsum(infections-recoveries),totalInfections=cumsum(infections))
  p = join(p,ddply(p,"day",summarize,totI=sum(I)))
  
  p2 = ddply(p,"day",summarize,totI=sum(I),totInfections=sum(totalInfections))
  SI = p2$totI*(params["initialNs"]-p2$totInfections)
  vF_ = data.frame(date=p2$day,SI=cumsum(SI))

  prevalence = ddply(p,"type",summarize,cnt=sum(I/totI > prop),birth=min(day),first=min(c(max(day),day[I/totI > prop])),last=max(c(-1,day[I/totI > prop])))
  prevalence = prevalence[prevalence$first< 1095-365,] # don't consider tail end of simulation
  transient = prevalence[prevalence$cnt>0 & prevalence$cnt<dur,]
  persistant = prevalence[prevalence$cnt>=dur,]
  
  types = rbind(transient,persistant)
  types = types[types$type!=0,]
  
  vF_$y = FALSE
  if(nrow(types)>0){
    for(i in 1:nrow(types)){
      vF_$y[vF_$date>=types$birth[i] & vF_$date<=types$first[i]] = TRUE 
    }
    vF_ = vF_[vF_$date<=max(types$first),]
    }
  vF_ = vF_[vF_$date< 1095-365,]
  vF_$type=NULL
  vF_$date=NULL
  vF_$block = d
  
  return( vF_ )
}))
cumulative_features$y = as.factor(cumulative_features$y)
df2 = cumulative_features
blocks = df2$block
df2$block = NULL
scores = cv(y~.,randomForest,rfpredict,df2,blocks = blocks)
perf <- performance(prediction(scores, df2$y=="TRUE"),"tpr","fpr")
plot(perf)
abline(0,1)

df2 = rbind(cumulative_features[cumulative_features$y=="TRUE",],
            cumulative_features[sample(which(cumulative_features$y=="FALSE"),sum(cumulative_features$y=="TRUE")),])
blocks = df2$block
df2$block = NULL
scores = cv(y~.,randomForest,rfpredict,df2,blocks=blocks)
perf <- performance(prediction(scores, df2$y=="TRUE"),"tpr","fpr")
plot(perf)
abline(0,1)

df2 = rbind(cumulative_features[cumulative_features$y=="TRUE",],
            cumulative_features[sample(which(cumulative_features$y=="FALSE"),sum(cumulative_features$y=="TRUE")),])
blocks = df2$block
df2$block = NULL
scores = cv(y~.,glmfit,glmpredict,df2,blocks=blocks)
slot(performance(prediction(scores, df2$y=="TRUE"),"auc"),'y.values')[[1]]
perf <- performance(prediction(scores, df2$y=="TRUE"),"tpr","fpr")
plot(perf)
abline(0,1)
summary(glm(y~.,family=binomial(logit),data=df2))

############################################

features = rbindlist(lapply(sims_finished,function(d){
  params = read_params(d)
  p = read.csv(paste0(d,"/out.infectionsByPhenotype"))
  p = ddply(p,"type",summarize,day=day,I=cumsum(infections-recoveries),totalInfections=cumsum(infections))
  p = join(p,ddply(p,"day",summarize,totI=sum(I)))
  
  vF = read.delim(paste0(d,"/out.viralFitnessSeries"))
  vF$date = floor((vF$date)*365)
  vF2 = vF[!is.na(vF$type),]
  vF_ = vF[is.na(vF$type),]
  
  p2 = ddply(p,"day",summarize,totI=sum(I),totInfections=sum(totalInfections))
  SI = p2$totI*(params["initialNs"]-p2$totInfections)
  vF_$SI = cumsum(SI)
  vF_$totI = p2$totI
  vF_$S = params["initialNs"]-p2$totInfections
  
  
  prevalence = ddply(p,"type",summarize,cnt=sum(I/totI > prop),birth=min(day),first=min(c(max(day),day[I/totI > prop])),last=max(c(-1,day[I/totI > prop])))
  prevalence = prevalence[prevalence$first< 1095-365,] # don't consider tail end of simulation
  transient = prevalence[prevalence$cnt>0 & prevalence$cnt<dur,]
  persistant = prevalence[prevalence$cnt>=dur,]
  
  types = rbind(transient,persistant)
  types = types[types$type!=0,]
  
  vF_$y = FALSE
  if(nrow(types)>0){
    for(i in 1:nrow(types)){
      vF_$y[vF_$date>=types$birth[i] & vF_$date<=types$first[i]] = TRUE 
    }
    vF_ = vF_[vF_$date<=min(types$first),]
  }
  vF_ = vF_[vF_$date< 1095-365,]
  vF_$type=NULL
  vF_$date=NULL
  vF_$block = d
  
  return( vF_ )
}))
features$y = as.factor(features$y)
features_orig = features
features = features[features$varSigma>0,]

df2 = features
blocks = df2$block
df2$block = NULL
scores = cv(y~.,randomForest,rfpredict,df2,blocks = blocks)
perf <- performance(prediction(scores, df2$y=="TRUE"),"tpr","fpr")
plot(perf)
abline(0,1)

df2 = rbind(features[features$y=="TRUE",],
            features[sample(which(features$y=="FALSE"),sum(features$y=="TRUE")),])
blocks = df2$block
df2$block = NULL
scores = cv(y~SI+totI+S+totIoS+d1_SI+d1_totI+d1_S,randomForest,rfpredict,df2,blocks=blocks)
slot(performance(prediction(scores, df2$y=="TRUE"),"auc"),'y.values')[[1]]
perf <- performance(prediction(scores, df2$y=="TRUE"),"tpr","fpr")
plot(perf)
abline(0,1)

df2 = rbind(features[features$y=="TRUE",],
            features[sample(which(features$y=="FALSE"),sum(features$y=="TRUE")),])
blocks = df2$block
df2$block = NULL
scores = cv(y~covBetaSigma,glmfit,glmpredict,df2,blocks=blocks)
slot(performance(prediction(scores, df2$y=="TRUE"),"auc"),'y.values')[[1]]
perf <- performance(prediction(scores, df2$y=="TRUE"),"tpr","fpr")
plot(perf)
abline(0,1)
summary(glm(y~.,family=binomial(logit),data=df2))

ggplot(features)+geom_density(aes(x=meanR,color=y=="TRUE"))



############################################

features = rbindlist(lapply(sims,function(d){
  params = read_params(d)
  p = read.csv(paste0(d,"/out.infectionsByPhenotype"))
  p = ddply(p,"type",summarize,day=day,I=cumsum(infections-recoveries),totalInfections=cumsum(infections))
  p = join(p,ddply(p,"day",summarize,totI=sum(I)))
  
  vF = read.delim(paste0(d,"/out.viralFitnessSeries"))
  vF$date = floor((vF$date)*365)
  vF2 = vF[!is.na(vF$type),]
  vF_ = vF[is.na(vF$type),]
  
  p2 = ddply(p,"day",summarize,totI=sum(I),totInfections=sum(totalInfections))
  SI = p2$totI*(params["initialNs"]-p2$totInfections)
  vF_$SI = cumsum(SI)
  vF_$totI = p2$totI
  vF_$N = params["initialNs"]
  vF_$S = params["initialNs"]-p2$totInfections
  vF_$S[vF_$S<0] = 0
  vF_$totIoS = vF_$totI/(vF_$S+1)
  
  
  prevalence = ddply(p,"type",summarize,cnt=sum(I/totI > prop),birth=min(day),first=min(c(max(day),day[I/totI > prop])),last=max(c(-1,day[I/totI > prop])))
  prevalence = prevalence[prevalence$first< 1095-365,] # don't consider tail end of simulation
  transient = prevalence[prevalence$cnt>0 & prevalence$cnt<dur,]
  persistant = prevalence[prevalence$cnt>=dur,]
  
  types = rbind(transient,persistant)
  types = types[types$type!=0,]
  
  if(nrow(vF_)<4) return(data.frame())
  delay1 = vF_[1:(nrow(vF_)-1),c("SI","totI","S","meanR","varR","meanBeta","varBeta","meanSigma","varSigma","covBetaSigma")]
  colnames(delay1) = paste("d1",colnames(delay1),sep="_")
  #delay2 = vF_[1:(nrow(vF_)-2),c("SI","totI","S","meanR","varR","meanBeta","varBeta","meanSigma","varSigma","covBetaSigma")]
  #colnames(delay2) = paste("d2",colnames(delay2),sep="_")
  vF_ = cbind(vF_[2:nrow(vF_),],delay1)
  #vF_[,colnames(delay1)] = vF_[,c("SI","totI","S","meanR","varR","meanBeta","varBeta","meanSigma","varSigma","covBetaSigma")]/(1e-6+vF_[,colnames(delay1)])
  #vF_ = cbind(vF_,delay2)
  
  vF_$y = FALSE
  if(nrow(types)>0){
    for(i in 1:nrow(types)){
      vF_$y[vF_$date>=types$birth[i] & vF_$date<=types$first[i]] = TRUE 
    }
    vF_ = vF_[vF_$date<=min(types$first),]
  }
  vF_ = vF_[vF_$date< 1095-365,]
  if(nrow(vF_)<1) return(data.frame())
  vF_$type=NULL
  vF_$date=NULL
  vF_$block = d
  
  return( vF_ )
}))
features$y = as.factor(features$y)
features_orig = features
features = features[features$varSigma>0,]

df2 = features
blocks = df2$block
df2$block = NULL
scores = cv(y~.,randomForest,rfpredict,df2,blocks = blocks)
perf <- performance(prediction(scores, df2$y=="TRUE"),"tpr","fpr")
plot(perf)
abline(0,1)

z=1.96
n = sum(df2$y=="FALSE")
roc$fpr2=(roc$fpr+0.5*z^2/n)/(1+z^2/n)
roc$fpr_ci = z/(1+z^2/n)*sqrt(roc$fpr*(1-roc$fpr)/n+z^2/4/n^2)
n = sum(df2$y=="TRUE")
roc$tpr2=(roc$tpr+0.5*z^2/n)/(1+z^2/n)
roc$tpr_ci = z/(1+z^2/n)*sqrt(roc$tpr*(1-roc$tpr)/n+z^2/4/n^2)
tmp = data.frame(x=roc$fpr2-roc$fpr_ci, y=roc$tpr2+roc$tpr_ci)
tmp = ddply(tmp,"x",summarize,y=max(y))
roc$u = approx(tmp$x,tmp$y,roc$fpr2)$y
tmp = data.frame(x=roc$fpr2+roc$fpr_ci, y=roc$tpr2-roc$tpr_ci)
tmp = ddply(tmp,"x",summarize,y=min(y))
roc$l = approx(tmp$x,tmp$y,roc$fpr2)$y
roc$l[is.na(roc$l)]=0
roc$u[is.na(roc$u)]=1
write.csv(roc,"new_variant_roc.csv",row.names=F)
p=ggplot(roc)+geom_ribbon(aes(x=fpr2,ymin=l,ymax=u),alpha=0.33)+
  geom_line(aes(x=fpr2,y=tpr2))+
  geom_line(aes(x=fpr2,y=fpr2),lty=2)+
  xlab("False positive rate")+ylab("True postive rate")+
  theme_minimal(14)+theme(aspect.ratio=1)
plot(p)
ggsave("new_variant_roc.png",p)
p2 = p+theme(panel.grid.minor = element_blank())+
       scale_x_log10(breaks = 10^seq(-3,0,1),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)))+
       annotation_logticks()
plot(p2)
ggsave("new_variant_roc_log.png",p2)
p2 = p+theme(panel.grid.minor = element_blank())+
  scale_x_log10(breaks = 10^seq(-3,0,1),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  annotation_logticks()
plot(p2)

df2 = rbind(features[features$y=="TRUE",],
            features[sample(which(features$y=="FALSE"),sum(features$y=="TRUE")),])
blocks = df2$block
df2$block = NULL
scores = cv(y~N+totI+S+SI+d1_totI+d1_S+d1_SI,randomForest,rfpredict,df2,blocks=blocks)
slot(performance(prediction(scores, df2$y=="TRUE"),"auc"),'y.values')[[1]]
perf <- performance(prediction(scores, df2$y=="TRUE"),"tpr","fpr")
plot(perf)
abline(0,1)

n=2000
roc2 = data.frame(fpr=slot(perf,'x.values')[[1]],tpr=slot(perf,'y.values')[[1]])
z=1.96
roc2$fpr2=(roc2$fpr+0.5*z^2/n)/(1+z^2/n)
roc2$fpr_ci = z/(1+z^2/n)*sqrt(roc2$fpr*(1-roc2$fpr)/n+z^2/4/n^2)
roc2$tpr2=(roc2$tpr+0.5*z^2/n)/(1+z^2/n)
roc2$tpr_ci = z/(1+z^2/n)*sqrt(roc2$tpr*(1-roc2$tpr)/n+z^2/4/n^2)
tmp = data.frame(x=roc2$fpr2-roc2$fpr_ci, y=roc2$tpr2+roc2$tpr_ci)
tmp = ddply(tmp,"x",summarize,y=max(y))
roc2$u = approx(tmp$x,tmp$y,roc2$fpr2)$y
tmp = data.frame(x=roc2$fpr2+roc2$fpr_ci, y=roc2$tpr2-roc2$tpr_ci)
tmp = ddply(tmp,"x",summarize,y=min(y))
roc2$l = approx(tmp$x,tmp$y,roc2$fpr2)$y
roc2$l[is.na(roc2$l)]=0
roc2$u[is.na(roc2$u)]=1
write.csv(roc2,"new_variant_roc2.csv",row.names=F)
p=ggplot(roc2)+geom_ribbon(aes(x=fpr2,ymin=l,ymax=u),alpha=0.33)+
  geom_line(aes(x=fpr2,y=tpr2))+
  geom_line(aes(x=fpr2,y=fpr2),lty=2)+
  xlab("False positive rate")+ylab("True postive rate")+
  theme_minimal(14)+theme(aspect.ratio=1)
plot(p)
ggsave("new_variant_roc2.png",p)
p2 = p+theme(panel.grid.minor = element_blank())+
  scale_x_log10(breaks = 10^seq(-3,0,1),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  annotation_logticks()
plot(p2)
ggsave("new_variant_roc2_log.png",p2)

tmp = roc; tmp$Features = "All"; roc3=tmp
tmp = roc2; tmp$Features = "Cases"; roc3=rbind(roc3,tmp)
p=ggplot(roc3)+geom_ribbon(aes(x=fpr2,ymin=l,ymax=u,fill=Features),alpha=0.33)+
  geom_line(aes(x=fpr2,y=tpr2,color=Features))+
  geom_line(aes(x=fpr2,y=fpr2),lty=3)+
  xlab("False positive rate")+ylab("True postive rate")+
  theme_minimal(14)+theme(aspect.ratio=1)
plot(p)
ggsave("new_variant_roc_both.png",p)
p2 = p+theme(panel.grid.minor = element_blank())+
  scale_x_log10(breaks = 10^seq(-3,0,1),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  annotation_logticks()
plot(p2)
ggsave("new_variant_roc_both_log.png",p2)



df2 = rbind(features[features$y=="TRUE",],
            features[sample(which(features$y=="FALSE"),sum(features$y=="TRUE")),])
blocks = df2$block
df2$block = NULL
scores = cv(y~.,glmfit,glmpredict,df2,blocks=blocks)
slot(performance(prediction(scores, df2$y=="TRUE"),"auc"),'y.values')[[1]]
perf <- performance(prediction(scores, df2$y=="TRUE"),"tpr","fpr")
plot(perf)
abline(0,1)
summary(glm(y~.,family=binomial(logit),data=df2))

ggplot(features)+geom_density(aes(x=meanR,color=y=="TRUE"))

##################################################
library(RColorBrewer)
d = "./sims/job145"
params = read_params(d)
p = read.csv(paste0(d,"/out.infectionsByPhenotype"))
p = ddply(p,"type",summarize,day=day,I=cumsum(infections-recoveries),totalInfections=cumsum(infections))
p = join(p,ddply(p,"day",summarize,totI=sum(I)))
types = unique(p$type); newtypes = paste0(runif(length(types)),types)
p$type = newtypes[match(p$type,types)]
vF = read.delim(paste0(d,"/out.viralFitnessSeries"))
vF$date = floor((vF$date)*365)
vF2 = vF[!is.na(vF$type),]
vF_ = vF[is.na(vF$type),]

prevalence = ddply(p,"type",summarize,cnt=sum(I/totI > prop),birth=min(day),first=min(c(max(day),day[I/totI > prop])),last=max(c(-1,day[I/totI > prop])))
prevalence = prevalence[prevalence$first< 1095-365,] # don't consider tail end of simulation
transient = prevalence[prevalence$cnt>0 & prevalence$cnt<dur,]
persistant = prevalence[prevalence$cnt>=dur,]
types = rbind(transient,persistant)

pg=ggplot(p)+geom_area(aes(x=day,y=I,fill=as.factor(type)),size=0.1,color='black')+
  theme_minimal(14)+theme(legend.position='none')+xlab("Day")+ylab("Infections")
plot(pg)
ggsave("Example1.png",pg)
pg=ggplot(p)+geom_area(aes(x=day,y=I/totI,fill=as.factor(type)),size=0.1,color="black")+
  theme_minimal(14)+theme(legend.position='none')+xlab("Day")+ylab("Infections")
plot(pg)
ggsave("Example1_prop.png",pg)
library(reshape)
tmp = melt(vF_,id.vars=c("date"))
tmp = tmp[tmp$variable%in%c("varSigma","varBeta","varR"),]
colnames(tmp) = c("date","Variable","value")
tmp$Variable = as.character(tmp$Variable)
tmp$Variable[tmp$Variable=="varSigma"] = "Susc. var."
tmp$Variable[tmp$Variable=="varBeta"] = "Beta var."
tmp$Variable[tmp$Variable=="varR"] = "R var."
pg=ggplot(tmp)+
  geom_line(aes(x=date,y=value,color=Variable))+
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))+theme_minimal(14)+ylab("Value")+xlab("Day")
plot(pg)
ggsave("Example1_var.png",pg)
