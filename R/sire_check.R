# Parameters: 1. genotest: SNP matrix (for test)
#             2. fullgeno: all geno data
#             3. metad: meta data (for test)
#             4. fullmeta: all meta 
#             5. resu: results calculated previously
sire_check<-function(genotest=data,fullgeno=add_geno,metad=meta_data,fullmeta=add_meta, resu=results,parent_id="Sire",joint=F){
  
  
  idx<-intersect(rownames(resu),rownames(metad)) # Overlap of our meta and meta_all
  metad<-metad[idx,] # Common part
  
  fullgeno<-fullgeno[,colnames(fullgeno)%in%colnames(genotest)] # Common SNPs in our geno test
  
  joint_geno<-t(rbind(fullgeno,genotest[,colnames(fullgeno)])) # All geno(ref panel) + geno data(new data set)
  
  
  prs<-parent_id # pass sire information to prs
  parent_id<-parent_id[1]
  
  fullmeta<-rbind(fullmeta,(metad[,colnames(fullmeta)])) # All meta (Ref) + meta data (new data set)
  all_RN<-paste(fullmeta$BREED,fullmeta$REGISTRATION.NUMBER,sep="_") # all_RN=Breed_REG.Number
  
  bsase_dt<-metad[,c("BARCODE","BREED",parent_id)] # Keep barcode, breed, and sire only from new meta
  colnames(bsase_dt)<-c("BARCODE","BREED",parent_id)
  
  
  meta_sire_RN<-paste(metad$BREED,metad[,parent_id],sep="_") # meta_sire_RN=breed_sire
  names(meta_sire_RN)<-metad$BARCODE # Search them by Barcode
  
  
  srs<-sapply(meta_sire_RN,grep,x = all_RN) # To find REG number in all_meta (by row.number)
  srs<-unlist(srs) ##(Try something to indicate here) 
  
  if(length(srs)==0)
  {
    stop("No parent's information matched")
  }
  
  idx<-srs[metad$BARCODE] # sire (dam, or both) in new_meta
  sire_dt<-fullmeta[idx,c("BREED","BARCODE","REGISTRATION.NUMBER","CASE.NUMBER")] # select these 4 columns
  colnames(sire_dt)<-paste("PARENT",colnames(sire_dt),sep = ".") 
  checkp<-cbind(bsase_dt,sire_dt) # combine new_meta and sire info
  
  checkp<-checkp[!is.na(checkp$PARENT.REGISTRATION.NUMBER),] # subjects with valid info only
  
  
  
  #sum(rownames(genotest)%in%rownames(fullgeno))
  percincomp<-rep(0,nrow(checkp))                  
  
  boarKIT<-sireKIT<-rep("",nrow(checkp))
  
  
  barc<-checkp$PARENT.BARCODE%in%colnames(joint_geno)
  regn<-checkp$PARENT.REGISTRATION.NUMBER%in%colnames(joint_geno)
  csn<-checkp$PARENT.CASE.NUMBER%in%colnames(joint_geno)
  
  # whonotmatch <- rep(FALSE, nrow(checkp))
  sname<-checkp$PARENT.BARCODE
  for (i in 1:nrow(checkp)){
    g1<-joint_geno[,checkp$BARCODE[i]]
    if (!barc[i]){
      if(!regn[i]){
        if(!csn[i]){
          # whonotmatch[i]=T
          print(checkp[i,])
          print("Can't find sire in genodata")
          next
        }
        sname[i]<-checkp$PARENT.CASE.NUMBER[i]
      }
      sname[i]<-checkp$PARENT.REGISTRATION.NUMBER[i]
    }
    
    g2<-joint_geno[,sname[i]]
    
    disc<-sum(abs(g1-g2)==2,na.rm = T)
    total<-sum(!is.na(g1-g2))
    percincomp[i]<-(disc/total*100)
    boarKIT[i]<-gsub("NA",".",paste(g1[breedTools:::kit_snps],sep="",collapse=""))
    sireKIT[i]<-gsub("NA",".",paste(g2[breedTools:::kit_snps],sep="",collapse=""))
  }
  boarKIT<-paste("|",boarKIT,"|",sep="")
  parentKIT<-paste("|",sireKIT,"|",sep="")
  checkp<-cbind(checkp,sname)
  Parent_report<-cbind(checkp,percincomp,boarKIT,parentKIT)
  
  
  if(length(prs)==2){
    parent_id<-prs[2]
    checkp1<-checkp
    
    bsase_dt<-metad[,c("BARCODE","BREED",parent_id)]
    colnames(bsase_dt)<-c("BARCODE","BREED","Parent")
    meta_sire_RN<-paste(metad$BREED,metad[,parent_id],sep="_")
    names(meta_sire_RN)<-metad$BARCODE
    
    
    srs<-sapply(meta_sire_RN,grep,x = all_RN)
    srs<-unlist(srs)
    
    
    if(length(srs)==0)
    {
      stop("No parent's information matched")
    }
    
    idx<-srs[metad$BARCODE]
    sire_dt<-fullmeta[idx,c("BREED","BARCODE","REGISTRATION.NUMBER","CASE.NUMBER")]
    colnames(sire_dt)<-paste("PARENT",colnames(sire_dt),sep = ".")
    checkp<-cbind(bsase_dt,sire_dt)
    
    checkp<-checkp[!is.na(checkp$PARENT.REGISTRATION.NUMBER),]
    
    
    
    #sum(rownames(genotest)%in%rownames(fullgeno))
    
    percincomp<-rep(0,nrow(checkp))                  
    
    boarKIT<-sireKIT<-rep("",nrow(checkp))
    
    
    barc<-checkp$PARENT.BARCODE%in%colnames(joint_geno)
    regn<-checkp$PARENT.REGISTRATION.NUMBER%in%colnames(joint_geno)
    csn<-checkp$PARENT.CASE.NUMBER%in%colnames(joint_geno)
    
    sname<-checkp$PARENT.BARCODE
    for (i in 1:nrow(checkp)){
      g1<-joint_geno[,checkp$BARCODE[i]]
      if (!barc[i]){
        if(!regn[i]){
          if(!csn[i]){
            print(checkp[i,])
            print("Can't find sire in genodata")
            next
          }
          sname[i]<-checkp$PARENT.CASE.NUMBER[i]
        }
        sname[i]<-checkp$PARENT.REGISTRATION.NUMBER[i]
      }
      
      g2<-joint_geno[,sname[i]]
      
      disc<-sum(abs(g1-g2)==2,na.rm = T)
      total<-sum(!is.na(g1-g2))
      percincomp[i]<-(disc/total*100)
      boarKIT[i]<-gsub("NA",".",paste(g1[breedTools:::kit_snps],sep="",collapse=""))
      sireKIT[i]<-gsub("NA",".",paste(g2[breedTools:::kit_snps],sep="",collapse=""))
    }
    boarKIT<-paste("|",boarKIT,"|",sep="")
    parentKIT<-paste("|",sireKIT,"|",sep="")
    checkp<-cbind(checkp,sname)
    Parent_report2<-cbind(checkp,percincomp,boarKIT,parentKIT)
    Parent_report<-rbind(Parent_report,Parent_report2)
    
  }
  if(joint&length(prs)==2){
    bothp<-intersect(rownames(checkp),rownames(checkp1))
    if (length(bothp)>0){
      bps<-checkp[bothp,,drop=F]
      bpd<-checkp1[bothp,,drop=F]
      Parent1<-as.character(bps$sname)
      Parent2<-as.character(bpd$sname)
      percer<-rep(0,nrow(bps))
      for (i in 1:nrow(bps)){
        gp<-joint_geno[,bps$BARCODE[i]]
        g1<-joint_geno[,Parent1[i]]  
        g2<-joint_geno[,Parent2[i]]
        err<-((g1==0)&(g2==0)&(gp!=0))+
          ((g1==2)&(g2==2)&(gp!=2))+
          ((g1==0)&(g2==2)&(gp!=1))+
          ((g1==2)&(g2==0)&(gp!=1))+
          ((g1==0)&(g2==1)&(gp==2))+
          ((g1==1)&(g2==0)&(gp==2))+
          ((g1==2)&(g2==1)&(gp==0))+
          ((g1==1)&(g2==2)&(gp==0))
        percer[i]<-100*sum(err,na.rm = T)/sum(!is.na(err))
      }
      
      Parent_report<-data.frame(BARCODE=bps$BARCODE,BREED=bps$BREED,Parent1,Parent2,percer)
    }
  }
  Parent_report <- Parent_report[-grep("sname", colnames(Parent_report))]
  Parent_report <- Parent_report[-grep("PARENT.CASE.NUMBER", colnames(Parent_report))]
  return(Parent_report)
  #write.csv(Parent_report,file="parent.csv")
}
