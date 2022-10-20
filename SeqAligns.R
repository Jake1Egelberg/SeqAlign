

#**********************************************

#--------------------INSTALL NECESSARY PACKAGES

#**********************************************

tryCatch(find.package("utils"),error=function(e){
  install.packages(utils,repos="http://lib.stat.cmu.edu/R/CRAN/") 
})
library(utils)

open_prog<-winProgressBar(title="SeqAlign Super Pro 9000",
                          label="Installing packages...",min=0,max=100,width=300,initial=0)

package.list<-c("this.path",
                "tcltk",
                "stringr",
                "BiocManager")
n_pack<-length(package.list)
for(i in 1:length(package.list)){
  tryCatch(find.package(package.list[i]),
           error = function(e){
             setWinProgressBar(open_prog,value=0,label=paste("Installing ",package.list[i],"...",sep=""))
             install.packages(package.list[i],repos="http://lib.stat.cmu.edu/R/CRAN/") 
           })
}

setWinProgressBar(open_prog,value=25,label="Loading this.path...")
library(this.path)
setWinProgressBar(open_prog,value=35,label="Loading tcltk...")
library(tcltk)
setWinProgressBar(open_prog,value=45,label="Loading stringr...")
library(stringr)
setWinProgressBar(open_prog,value=55,label="Loading BiocManager...")
library(BiocManager)

#Install BiocManager packages
setWinProgressBar(open_prog,value=65,label="Installing BiocManager packages...")
bioc.packages<-c("Rsubread")
for(i in 1:length(bioc.packages)){
  tryCatch(find.package(bioc.packages[i]),
           error=function(e) {
             setWinProgressBar(open_prog,value=65,label=paste("Installing ",bioc.packages[i],"...",sep=""))
             BiocManager::install(bioc.packages[i])
           })
}

setWinProgressBar(open_prog,value=75,label="Loading Rsubread...")
library(Rsubread)

close(open_prog)

#**********************************************

#--------------------------------START ANALYSIS

#**********************************************

.GlobalEnv$dir<-this.dir()

align_fun<-function(){
  
  if(is.na(ref_seq)==FALSE&&length(sel_reads)>=1&&is.na(sel_reads)==FALSE){
  
    prog_inc<-100/(length(sel_reads)+1)
    prog_val<-0
    prog<-winProgressBar(title="SeqAlign Super Pro 9000",
                         label=paste("Aligning...", sep=""),
                         min=0,max=100,width=300,initial=prog_val)
    
    alignment_output<-data.frame()
    for(i in 1:length(sel_reads)){
      .GlobalEnv$x<-sel_reads[i]
      
      prog_val<-prog_val+prog_inc
      setWinProgressBar(prog,value=prog_val,label=paste("Aligning ",gsub(paste(exp_dir,"/",sep=""),"",x),"...", sep=""))
      
      aligned_reads<-align(index=paste(exp_dir,"/refseq",sep=""),readfile1=x)
      df<-data.frame(Read=gsub(".subread.BAM","",colnames(aligned_reads)),
                     Total_reads=aligned_reads[1,],
                     Mapped_reads=aligned_reads[2,],
                     Uniquely_mapped_reads=aligned_reads[3,],
                     Multi_mapping_reads=aligned_reads[4,],
                     Unmapped_reads=aligned_reads[5,],
                     Indels=aligned_reads[6,])
      alignment_output<-rbind(alignment_output,df)
    }
    .GlobalEnv$alignment_output<-alignment_output
    
    #Save alignment output
    setwd(exp_dir)
    write.csv(alignment_output,"AlignmentOutput.csv",row.names=FALSE)
    
    close(prog)
    tk_messageBox(message="Alignment complete and saved to experiment directory!")
    
    edit(alignment_output)
    
  } else{
    tk_messageBox(message="Ensure you have selected one reference sequence and at least one read to align!")
  }
  
}

open_exp_fun<-function(){
  .GlobalEnv$exp_opt<-tk_select.list(choices=c("Create new experiment","Open existing experiment"),multiple=FALSE)
  
  if(exp_opt=="Create new experiment"||exp_opt==""){
    
    date_time<-str_replace_all(Sys.time(),"[:punct:]| ","-")
    tkconfigure(exp_display,text=date_time)
    
    .GlobalEnv$exp_dir<-paste(dir,"/",date_time,sep="")
    dir.create(exp_dir)
    
  } else if(exp_opt=="Open existing experiment"){
    
    cur_dirs<-list.dirs(dir,recursive = FALSE)[-which(list.dirs(dir,recursive = FALSE,full.names = FALSE)=="R-4.2.1")]
    cur_dirs_cur<-gsub(paste(dir,"/",sep=""),"",cur_dirs)
    .GlobalEnv$selected_exp<-tk_select.list(choices=cur_dirs_cur,multiple=FALSE)
    .GlobalEnv$exp_dir<-paste(dir,"/",selected_exp,sep="")
    tkconfigure(exp_display,text=selected_exp)
    
    #Check for parms
    if(file.exists(paste(exp_dir,"/Parms.txt",sep=""))==TRUE){
      
      .GlobalEnv$parms<-read.table(paste(exp_dir,"/Parms.txt",sep=""),sep=",")
      print(parms$RefSeq)
      print(str_split(parms$SelReads," ")[[1]])
      
      
    } else{
      .GlobalEnv$parms<-data.frame(SelReads=NA,
                                   RefSeq=NA)
      setwd(exp_dir)
      write.table(parms,file="Parms.txt",quote=FALSE,sep=",",eol="\n")
    }
    
    .GlobalEnv$sel_reads<-str_split(parms$SelReads," ")[[1]]
    .GlobalEnv$ref_seq<-parms$RefSeq
    
    #Update display
    display_ref()
    display_reads()
    
  }
  
  tkconfigure(sel_ref_seq,state="normal")
  tkconfigure(sel_reads_var,state="normal")
  tkconfigure(align_reads,state="normal")
  
  print(exp_dir)
  
}

display_ref<-function(){
  ref_seq_name<-gsub(paste(exp_dir,"/",sep=""),"",ref_seq)
  tkconfigure(ref_display,text=ref_seq_name)
}

select_ref<-function(){
  .GlobalEnv$ref_seq<-tk_choose.files(default=paste(exp_dir,"/REFERENCE_SEQUENCE",sep=""),multi=FALSE,
                                      filters=matrix(c("FASTA",".fa","TXT",".txt","All files","*"),2,2,byrow=TRUE))
  display_ref()
  
  #Build index for reference sequence
  if(file.exists(paste(exp_dir,"/refseq.reads",sep=""))==FALSE){
    tk_messageBox(message="Reference index not detected. Builing index now.")
    
    prog<-winProgressBar(title="SeqAlign Super Pro 9000",
                         label="Building reference index...",
                         min=0,max=100,width=300,initial=30)
    
    setwd(exp_dir)
    buildindex(basename="refseq",ref_seq) 
    
    #Update parms
    .GlobalEnv$parms$RefSeq<-ref_seq
    setwd(exp_dir)
    write.table(parms,file="Parms.txt",quote=FALSE,sep=",",eol="\n")
    
    close(prog)
    tk_messageBox(message="Index completed!")
    
  } else{
    
  }
  
  print(ref_seq)
}

display_reads<-function(){
  disp_reads<-gsub(paste(exp_dir,"/",sep=""),"",sel_reads[1:3])
  if(length(which(is.na(disp_reads)==TRUE))>0){
    disp_reads[is.na(disp_reads)]<-"NA" 
  }
  if(length(sel_reads)>3){
    disp_reads[4]<-"Additional reads hidden"
  } else{
    disp_reads[4]<-"NA"
  }
  tkconfigure(sel_reads_dis,text=paste(disp_reads,collapse="\n"))
}

select_reads<-function(){
  .GlobalEnv$sel_reads<-tk_choose.files(default=paste(exp_dir,"/READS",sep=""),
                                        filters=matrix(c("TXT",".txt","FASTA",".fa","All files","*"),3,2,byrow=TRUE))
  display_reads()
  
  #Update parms
  .GlobalEnv$parms$SelReads<-paste(sel_reads,collapse=" ")
  setwd(exp_dir)
  write.table(parms,file="Parms.txt",quote=FALSE,sep=",",eol="\n")
  
  print(sel_reads)
}

gui<-tktoplevel()
tkwm.geometry(gui,"300x330")
tkwm.title(gui,"SeqAlign Super Pro 2000")
frm<-tkframe(gui)
tkgrid(frm,padx=30)

ttl<-tklabel(frm,text="SeqAlign Super Pro 9000",font=tkfont.create(weight="bold",size=15))
tkgrid(ttl,row=1,column=1,columnspan=2)
open_exp<-tkbutton(frm,text="Open Experiment",command=open_exp_fun)
tkgrid(open_exp,row=2,column=1,columnspan=2,pady=5)
exp_display<-tklabel(frm,text="")
tkgrid(exp_display,row=3,column=1,columnspan=2,pady=5)

sel_ref_seq<-tkbutton(frm,text="Select reference sequence (FASTA)",command=select_ref,state="disabled")
tkgrid(sel_ref_seq,row=4,column=1,pady=5,columnspan=2)
ref_display<-tklabel(frm,text="")
tkgrid(ref_display,row=5,column=1,pady=6,columnspan=2)
sel_reads_var<-tkbutton(frm,text="Select reads to be aligned (FASTA)",command=select_reads,state="disabled")
tkgrid(sel_reads_var,row=6,column=1,pady=5,columnspan=2)
sel_reads_dis<-tklabel(frm,text="\n\n\n")
tkgrid(sel_reads_dis,row=7,column=1,pady=5,columnspan=2)
align_reads<-tkbutton(frm,text="Align reads",state="disabled",command=align_fun)
tkgrid(align_reads,row=8,column=1,pady=5,columnspan=2)

tkwait.window(gui)
