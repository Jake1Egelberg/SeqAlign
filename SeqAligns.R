

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
                "BiocManager",
                "ggplot2",
                "dplyr")
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
setWinProgressBar(open_prog,value=65,label="Loading ggplot2...")
library(ggplot2)
setWinProgressBar(open_prog,value=75,label="Loading dplyr...")
library(dplyr)

#Install BiocManager packages
setWinProgressBar(open_prog,value=85,label="Installing BiocManager packages...")
bioc.packages<-c("Rsubread")
for(i in 1:length(bioc.packages)){
  tryCatch(find.package(bioc.packages[i]),
           error=function(e) {
             setWinProgressBar(open_prog,value=85,label=paste("Installing ",bioc.packages[i],"...",sep=""))
             BiocManager::install(bioc.packages[i])
           })
}

setWinProgressBar(open_prog,value=95,label="Loading Rsubread...")
library(Rsubread)

close(open_prog)

#**********************************************

#--------------------------------START ANALYSIS

#**********************************************

.GlobalEnv$dir<-this.dir()
.GlobalEnv$skip_subread<-tclVar(1)
.GlobalEnv$align_method_var<-tclVar("Midpoint")
.GlobalEnv$n_cross_vals<-10
.GlobalEnv$n_samples<-10

align_fun<-function(){
  
  .GlobalEnv$skip_subread_val<-tclvalue(skip_subread)
  
  if(is.na(ref_seq)==FALSE&&length(sel_reads)>=1&&is.na(sel_reads)==FALSE){
  
    prog_inc<-100/(length(sel_reads)+1)
    prog_val<-0
    prog<-winProgressBar(title="SeqAlign Super Pro 9000",
                         label=paste("Aligning...", sep=""),
                         min=0,max=100,width=300,initial=prog_val)
    
    if(skip_subread_val=="0"){
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
      write.csv(alignment_output,"ALIGNMENT_OUTPUT.csv",row.names=FALSE)
      
      #Get unaligned/aligned reads
      .GlobalEnv$unaligned_reads<-sel_reads[which(alignment_output$Unmapped_reads>0)]
      .GlobalEnv$aligned_reads<-sel_reads[-which(sel_reads%in%unaligned_reads)]
      
    } else{
      .GlobalEnv$aligned_reads<-sel_reads
    }

    #Generate plots
    setWinProgressBar(prog,value=90,label="Generating alignment metadata...")
    get_alignment()
    
    close(prog)
    tk_messageBox(message="Alignment complete and saved to experiment directory!")
    
    edit(align_meta)
    
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
  tkconfigure(skip_sub,state="normal")
  
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

get_alignment<-function(){
  
  #Get alignment method
  .GlobalEnv$align_method_val<-tclvalue(align_method_var)
  
  align_meta<-data.frame()
  for(i in 1:length(aligned_reads)){
    
    dir_ref_seq<-toupper(paste(read.delim(ref_seq)[,1],collapse=""))
    read_cur<-paste(toupper(read.delim(aligned_reads[i])[,1]),collapse="")
    
    ref_length<-str_count(dir_ref_seq)
    read_length<-str_count(read_cur)
    
    #Align midpoint of read to ref seq
    if(align_method_val=="Midpoint"){
      
        #Get shorter/longer sequence
        if(read_length>ref_length){
          .GlobalEnv$chosen_seq<-read_cur
          .GlobalEnv$chosen_seq_name<-"Read"
          
          .GlobalEnv$others<-dir_ref_seq
          .GlobalEnv$others_name<-"Reference"
        } else{
          .GlobalEnv$chosen_seq<-dir_ref_seq
          .GlobalEnv$chosen_seq_name<-"Reference"
          
          .GlobalEnv$others<-read_cur
          .GlobalEnv$others_name<-"Read"
        }
        
        #Get 10 nucleotides from midpoint of shorter sequence
        midpoint<-str_count(others)/2
        nucs<-str_sub(others,start=midpoint,end=midpoint+9)
        
        #Align to other seq
        chosen_loc<-str_locate(chosen_seq,nucs)[[1]]

        #If alignment is successul
        if(is.na(chosen_loc)==FALSE){
          align_succcess<-"Y"
          
          shift<-abs(chosen_loc-midpoint)
          
          #Create alignment matrix from longer sequence
          alignment_matrix<-matrix(nrow=3,ncol=str_count(chosen_seq))
          colnames(alignment_matrix)<-c(1:ncol(alignment_matrix))
          rownames(alignment_matrix)<-c(chosen_seq_name,others_name,"Alignment")
          
          alignment_matrix[1,]<-unlist(strsplit(chosen_seq,""))
          alignment_matrix[2,(shift+1):(str_count(others)+shift)]<-unlist(strsplit(others,""))
          
          #Get alignments
          aligned_status<-lapply(1:ncol(alignment_matrix),function(p){
            sub_m<-alignment_matrix[,p]
            if(is.na(sub_m[[2]])==TRUE){
              aligned<-"NA"
            } else if(is.na(sub_m[[2]])==FALSE&&sub_m[[2]]!=sub_m[[1]]){
              aligned<-"N"
            } else if(is.na(sub_m[[2]])==FALSE&&sub_m[[2]]==sub_m[[1]]){
              aligned<-"Y"
            }
          })
          aligned_status_cur<-unlist(aligned_status)
          alignment_matrix[3,]<-aligned_status_cur
          
          #Generate aligment_df
          .GlobalEnv$alignment_df<-data.frame(Nucleotide=1:ncol(alignment_matrix),
                                   Align_Success=alignment_matrix[3,],
                                   Seq=nucs)
          
          
          #Calculate alignment success relative to shorter sequence
          perc_aligned<-100*(length(which(alignment_matrix[3,]=="Y"))/str_count(others))
          
          #END IF ALIGNMENT SUCCESSFUL
        } else{
          align_succcess<-"N"
          perc_aligned<-0
          shift<-"NA"
        }
        
      #Save output
      results<-data.frame(Read=gsub(paste(exp_dir,"/",sep=""),"",aligned_reads[i]),
                          Alignment_Success=align_succcess,
                          Perc_Align=perc_aligned)
      align_meta<-rbind(align_meta,results)
      
      if(align_succcess=="Y"){
        
        #Generate graph
        best_df<-alignment_df
        names(best_df)[2]<-"Aligned"
        
        num_aligned<-length(which(best_df$Aligned=="Y"))
        perc_aligned_long<-round(100*num_aligned/str_count(chosen_seq),1)
        perc_aligned_short<-round(100*num_aligned/str_count(others),1)
        ggplot(best_df,aes(x=as.numeric(Nucleotide),y=1))+
          geom_point(size=3,col=ifelse(best_df$Aligned=="Y","darkgreen","red"))+
          geom_text(aes(x=as.numeric(min(Nucleotide)),y=1),label=paste("% Aligned (",chosen_seq_name, ") = ",perc_aligned_long,"%\n",num_aligned,"/",str_count(chosen_seq)," nucleotides","\n% Aligned (",others_name,")= ",perc_aligned_short,"%\n",num_aligned,"/",str_count(others)," nucleotides\nMETHOD=MIDPOINT",sep=""),hjust=0,vjust=-1)+
          scale_x_continuous(n.breaks=10)+
          xlab(paste("Nucleotide Position of ", chosen_seq_name,sep=""))+
          ylab(paste("Alignment of ", chosen_seq_name,sep=""))+
          theme_bw()+
          theme(axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor=element_blank())+
          ggtitle(paste(gsub(paste(exp_dir,"/",sep=""),"",aligned_reads[i])," Alignment to Reference",sep=""))
        setwd(exp_dir)
        ggsave(paste(gsub(paste(exp_dir,"/",sep=""),"",aligned_reads[i]),"_PLOT.png",sep=""),width=10,height=5)
        
      } else{
        print(paste("No successful alignments for ",gsub(paste(exp_dir,"/",sep=""),"",aligned_reads[i]),sep=""))
      }
      
    } else if(align_method_val=="Random"){
      
      #Define parameters
      all_alignments<-data.frame()
      alignment_df_all<-data.frame()
      for(q in 1:n_cross_vals){
        set.seed(q)
        
        #Get shorter/longer sequence
        if(read_length>ref_length){
          seq_samps<-sample(1:(read_length-11),n_samples)
          .GlobalEnv$chosen_seq<-read_cur
          .GlobalEnv$chosen_seq_name<-"Read"
          
          .GlobalEnv$others<-dir_ref_seq
          .GlobalEnv$others_name<-"Reference"
        } else{
          seq_samps<-sample(1:(ref_length-11),n_samples)
          .GlobalEnv$chosen_seq<-dir_ref_seq
          .GlobalEnv$chosen_seq_name<-"Reference"
          
          .GlobalEnv$others<-read_cur
          .GlobalEnv$others_name<-"Read"
        }
        
        #Get 10 nucleotides for each seq samp in list
        seq_samp_df<-data.frame()
        for(m in seq_samps){
          #Get random 10 nucs in chosen seq
          nucs<-str_sub(chosen_seq,start=m,end=m+9)
          
          #Align to other seq
          others_loc<-str_locate(others,nucs)[[1]]
          
          #If alignment successful
          if(is.na(others_loc)==FALSE){
            align_succcess<-"Y"
            
            #Get shift to align shorter sequence
            shift<-m-others_loc
            
            #Create alignment matrix from longer sequence
            alignment_matrix<-matrix(nrow=3,ncol=str_count(chosen_seq))
            colnames(alignment_matrix)<-c(1:ncol(alignment_matrix))
            rownames(alignment_matrix)<-c(chosen_seq_name,others_name,"Alignment")
            
            alignment_matrix[1,]<-unlist(strsplit(chosen_seq,""))
            alignment_matrix[2,(shift+1):(str_count(others)+shift)]<-unlist(strsplit(others,""))
            
            #Get alignments
            aligned_status<-lapply(1:ncol(alignment_matrix),function(p){
              sub_m<-alignment_matrix[,p]
              if(is.na(sub_m[[2]])==TRUE){
                aligned<-"NA"
              } else if(is.na(sub_m[[2]])==FALSE&&sub_m[[2]]!=sub_m[[1]]){
                aligned<-"N"
              } else if(is.na(sub_m[[2]])==FALSE&&sub_m[[2]]==sub_m[[1]]){
                aligned<-"Y"
              }
            })
            aligned_status_cur<-unlist(aligned_status)
            alignment_matrix[3,]<-aligned_status_cur
            
            #Generate aligment_df
            alignment_df<-data.frame(Nucleotide=1:ncol(alignment_matrix),
                                     Aligned=alignment_matrix[3,],
                                     Seed=q,
                                     Seq=nucs)
            alignment_df_all<-rbind(alignment_df_all,alignment_df)
            
            #Calculate alignment success relative to shorter sequence
            perc_aligned<-100*(length(which(alignment_matrix[3,]=="Y"))/str_count(others))
            
            #End if alignment successful
          } else{
            align_succcess<-"N"
            perc_aligned<-0
            shift<-"NA"
          }
          
          #Format as df
          seq_samp_df<-data.frame(Read=gsub(paste(exp_dir,"/",sep=""),"",aligned_reads[i]),
                                  Nucleotide_Selected=m,
                                  Sample_Sequence=nucs,
                                  Align_Success=align_succcess,
                                  Perc_Aligned=perc_aligned,
                                  Seed=q,
                                  Shift=shift)
          
          #End cycle thru random samples
          all_alignments<-rbind(all_alignments,seq_samp_df)
        }
  
      #End cross validations
      }
      .GlobalEnv$all_alignments<-all_alignments
      .GlobalEnv$alignment_df_all<-alignment_df_all
      
      #------------------------------------
      
      #Remove duplicates
      if(length(which(duplicated(all_alignments$Sample_Sequence)==TRUE))>0){
        all_alignments<-all_alignments[-which(duplicated(all_alignments$Sample_Sequence)==TRUE),]
      }
      
      #Rank alignments
      .GlobalEnv$ranked_al<-all_alignments[order(all_alignments$Perc_Aligned,decreasing=TRUE),]
      
      #Obtain best alignment
      best_al<-ranked_al[1,]
      best_al$Unique_Alignments_Attempted<-nrow(ranked_al)
      best_al$Perc_Successful_Alignments<-round(100*length(which(ranked_al$Align_Success=="Y"))/nrow(ranked_al),1)
      .GlobalEnv$best_al<-best_al
      
      align_meta<-rbind(align_meta,best_al)
      
      if(align_succcess=="Y"){
        
        #Generate graph
        best_df<-alignment_df_all[which(alignment_df_all$Seq==best_al$Sample_Sequence),]
        num_aligned<-length(which(best_df$Aligned=="Y"))
        perc_aligned_long<-round(100*num_aligned/str_count(chosen_seq),1)
        perc_aligned_short<-round(100*num_aligned/str_count(others),1)
        ggplot(best_df,aes(x=as.numeric(Nucleotide),y=1))+
          geom_point(size=3,col=ifelse(best_df$Aligned=="Y","darkgreen","red"))+
          geom_text(aes(x=as.numeric(min(Nucleotide)),y=1),label=paste("% Aligned (",chosen_seq_name, ") = ",perc_aligned_long,"%\n",num_aligned,"/",str_count(chosen_seq)," nucleotides","\n% Aligned (",others_name,")= ",perc_aligned_short,"%\n",num_aligned,"/",str_count(others)," nucleotides\nMETHOD=RANDOM",sep=""),hjust=0,vjust=-1)+
          scale_x_continuous(n.breaks=10)+
          xlab(paste("Nucleotide Position of ", chosen_seq_name,sep=""))+
          ylab(paste("Alignment of ", chosen_seq_name,sep=""))+
          theme_bw()+
          theme(axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor=element_blank())+
          ggtitle(paste(gsub(paste(exp_dir,"/",sep=""),"",aligned_reads[i])," Alignment to Reference",sep=""))
        setwd(exp_dir)
        ggsave(paste(gsub(paste(exp_dir,"/",sep=""),"",aligned_reads[i]),"_PLOT.png",sep=""),width=10,height=5)
        
      } else{
        print(paste("No successful alignments for ",gsub(paste(exp_dir,"/",sep=""),"",aligned_reads[i]),sep=""))
      }
      
      #End random option
    }
    #End cycle thru reads
  }
  
  rownames(align_meta)<-NULL
  .GlobalEnv$align_meta<-align_meta
  setwd(exp_dir)
  write.csv(align_meta,"ALIGNMENT_METADATA.csv",row.names = FALSE)
 
}

gui<-tktoplevel()
tkwm.geometry(gui,"300x380")
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
sel_align_method<-tklabel(frm,text="Alignment method:")
tkgrid(sel_align_method,row=8,column=1,pady=0)
align_method<-ttkcombobox(frm,values=c("Midpoint","Random"),width=12,textvariable=align_method_var)
tkgrid(align_method,row=8,column=2,pady=0)
align_reads<-tkbutton(frm,text="Align reads",state="disabled",command=align_fun)
tkgrid(align_reads,row=9,column=1,pady=5,columnspan=2)
skip_sub<-tkcheckbutton(frm,text="Skip RSubread alignment?",variable=skip_subread,state="disabled")
tkgrid(skip_sub,row=10,column=1,pady=5,columnspan=2)


tkwait.window(gui)
