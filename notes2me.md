########NEXT#######
```R packages: reshape2, ggplot2, org.Hs.eg.db,

**analysis**:  #For views
  average_scross.R :  
    input : robustness.rds, simu_subset.rds, simu_noise.rds
    outer: ggplot plots(geom_violin, geom_point)
  compareperf.R:
    input: anno/compiled/all.csv
    outer:  ggplot plots
  gpt4version.R: # gpt4aug3 vs gpt4mar23
    input: anno/compiled/all.csv
    outer: ggplot plots  
  gptperf.R: #
    input: anno/compiled/all.csv, numcell/numcell.rds, anno/subtype/list.csv, analysis/method_color.rds
    outer: ggplot plots
  method_color.R:
    input: NA 
    outer: analysis/method_color.rds
  robustness.R:  #simu data
    input: simu/gpt4aug3/res/mixanno.csv, simu/gpt4aug3/res/null.csv
    outer: plot/robustness.rds , ggplot plots
  runtime.R: #runtime compare
    input: analysis/method_color.rds, runtime/SingleR.rds, runtime/sctype.rds, runtime/gpt4aug3.rds, runtime/gpt3.5aug3.rds
    outer: ggplot plots
  simu_noise.R:
    input:  simu/gpt4aug3/compiled/noise.csv, anno/compiled/all.csv
    outer: plot/simu_noise.rds, ggplot plots
  simu_subset.R:
    input: simu/gpt4aug3/compiled/subset.csv, anno/compiled/all.csv
    outer: plot/simu_subset.rds, ggplot plots
  stromal.R:
    input: anno/compiled/all.csv, analysis/stromalexpr.rds
    outer: ggplot plots
  top.R:
    input: analysis/method_color.rds, anno/compiled/gpt4topgenenumber.csv, 
    outer: ggplot plots

**numcell**:  # mult cancer 
  make.R: 
    input: celltype/cancer/bcl/data/proc/ct.rds, celltype/cancer/coloncancer/data/proc/ct.rds, celltype/hca/data/proc/*/ct.rds, celltype/hcl/data/proc/anno.rds,celltype/mca/data/proc/anno.rds, celltype/tabulasapiens/data/meta/*
    outer: celltype/numcell/numcell.rds 

**runtime**:
  SingleR.R: #from log
  gpt3.5aug3.R #when run gpt4 , system.time() function return time
  gpt4aug3.R #
  sctype.R #from log

**simu/gpt4aug3**:
  **compiled**:
    noise.R: 
      input: simu/gpt4aug3/res/noise.csv  #?? ,anno/cl/compiled.csv, anno/cl/relation.csv
      outer: simu/gpt4aug3/compiled/noise.csv
    noise.csv
    subset.R:
      input: simu/gpt4aug3/res/subset.csv,anno/cl/compiled.csv, anno/cl/relation.csv
      outer: simu/gpt4aug3/compiled/subset.csv
    subset.csv
    
  **res**:
    annonoise.R: #调用gptcelltype得到annotaion, d$annotation[id] <- gptcelltype(input=gl[id],openai_key='',model='gpt-4')
      input: simu/gpt4aug3/noise.csv
      outer: simu/gpt4aug3/noise.csv (add gptcelltype annotation)
    annsubset.R:
      input: simu/gpt4aug3/subset.csv
      outer: simu/gpt4aug3/subset.csv (add gptcelltype annotation)
    mixanno.R:
      input: anno/GPT4_mar23/tabulasapiens_nuc_lit/proc.csv
      outer: simu/mixanno.csv
    mixanno.csv
    
    noise.R:  # 随机掺入 n 个基因， #set.seed
      input: anno/source/tabulasapiens_nuc_lit/proc.csv
      am <- do.call(rbind,sapply(1:5,function(simuround) {  ##new rows = row(d)* 5 * 4 #扩大20倍 
	            do.call(rbind,sapply(c(0.25,0.5,0.75,1),function(per) {  ## noise: 0.25,0.5,0.75,1 
		              tmpd <- d
		              tmpd$gene <- sapply(d$gene,function(i) {
			            gl <- strsplit(i,', ')[[1]]  # genelist 
			            snum <- floor(length(gl)*per) # 根据noise 确定gene 个数：eg. 20 * 0.5 =10
			            rnum <- length(gl)-snum  #
			            paste0(c(gl,sample(gn,snum)),collapse = ',') #gn : org.Hs.eg.db的所有基因名字（77510个基因） ； 在gl 原来基础上随机掺入snum个 基因
		              })
		        data.frame(simuround=simuround,noise=per,celltype=tmpd$celltype,gene=tmpd$gene)
	          },simplify = F))
            },simplify = F))
	write.csv(am,'Dropbox/research/gptcelltype/simu/gpt4aug3/noise.csv',row.names = F)
      	outer: simu/gpt4aug3/noise.csv #扩大20倍的模拟数据
    noise.csv
    
    null.R：#人为混入celltype未知的数据
      input: anno/tabulasapiens_nuc_lit/proc.csv
          am <- do.call(rbind,sapply(1:10,function(simuround) { # 混合数据：10*20 行 
                  m <- t(sapply(1:10,function(id) {
                  c('random',paste0(sample(gn,10),collapse = ',')) #混入 celltype="random", gene= gn中的随机10个基因组成
                }))
                colnames(m) <- c('celltype','gene')
                m <- rbind(m,d[sample(1:nrow(d),10),2:3]) #从已知中提取10行，将已知+人为的20行合并
                m <- cbind(simuround,m) 
            },simplify = F))
      outer: simu/null.csv
    null.csv
    
    subset.R：#和 noise.R 区别: (1)noise c(0.25,0.5,0.75,1) 调整为 c(0.25,0.5,0.75), (2) gene ：  paste0(c(gl,sample(gn,snum)),collapse = ',') 调整为 paste0(c(sample(gl,snum)),collapse = ',')， 只包含随机的基因
      input: anno/source/tabulasapiens_nuc_lit/proc.csv
      outer: simu/gpt4aug3/subset.csv
    subset.csv
