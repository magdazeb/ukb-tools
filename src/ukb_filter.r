# 2019-12-26 First version in UKB-menopause analysis
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
ukb_exclude = function(
    ukb_pheno = NULL,  # pruned UKB RDS file
    fid       = NULL,  # interesting field ids
    cutoff    = NULL,  # option: data type is continuous or integer,
                       #     and cutoff criteria such as (<9 or >20) years
    cat_no    = FALSE, # option: data type is categorical data and select 'No' answers
    verbose   = FALSE
) {
    paste0('Excluding/filtering data...') %>% cat
    if(length(cutoff)>0 & cat_no) {
        paste0('Error: You cannot use both options cutoff and cat_no.') %>% cat
        stop() }
    if(length(fid)==1) {
        ukb_fid  = ukb_pheno %>% select(eid,starts_with(fid))
        fid_name_n = colnames(ukb_fid) %>% length
        if(fid_name_n==1) {
            paste0('Error: "',fid,'" is not found.\n') %>% cat
            stop()
        } else {
            fid_name = colnames(ukb_fid)[2]
            paste0('\n  ',fid_name,'\n') %>% cat
        }
        
        eid_fid = na.omit(ukb_fid)$eid
        out = eid_fid
        if(verbose) {
            paste0('\n  na.omit = ',length(eid_fid),'\n') %>% cat
            summary(ukb_fid[,2]) %>% print
        }
        if(length(cutoff)>0) {
            which_fid = which(
                ukb_fid[,2]>=cutoff[1]&
                ukb_fid[,2]<=cutoff[2] )
            out = ukb_fid[which_fid,]$eid
            if(verbose) {
                cutoff_ = paste(cutoff,collapse=",")
                paste0(
                    '  filtered by (',cutoff_,') = ',
                    length(which_fid),'\n') %>% cat
                summary(ukb_fid[which_fid,2]) %>% print
            }
        }
    } else if(cat_no) {
        n = length(fid)
        eid_fid_li = list()
        for(i in 1:n) {
            ukb_fid = ukb_pheno %>% select(eid,starts_with(fid[i]))
            fid_name_n = colnames(ukb_fid) %>% length
            if(fid_name_n==1) {
                paste0('Error: "',fid,'" is not found.\n') %>% cat
                stop() }
            fid_name = colnames(ukb_fid)[2]
            paste0('\n[',i,'/',n,'] ',fid_name) %>% cat
            
            eid_fid = na.omit(ukb_fid)$eid
            if(i==1) {
                eid_fid_no = ukb_fid[which(ukb_fid[,2]==0),]$eid
                eid_fid_li[[i]] = ukb_fid[which(ukb_fid[,2]==1),]$eid
            } else eid_fid_li[[i]] = eid_fid
            if(verbose) {
                '\n' %>% cat
                if(i==1) {
                    table(ukb_fid[,2]) %>% print
                    paste0('  eid_fid_no = ',length(eid_fid_no),'\n') %>% cat
                } else summary(ukb_fid[,2]) %>% print
                paste0(
                    '  total = ',nrow(ukb_fid),
                    ', na.omit = ',length(eid_fid),'\n') %>% cat
            }
        }
        #length(eid_fid_li) %>% print
        eid_uni = Reduce(union,eid_fid_li)
        eid_set = setdiff(eid_fid_no,eid_uni)
        if(verbose) {
            paste0('\n- union(yes, start, ...) = ',length(eid_uni),'\n',
                '- setdiff(no, uion)      = ',length(eid_set),'\n' ) %>% cat
        }
        out = eid_set
    }
    return(out)
}

ukb_prune = function(
    dir        = NULL,  # directory for figures
    ukb_pheno  = NULL,  # download file from ukb
    fid        = NULL,  # interesting field id
    except_fid = NULL,  # exceptional fids. No filtering
    fid_info   = NULL,  # UKB field id information table
    fid_code   = NULL,  # UKB field code information table
    verbose    = FALSE
) {
    paste0('1. Prepare...') %>% cat
    except_fid = except_fid
    fid_which  = which(fid_info$FieldID==fid)
    fid_field  = fid_info$Field[fid_which]
    fid_unit   = fid_info$Units[fid_which]
    fid_vtype  = fid_info$ValueType[fid_which]
    paste0(' [',fid,'; ',fid_vtype,'] ',fid_field,'\n') %>% cat
    cat_type   = c('Categorical single','Categorical multiple','Date')
    con_type   = c('Integer','Continuous')
    if(!fid_vtype%in%c(cat_type,con_type)) {
        paste0('Error: "',fid_vtype,'" is not considered yet.\n') %>% cat
        stop() }
    
    fid_code_ = fid_code %>% filter(FieldID==fid)
    fid_code_ = na.omit(fid_code_) # 19.12.20 debug
    if(verbose) {
        if(nrow(fid_code_)==0) {'  There is no fid code information.\n' %>% cat
        } else print(fid_code_[,5:6])
    }
    #table(fid_info[,1:2]) %>% addmargins %>% print # show contingency table of field categories
    
    paste0('2. Calculate overlap fields...') %>% cat
    ukb_fid = ukb_pheno %>% select(eid,starts_with(fid))
    paste0(' (cols= ',dim(ukb_fid)[1],', rows= ',dim(ukb_fid)[2],')') %>% cat
    #head(ukb_fid) %>% print # for debugging
    ukb_fid[ukb_fid==''] = NA # 19.12.23 debug
    if(nrow(fid_code_)>0) { code = fid_code_$Code[which(fid_code_$Code<0)]
    } else code = NULL
    if(length(code)>0 & ncol(ukb_fid)>1 & !is.null(code)) {
        n = length(code)
        for(j in 1:n) {
            if(verbose) { paste0('\n  [',j,'/',n,'] Converting: ',fid_code_$Code[j],' >> ',
                fid_code_$Meaning[j]) %>% cat }
            ukb_fid[ukb_fid==code[j]] = NA }
        #ukb_fid[ukb_fid==-1] = NA # -1: Do not know
        #ukb_fid[ukb_fid==-3] = NA # -3: Prefer not to answer
        paste0(' done\n') %>% cat
    } else {
        paste0(', stop\n') %>% cat
        if(verbose) paste0('  There is no convertable info/data.\n') %>% cat
    }
    ukb_eid   = ukb_fid[,1] %>% unlist %>% as.character
    ukb_fid_c = ukb_fid[,-1] # only contents
    if(ncol(ukb_fid)>1) {
        if(ncol(ukb_fid)==2) ukb_fid_c = as.data.frame(ukb_fid_c)
        fid_over = apply(ukb_fid_c,1,function(row){ length(row)-sum(is.na(row)) })
    } else fid_over = NA
    if(verbose & ncol(ukb_fid)>1) table(fid_over) %>% addmargins %>% print
    #return(ukb_fid) # debugging
    
    paste0('3. Calculate answer difference...') %>% cat
    if(ncol(ukb_fid)==1) {
        paste0(' No data, pass\n') %>% cat
    } else if(fid_vtype %in% cat_type) {
        paste0(' Categorical data, pass\n') %>% cat
    } else {
        ukb_diff = apply(ukb_fid_c,1,function(row) {
            row_ = na.omit(row)
            n = length(row_)
            if(n>1) { return(row_[1]-row_[n]) # ?! max(row_)-min(row_)
            } else return(0)
        }) 
        paste0(' done\n') %>% cat
    }
    #length(ukb_diff) %>% print
    
    paste0('4-a. Find cutoff decile thresholds...') %>% cat
    if(ncol(ukb_fid)==1) {
        paste0(' No data, pass\n') %>% cat
    } else if(length(which(fid_over>1))==0) {
        paste0(' No duplicated value, pass\n') %>% cat
    } else if(fid_vtype %in% cat_type) {
        paste0(' Categorical data, pass\n') %>% cat
    } else if(fid %in% except_fid) {
        paste0('fid:',fid,', pass\n') %>% cat
        ukb_diff_over = ukb_diff[which(fid_over>1)]
    } else {
        ukb_eid_over    = ukb_eid[which(fid_over>1)]
        ukb_diff_over   = ukb_diff[which(fid_over>1)]
        ukb_diff_decile = ntile(ukb_diff_over,10)
        cutoff = c(
            max(ukb_diff_over[which(ukb_diff_decile==1)]),
            min(ukb_diff_over[which(ukb_diff_decile==10)]) )
        if(fid_vtype %in% c('Integer')) {
            cutoff = c(cutoff[1]+1,cutoff[2]-1) }
        diff_over_decile = data.frame(
            eid = ukb_eid_over,
            diff_decile = ukb_diff_decile )
        paste0(' done\n') %>% cat
        if(verbose) paste0(
            '  length= ',length(ukb_diff_over),
            ', left= ',cutoff[1],', right= ',cutoff[2],'\n') %>% cat
    }
    
    paste0('4-b. Draw answer diff histogram plot...')%>%cat
    if(ncol(ukb_fid)==1) {
        paste0(' No data, pass\n') %>% cat
    } else if(length(which(fid_over>1))==0) {
        paste0(' No duplicated value, pass\n') %>% cat
    } else if(fid_vtype %in% cat_type) {
        paste0(' pass\n') %>% cat
    } else {
        if(fid %in% except_fid) {
            ukb_diff_over_ = rep(TRUE,length(ukb_diff_over))
        } else {
            left = cutoff[1]; right = cutoff[2]
            if(fid_vtype %in% c('Integer')) {
                ukb_diff_over_ = ifelse(
                ukb_diff_over %in% c(left:right),TRUE,FALSE) # Set Cutoff criteria
            } else {
                ukb_diff_over_ = rep(TRUE,length(ukb_diff_over))
                ukb_diff_over_[ukb_diff_over<=left]  = FALSE
                ukb_diff_over_[ukb_diff_over>=right] = FALSE
        }   }
        
        ukb_diff_over_df = data.frame(
            answer_diff = ukb_diff_over,
            cutoff      = factor(ukb_diff_over_,levels=c(TRUE,FALSE)) )
        #summary(ukb_diff_over_df) %>% print
        p_title = paste0(
            '[',fid,'; ',fid_vtype,'] ',fid_field,'\n',
            'Total n = ',which(fid_over>0) %>% length,
            '; Multiple answer n = ',nrow(ukb_diff_over_df))
        x_label = paste0('answer diff (',fid_unit,')')
        #y_label = paste0('count (log10 scale)')
        p_colors = c('#56B4E9','Red') #'#999999'
        if(fid_vtype %in% c('Integer')) { # 19.12.23 debug
            p1=ggplot(ukb_diff_over_df,aes(answer_diff,fill=cutoff))+theme_bw()+
            geom_histogram(binwidth=1)+
            xlab(x_label)+
            ggtitle(p_title)+
            #scale_y_continuous(trans='log10')+
            scale_fill_manual(values=p_colors)
        } else {
            p1=ggplot(ukb_diff_over_df,aes(answer_diff,fill=cutoff))+theme_bw()+
            geom_histogram(bins=100)+
            xlab(x_label)+
            ggtitle(p_title)+
            #scale_y_continuous(trans='log10')+
            scale_fill_manual(values=p_colors) }
        if(verbose) print(p1)
        f_name1 = paste0(dir,'/fid',fid,'_diff.png')
        ggsave(f_name1,width=8,height=9)
        paste0(' done\n') %>% cat
        if(verbose) paste0('  Fig generated: ',f_name1,'\n') %>% cat
    }
    
    paste0('5. Prepare to return pruned values...') %>% cat
    ukb_fid_c2 = ukb_fid_c
    if(ncol(ukb_fid)==1) {
        which_id = NULL
        paste0(' No data, pass\n') %>% cat
    } else if(length(which(fid_over>1))==0) {
        which_id = which(fid_over>0)
    } else if(fid_vtype %in% c('Categorical single','Categorical multiple')) {
        which_id = which(fid_over>0)
        code_ = fid_code_$Code[which(fid_code_$Code>=0)]; #print(code_)
        meaning_ = fid_code_$Meaning[which(fid_code_$Code>=0)]; #print(meaning_)
        if(length(code_)>=0) {
            for(k in 1:length(code_)) {
                #if(verbose) paste0(code_[k],':',meaning_[k]) %>% print # for debug
                ukb_fid_c2[ukb_fid_c==code_[k]] = meaning_[k] } }
    } else if(fid_vtype %in% c('Date') | fid %in% except_fid ) {
        which_id = which(fid_over>0)
    } else if(fid_vtype %in% c('Integer')) { # 19.12.23 debug
        which_id = intersect(
            which(ukb_diff %in% c(left:right)),
            which(fid_over>0) )
    } else {
        which_id = Reduce(intersect,list(
            which(ukb_diff>left),
            which(ukb_diff<right),
            which(fid_over>0) )) }
    if(verbose) paste0('  Pruned n= ',length(which_id)) %>% cat
    
    # Get the latest answer for return
    if(ncol(ukb_fid)==1) {
        out = ukb_fid
        if(verbose) paste0(', and Return is only eid.') %>% cat
    } else {
        if(ncol(ukb_fid)==2) { # 19.12.26 debug
            ukb_fid_c2 = ukb_fid_c2 %>% as.data.frame
            ukb_fid_c  = ukb_fid_c  %>% as.data.frame
        }
        ukb_fid_c_  = ukb_fid_c2[which_id,] %>% as.data.frame
        ukb_fid_c_2 = ukb_fid_c[which_id,]  %>% as.data.frame
        out_c = apply(ukb_fid_c_,1,function(row) {
            row_ = na.omit(row)
            n = length(row_)
            return(row_[n]) })
        out_c2 = apply(ukb_fid_c_2,1,function(row) {
            row_ = na.omit(row)
            n = length(row_)
            return(row_[n]) })
        if(fid_vtype %in% c('Date')) out_c = as.POSIXct(out_c)
        out = data.frame( # categories are converted.
            eid = ukb_eid[which_id],
            fid = out_c )
        out2 = data.frame( # categories are not converted.
            eid = ukb_eid[which_id],
            fid = out_c2)
    }
    if(!fid %in% except_fid & 
       fid_vtype %in% c('Integer','Continuous') & 
       length(which(fid_over>1))>1) {
        out_decile = merge(out,diff_over_decile,by='eid',all.x=T)
        out_decile$diff_decile = as.factor(out_decile$diff_decile)
    }
    paste0(' done\n') %>% cat
    
    paste0('6. Draw a result plot...') %>% cat
    # Draw result plot
    if(!ncol(ukb_fid)==1) {
        f_name2 = paste0(dir,'/fid',fid,'.png')
        p_title = paste0(
            '[',fid,'; ',fid_vtype,'] ',fid_field,'\n',
            'Total n = ',which(fid_over>0) %>% length,
            '; Pruned n = ',nrow(out)) }
    
    if(ncol(ukb_fid)==1) { paste0(' >> Cannot draw a plot.\n') %>% cat
    } else if(length(which(fid_over>1))==0) {
        # Draw integer-histogram
        p2 = ggplot(out,aes(fid,..count..))+theme_bw()+
            geom_histogram(binwidth=1)+
            xlab(fid_unit)+
            ggtitle(p_title)
        ggsave(f_name2,width=7,height=6)
    } else if(fid_vtype %in% c('Categorical single','Categorical multiple')) {
        p_title = paste0(
            fid_vtype,'\n',
            '[',fid,'] ',fid_field,'\n',
            'Total n = ',which(fid_over>0) %>% length,'\n',
            'Pruned n = ',nrow(out))
        out_df = table(out_c) %>% data.frame
        colnames(out_df) = c('Answer','Hits')
        #if(verbose) print(out_df)
        # Draw bar plot
        p2 = ggplot(out_df,aes(x=Answer,y=Hits))+theme_bw()+
            geom_bar(stat='identity',fill='steelblue')+
            theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
            xlab(fid_unit)+
            ggtitle(p_title)
        ggsave(f_name2,width=4,height=5)
    } else if(fid_vtype %in% c('Date')) {
        # Draw date-barplot
        p2 = ggplot(out,aes(fid))+theme_bw()+
            geom_histogram(bins=100)+
            scale_x_datetime(
                breaks = date_breaks('24 month'), # scales library
                labels = date_format('%Y-%b') )+  # scales library
            xlab(fid_unit)+
            ggtitle(p_title)
        ggsave(f_name2,width=7,height=6)
    } else if(fid_vtype %in% c('Integer')) {
        # Draw integer-histogram
        if(fid %in% except_fid | length(which(fid_over>1))==0) {
            p2 = ggplot(out,aes(fid,..count..))+theme_bw()+
                geom_histogram(binwidth=1)+
                xlab(fid_unit)+
                ggtitle(p_title)
            ggsave(f_name2,width=7,height=6)
        } else {
            p2 = ggplot(out_decile,aes(fid,..count..,fill=diff_decile))+theme_bw()+
                geom_histogram(binwidth=1)+
                xlab(fid_unit)+
                ggtitle(p_title)
            ggsave(f_name2,width=7,height=6) }
    } else {
        # Draw continuous histogram
        if(fid %in% except_fid | length(which(fid_over>1))==0) {
            p2 = ggplot(out,aes(fid))+theme_bw()+
                geom_histogram(bins=100)+
                xlab(fid_unit)+
                ggtitle(p_title)
            ggsave(f_name2,width=7,height=6)
        } else {
            p2 = ggplot(out_decile,aes(fid,fill=diff_decile))+theme_bw()+
                geom_histogram(bins=100)+
                xlab(fid_unit)+
                ggtitle(p_title)
            ggsave(f_name2,width=7,height=6) }
    }
    if(exists('p2')) paste0(' done\n') %>% cat
    if(verbose&exists('p2')) print(p2)
    if(ncol(ukb_fid)>1) colnames(out2) = c('eid',paste0('fid',fid,', ',fid_field))
    if(verbose&exists('p2')) paste0('  Fig generated: ',f_name2,'\n') %>% cat
    if(verbose) paste0('  Return table nrow= ',nrow(out),'\n') %>% cat
    if(ncol(ukb_fid)==1) { ukb_fid
    } else return(out2)
}