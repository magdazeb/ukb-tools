# 2020-01-08 code minor edit, add dir address
suppressMessages(library(limma))
suppressMessages(library(eulerr))
venn_analysis = function(grouplist,title='',dir='') {
    g_names = names(grouplist)
    unionlist = Reduce(union,grouplist)
    
    # Collapse the list to one dataFrame list
    unionPr=NULL; g_titles=NULL
    for(i in 1:length(grouplist)) {
        unionPr = cbind(unionPr,grouplist[[i]][match(unionlist,grouplist[[i]])])
        g_titles = c(g_titles,paste0(g_names[i],'\n',length(grouplist[[i]])))
    }
    rownames(unionPr) = unionlist
    
    # Generate binary table to match with ID list
    union = (unionPr != '')			# Transcform values to TRUE, if ID exists.
    union[is.na(union)] = FALSE		# Transform NA to FALSE value
    union = as.data.frame(union)	# Make 'union' to data.frame from
    colnames(union) = g_titles		# Names attach to venn diagram
    
    # Draw Venn Diagram
    venn_counts = vennCounts(union)
    venn_c = unlist(venn_counts[,4])
    f.name1 = paste0(dir,'/','venn_',g_names[1],'_',g_names[2],'.png')
    if(title=='') title = paste0('Venn analysis of ',nrow(union),' subjects')
    png(f.name1,width=10,height=10,units='in',res=100)
    if(length(grouplist)==3) {
        # ref: https://creately.com/blog/diagrams/venn-diagrams-vs-euler-diagrams/
        venn_fit = euler(c(
            'A'    =as.numeric(venn_c[5]), # How to automate this?!
            'B'    =as.numeric(venn_c[3]),
            'C'    =as.numeric(venn_c[2]),
            'A&B'  =as.numeric(venn_c[7]),
            'A&C'  =as.numeric(venn_c[6]),
            'B&C'  =as.numeric(venn_c[4]),
            'A&B&C'=as.numeric(venn_c[8])
                     )) # eulerr
        print(venn_fit)
        cat(paste0('> Euler fit is done.'))
        p=plot(
            venn_fit, quantities=T,
            labels=colnames(venn_counts)[1:3],
            edges=list(col=c('red','green','blue'),lwd=3),
            fills=list(fill=c('red','green','blue'),alpha=0.2),
            main=''
            ) # eulerr
        print(p)
    } else vennDiagram(union, main=title, circle.col=rainbow(length(g_titles))) # limma
    #replicate(sample(10,1), dev.new())
    graphics.off() # killing all devices
    cat(paste0('\nFigure draw: ',f.name1,'\n'))
    
    # Write files
    colnames(union) = g_names
    union.df = cbind(union,ids=rownames(union))
    venn.li  = list(union=union.df,vennCounts=venn_counts)
    f.name2  = paste0(dir,'/','venn_',g_names[1],'_',g_names[2],'.tsv')
    write.table(venn.li[[1]],f.name2,row.names=F,col.names=T,quote=F,sep='\t')
    cat(paste0('File write: ',f.name2,'\n'))
    return(venn.li)
}