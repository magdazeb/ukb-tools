# command.r
# 2019-12-27 First version in UKB-menopause analysis
# This function is developed for capturing arguments from command line.
command = function(
    commandLine = NULL,
    help_message = NULL
) {
    ## Capture the call for Help Messages.
    help = (sum(c('--help','-h') %in% commandLine) >= 1)
    prun = (sum(c('--ukb_prune','-pr') %in% commandLine) >=1)
    excl = (sum(c('--ukb_excld','-ex') %in% commandLine) >=1)
    glob = (sum(c('-b','-t','--verbose','-v') %in% commandLine) >=1)
    n = length(help_message)
    if(help & prun) { message(help_message[c(1:2, 3:4 ,n)]); quit()
    } else if(help & excl) { message(help_message[c(1:2, 5:6 ,n)]); quit()
    } else if(help & glob) { message(help_message[c(1:2, 3,5, n)]); quit() 
    } else message(help_message); quit()
    
    ## If help is not invoked, we can start processing.
    comm = paste(unlist(commandLine),collapse=' ')
    listoptions = unlist(strsplit(comm,'--|-'))[-1]
    args = sapply(listoptions,function(x) {
        arg = unlist(strsplit(x,' '))[-1]
        if(length(arg)==0) { return(TRUE)
        } else return(arg) })
    args.names = sapply(listoptions,function(x) {
        option = unlist(strsplit(x,' '))[1] })
    names(args) = unlist(args.names)
    return(args)
}