ComEnvir<-function(is_Unix = TRUE, core_num = 1, clusterExport = list('py','Y','Npop','dmvnorm')){
  if(is_Unix){
    cl<-makeCluster(core_num,type = 'FORK')
    registerDoParallel(cl)
    return(list(is_Unix=is_Unix,cl=cl))
  }else{
    #cl<-makeCluster(detectCores( ))
    cl<-makeCluster(core_num)
    registerDoParallel(cl)
    return(list(is_Unix=is_Unix, cl=cl, clusterExport=clusterExport))
  }
}
