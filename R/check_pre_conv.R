check_pre_conv<-function(pre_values,PhiC){
  Phi0_visit <- pre_values$st.I
  num_Phi0 <- dim(PhiC)[1]/10
  hist(Phi0_visit,breaks = num_Phi0,main = 'Histogram of visited auxiliary phi')
}
