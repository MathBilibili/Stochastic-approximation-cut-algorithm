CutModel<-function(px = px, py = py, prox = prox, rprox = rprox, proy = proy, rproy = rproy, Z = Z, Y = Y, d_x = d_x, d_y = d_y, cpp_yes = FALSE, cpp_package = ''){
  if(cpp_yes){
    rm(px)
    rm(py)
    library(cpp_package,character.only=TRUE)
    return(list(px=px,py=py, prox = prox, rprox = rprox, proy = proy, rproy = rproy, Z = Z, Y = Y, d_x = d_x, d_y = d_y))
  }else{
    return(list(px=px,py=py, prox = prox, rprox = rprox, proy = proy, rproy = rproy, Z = Z, Y = Y, d_x = d_x, d_y = d_y))
  }
}
