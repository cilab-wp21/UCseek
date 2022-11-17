UCseek_model_value <- function(weight,cnv_model_value,methy_model_score,cut-off)
{
  
  UCseek_model_value <- weight*cnv_model_value+(1-weight)* methy_model_score
  if(UCseek_model_value > cut-off){
    diagnostic_result <- "tumor"
  }else{
    diagnostic_result <- "normal"
  }
}
