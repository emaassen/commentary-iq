### This is simulation code for the condition 'Idealworld'
### OSF:
### Email: esther@emaassen.com 
### This is a TryCatch function which also saves the warnings and errors
### Credit to user2161065: https://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function

myTryCatch <- function(expr) {
  
  warn <- err <- NA                                                         # if no errors occur, let warnings and errors be NA
  value <- withCallingHandlers(tryCatch(expr, error=function(e) {
      
      err <<- e                                                             # if an error occurs, save it in err
      
    }), warning=function(w) {
      
      warn <<- w                                                            # if a warning occurs, save it in warn
      invokeRestart("muffleWarning")
      
    })
  
  list(value=value, warning=warn, error=err)                                
}



