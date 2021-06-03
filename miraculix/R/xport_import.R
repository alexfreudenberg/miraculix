
# optionsDelete <- function(RFopt, register=NULL) {}

## BEIM AUFRUF: entweder mit COPY = FALSE oder nachfolgend mit !hasArg("COPY")
getRFoptions <- function(...) RFoptions(local_=TRUE, ..., no.class=TRUE)
setRFoptions <- function(...) {
  RFoptions(local_=TRUE, ..., no.class=TRUE)
}


internalRFoptions <- function(..., local_=TRUE, getoptions_=NULL,
                              saveoptions_=NULL, FORMER = FALSE, RETURN = TRUE){
  ## saveoptions_ darf nicht in ... auftauchen, da dann ...length() > 0
  ## um saveoptions_ nicht als Extrafall zu behalten, verschmelzung mit
  ## getoptions_:
   if (hasArg("saveoptions_"))
     getoptions_ <- unique(c(getoptions_, saveoptions_, "basic", "general"))

  if (FORMER)
    former <-
      if (length(getoptions_) == 0) RFoptions(local_=local_, no.class=TRUE)
      else RFoptions(local_=local_, getoptions_=getoptions_, no.class=TRUE)
  
  
  setoptions <- ...length() > 0
  if (TRUE || setoptions) { ## 28.1.21 setoptions ausgeblended,
    ## da 1. eh alle RFfoo mit internalRFoptions() o.ae. versehen sein muessen
    ## 2. es Aerger gibt, bei RFfoo(...) wenn weder Benutzer noch
    ## Programmiere eine Option uebergeben
    
    ## verhindert, das eine RFfoo ohne ... parameter aufgerufen wird,
    ## im Glauben, RFfoo setzt keine Optionen. Und dann wird hinterruecks
    ## die gesetzten optionen ueberschrieben.
    
    if (hasArg("COPY")) {
      L <- list(...)
      idx <- names(L) == "COPY"
      if (all(unlist(L[idx]))) {
##        Print("Copy option")
        .Call(C_copyoptions)
      }
      L <- L[!idx]
      if (length(L) > 0) do.call("RFoptions", c(list(local_=local_), L))
      setoptions <- FALSE
    } else {
##      Print("Copy option")
      .Call(C_copyoptions) ## from global to KT (local)
    }
  }
  
  ## sequence: local_, warnUnknown_, LIST / saveoptions_ / GETOPTION
   if (setoptions) {
    RFoptions(local_=local_, ...)
  } 

  if (RETURN) {
    new <- if (length(getoptions_) == 0) RFoptions(local_=local_, no.class=TRUE)
           else RFoptions(local_=local_, getoptions_=getoptions_, no.class=TRUE)
    if (FORMER) list(former, new) else new
  } else if (FORMER) former
  
}

