#### compres probability matrix} 
# Create the PMatrix class
#
# This is used to represent the probabilities using the DISCOVER
# framework.
PMatrix <- setClass(
  # Set the name for the class
  "PMatrix",
  
  # Define the slots
  slots = c(
    rowExps = "numeric",
    colExps = "numeric"
  )
)

# TODO: decide which are the data stored on the matrix (logs or plain) and change the
# accesors accordingly. I think that plain is better: exponentials are avoid and
# we can use %*% instead of outer "+".

# TODO: implement row and colnames in the constructor and accessing to the matrix.

## Select rows
setMethod("[", signature(x = "PMatrix", i = "numeric", j = "missing"),
          function (x, i, j, ..., drop) { ## select rows
            na <- nargs()
            if(na == 3) {
              elog5 <- matrix(x@rowExps[i], ncol = 1) %*% matrix(x@colExps, nrow = 1)
              return(elog5 /(1+elog5))}
            else ## should not happen
              print(na)
              stop("PMatrix-internal error in <TsparseM>[i,,d]; please report")
          })

## Select columns
setMethod("[", signature(x = "PMatrix", i = "missing", j = "numeric"),
          function (x, i, j, ..., drop) { ## select columns
            na <- nargs()
            if(na == 3) {
              elog5 <- matrix(x@rowExps, ncol = 1) %*% matrix(x@colExps[j], nrow = 1)
              return(elog5 /(1+elog5))}
            else ## should not happen
              stop("PMatrix-internal error in <TsparseM>[i,,d]; please report")
          })
## Select columns
setMethod("[", signature(x = "PMatrix", i = "numeric", j = "numeric"),
          function (x, i, j, ..., drop) { ## select columns
            na <- nargs()
            if(na == 3) {
              elog5 <- matrix(x@rowExps[i], ncol = 1) %*% matrix(x@colExps[j], nrow = 1)
              return(elog5 /(1+elog5))}
            else ## should not happen
              stop("Matrix-internal error in <TsparseM>[i,,d]; please report")
          })

## A[ ij ]  where ij is (i,j) 2-column matrix :
setMethod("[", signature(x = "PMatrix",i = "matrix", j = "missing"),
            function (x, i, j, ..., drop){
              elog5 <- matrix(x@rowExps[i[1,]], ncol = 1) %*% matrix(x@colExps[i[,1]], nrow = 1)
              return(elog5 /(1+elog5))
            })

setMethod("as.matrix", signature(x = "PMatrix"),
          function (x, i, j, ..., drop){
            elog5 <- matrix(x@rowExps, ncol = 1) %*% matrix(x@colExps, nrow = 1)
            return(elog5 /(1+elog5))
          })


setMethod("dim", signature(x = "PMatrix"),
          function (x) { ## select columns
            return(c(length(x@rowExps),length(x@colExps)))
          })

# TODO: show only part of the matrix to avoid memory problems.
setMethod("print", signature(x = "PMatrix"),
          function (x, ...) { ## select columns
            print(x[1:length(x@rowExps),],...)
          })

# TODO: show only part of the matrix to avoid memory problems.
setMethod("show", "PMatrix",
          function (object) { ## select columns
            print(object[1:length(object@rowExps),])
          })
