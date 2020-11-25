#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


ss2mrf = function(ss) {
  j_pair = structure(
    c(
      -1.30412,
      3.96785,-1.30474,-1.31304,-0.00854209,
      3.91698,-1.30412,-1.29931,-1.2937,-0.00849558,-1.2937,-1.31304,
      -1.28455,
      3.86993,-0.00847145,-1.29931,-1.30474,
      3.91198,
      -1.28455,-0.00840413,-0.00849558,-0.00854209,-0.00840413,
      -0.00847145,
      0
    )
  )
  
  h_single = c(1, 1, 1, 1,-4)
  
  
  df.pair = RNAmrf:::ss2pairs(ss)
  line_h = sprintf("V[%d] %s", df.pair$id - 1, paste0(h_single, collapse =
                                                        " "))
  
  df.pair.tmp = df.pair[df.pair$pair > df.pair$id, ]
  
  line_v = sprintf("W[%d][%d] %s",
                   df.pair.tmp$id - 1,
                   df.pair.tmp$pair - 1,
                   paste0(j_pair, collapse = " "))
  
  mrflines = c(line_h, line_v)
  
  return(mrflines)
}

#### main

# ss.file=args[1]
fsto=args[1]

dfref=RNAmrf:::read_dfref(fsto)

mrflines=ss2mrf(dfref$ss)

# f=readLines(ss.file)
writeLines(mrflines)
