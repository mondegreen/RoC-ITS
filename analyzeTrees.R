library(ape)
library(cluster)

args <- commandArgs(TRUE)

files = list.files(path=args[1],pattern="best.dnd$",full.names=TRUE,recursive=TRUE)
reps = 10

for (file in files) {
  tree<-read.tree(file=file)
  distance<-cophenetic(tree)
  mdist<-max(distance)
  cnts<-2:min(c(5,length(tree$tip.label)-1))
  if (length(tree$tip.label) > 2) {
    for (iter in cnts) {
      total = 0
      for (test in 1:reps) {
        clusters<-pam(distance,k=iter)
        successes = -1
        expect = -1
        last = -1
        pos = 0
        for (seq in sort(names(clusters$cluster))) {
          c = as.integer(clusters$cluster[seq])
          v = paste(last,c,sep="")
          if (expect == -1) {
            expect = c
          }
          print(paste(iter,seq,test,pos,c,v,expect,successes,mdist,sep=" "))
          if (c == expect) {
            successes = successes + 1
          } else {
            expect = c
          }
          expect = (c %% iter) + 1
          last = c
          pos = pos + 1
        }
        prob = 1 / iter
        bin = 1-pbinom(successes,length(clusters$cluster),prob)
        total = total + bin
      }
      bin = total / reps
    }
  }
}
