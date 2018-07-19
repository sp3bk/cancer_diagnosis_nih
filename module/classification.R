draw_class<-function(train_path,test_path,classtype){
  library(stringr)
  
  library(h2o)
  library(ROCR)
  
  train_data<-read.table(train_path,header=F,sep=",")
  test_data<-read.table(test_path,header=F,sep=",")
  
  if(classtype=="Naive bayes"){
    library(e1071)
    #View(train_data)
    #cat(ncol(train_data))
    classifier<-naiveBayes(train_data[,-ncol(train_data)], train_data[,ncol(train_data)])
    test_pred <- predict(classifier, test_data[,-ncol(train_data)],type="raw")
    #View(test_pred )
    brain_score <- test_pred[, c("brain")]
    brain_actual_class <- test_data[,ncol(train_data)] == 'brain'
    #View(brain_score)
    #View(brain_actual_class)
    kidney_score <- test_pred[, c("kidney")]
    kidney_actual_class <- test_data[,ncol(train_data)] == 'kidney'
    ovar_score <- test_pred[, c("ovarian")]
    ovar_actual_class <- test_data[,ncol(train_data)] == 'ovarian'
    
    brain_pred <- prediction(brain_score, brain_actual_class)
    brain_perf <- performance(brain_pred , measure = "tpr", x.measure = "fpr")
    
    #tiff(file="./naive.tiff")
    plot(brain_perf,col="red",main="ROC Curve (Naive bayes)")
    
    kidney_pred <- prediction(kidney_score, kidney_actual_class)
    kidney_perf <- performance(kidney_pred , measure = "tpr", x.measure = "fpr")
    lines(kidney_perf@x.values[[1]], kidney_perf@y.values[[1]], col = "blue")
    
    ovar_pred <- prediction(ovar_score, ovar_actual_class)
    ovar_perf <- performance(ovar_pred , measure = "tpr", x.measure = "fpr")
    lines(ovar_perf@x.values[[1]], ovar_perf@y.values[[1]], col = "green")
    
    legend('bottomright', c("brain","kidney","ovarian"),lty=1, col=c('red', 'blue','green'), bty='n', cex=.75)
    #dev.off()
    
  }
  
  if(classtype=="Deep learning(CNN)"){
    library(h2o)
    library(stringr)
    localH2O<-h2o.init(ip = "localhost", startH2O = TRUE,max_mem_size = '7g')
    train.hex = h2o.importFile(path = train_path,destination_frame = "train.hex")
    test.hex = h2o.importFile(path = test_path,destination_frame = "test.hex")
    cancer.dl <- h2o.deeplearning(x = 1:(ncol(train_data)-1),y=ncol(train_data), training_frame= train.hex)
    deep_pred<-as.data.frame(h2o.predict(cancer.dl, newdata=test.hex))
    
    brain_score <- deep_pred[, c("brain")]
    brain_actual_class <- test_data[,ncol(train_data)] == 'brain'
    #View(brain_score)
    #View(brain_actual_class)
    kidney_score <- deep_pred[, c("kidney")]
    kidney_actual_class <- test_data[,ncol(train_data)] == 'kidney'
    ovar_score <- deep_pred[, c("ovarian")]
    ovar_actual_class <- test_data[,ncol(train_data)] == 'ovarian'
    
    brain_pred <- prediction(brain_score, brain_actual_class)
    brain_perf <- performance(brain_pred , measure = "tpr", x.measure = "fpr")
    
    #tiff(file="./naive.tiff")
    plot(brain_perf,col="red",main="ROC Curve (Naive bayes)")
    
    kidney_pred <- prediction(kidney_score, kidney_actual_class)
    kidney_perf <- performance(kidney_pred , measure = "tpr", x.measure = "fpr")
    lines(kidney_perf@x.values[[1]], kidney_perf@y.values[[1]], col = "blue")
    
    ovar_pred <- prediction(ovar_score, ovar_actual_class)
    ovar_perf <- performance(ovar_pred , measure = "tpr", x.measure = "fpr")
    lines(ovar_perf@x.values[[1]], ovar_perf@y.values[[1]], col = "green")
    
    legend('bottomright', c("brain","kidney","ovarian"),lty=1, col=c('red', 'blue','green'), bty='n', cex=.75)
    
    
  
  }
  
}

#---functions
draw.roc<-function(cnv_score,cnv_actual_class,add_score,add_actual_class,clf,indx){
  
  cnv_pred <- prediction(cnv_score, cnv_actual_class)
  cnv_perf <- performance(cnv_pred , measure = "tpr", x.measure = "fpr")
  title<-paste("ROC Curve",clf,indx,sep="_")
  tiff(file="./temp.tiff")
  plot(cnv_perf,col="red",main=title)  
  
  
  add_pred <- prediction(add_score, add_actual_class)
  add_perf <- performance(add_pred , measure = "tpr", x.measure = "fpr")
  lines(add_perf@x.values[[1]], add_perf@y.values[[1]], col = "blue")
  
  legend('bottomright', c("CNV","CNV+DNA methylation"),lty=1, col=c('red', 'blue'), bty='n', cex=.75)
  dev.off()
}