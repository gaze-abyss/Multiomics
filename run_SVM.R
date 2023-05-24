####################update 2023/5/24
##########SVM
setwd("/path/test/")           #set path

library(data.table)
library(e1071)
library(caret)
set.seed(12345)

train <- fread('1ibaq_rm0log_genesymbol_T_50_knn_scale.txt',header = T)
test <- fread('D2_T_scale.txt',header = T)

cli <- read.table('group.txt',header = T)
train <- as.data.frame(train)
test <- as.data.frame(test)
rownames(train) <- train$id
train <- train[,-1]
rownames(test) <- test$id
test <- test[,-1]
gene <- intersect(rownames(test),rownames(train))
train <- train[gene,]
test <- test[gene,]

cli$group <- factor(cli$group)
train <- t(train)
train <- cbind(train,group=cli$group)
train <- as.data.frame(train)
x <- train[,1:(ncol(train)-1)]
y <- factor(train[,ncol(train)])
test <- t(test)

tc <- tune.control(cross = 5)
model <- tune.svm(x = x,
                    y = y,
                    kernel = 'linear',
                    cost = c(0.001,0.01,0.1,1,5,10),
                    tunecontrol = tc)
n=model$best.parameters
model <- tune.svm(x = x,
y = y,
kernel = 'linear',
cost = n)
summary(model)
svm_model <- model$best.model
svm_train <- predict(svm_model,x)
result_train <- cbind(cli,svm_train)
head(result_train)
svmtest<- predict(svm_model,test)
output <- cbind(row.names(test),svmtest)
colnames(output)[2] <- "group"
output  # classification results