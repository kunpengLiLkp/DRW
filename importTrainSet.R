#start:2023-12-15
#finish:2024-03-02


library(org.Hs.eg.db)
library(dplyr)
#导入数据----------------------------------------------------------------------

#M <-  read.csv("C:/Users/admin/Desktop/实验记录/MtNormal.csv",header = TRUE,row.names = 1)
#colnames(M) <- rownames(M)
#M <- as.matrix(M)

M <- readRDS("D:/通路数据/MtNormal.rds")


read_directory <- "D:/下载/GEO数据/非小细胞肺癌/肺腺癌/"

read_file <- paste0(read_directory, "GSE19188_GPL570_1_12_0_XML.matrix.csv")

dataFilter <-  read.csv(read_file,header = TRUE)
dataFilter <-  na.omit(dataFilter)

# 提取t统计量和p值-------------------------------------------------------------------
#将case、control组计算t值-----------------------------------------Cancer\Normal-----

t_test_results <- data.frame()
t_test_p.value <- c()
for (i in 1:nrow(set)) {
  t_test_results <-  append(t_test_results,t.test(cancer[i,],normal[i,])$statistic)
  t_test_p.value <- append(t_test_p.value,t.test(cancer[i,],normal[i,])$p.value)
  # t.test(TU,NO,var.equal = FALSE)
  t_Statistic <- unlist(t_test_results)
}
t_result <- cbind(t_Statistic,t_test_p.value)

# t_result <- tTest(dataFilter)


start_time <- proc.time()
t_result <- tTest(dataFilter)
end_time <- proc.time()
run_time <- end_time - start_time
print(run_time)



#TotalData为geneID、t统计量、p值三列的矩阵
TotalData <- cbind(dataFilter$GeneId,t_result)
TotalData <-  na.omit(TotalData)
any_greater_than_0.05 <- any(TotalData[, 3] < 0.05)
#获取p值小于0.05的行

GeneID <- TotalData[which(TotalData[, 3] < 0.05),1]
DataF1 <- TotalData[which(TotalData[, 3] < 0.05),]

dataF2 <- DataF1
#计算W∞----------------
#将t统计量取绝对值

dataF2[,2] <- abs(as.numeric(dataF2[,2]))
dataF2[,1] <- paste("hsa:", dataF2[,1], sep = "")
hsacolName <- colnames(M)

W0 <- c()
for (i in 1:ncol(M)) {
  matching_indices <- which(hsacolName[i] == dataF2[,1])
  if(length(matching_indices) > 0){
    W0[i] <- which(hsacolName[i] == dataF2[,1]) %>% dataF2[.,2]
  }
}
# 将NA值替换为0---------------------------------------------
W0[is.na(W0)] <- 0
W0[(length(W0)+1):1676] <- 0
w0_test <- as.numeric(W0)

r <- 0.7
# 计算向量的模
vector_length <- sqrt(sum(w0_test^2))
# 归一化向量为单位向量
W0 <- w0_test / vector_length
# 计算随机游走后的W∞值
W1 <- wzhi(W0)

class(GeneID)

#计算通路活性-------------------------------------------------------------------
GeneID <- as.character(GeneID)
#找到基因对应的通路-------------------------------------------------------------
kegg_pathways <- select(org.Hs.eg.db, keys = GeneID, columns = c("ENTREZID", "PATH"), keytype = "ENTREZID")
#删除NA值的行-------------------------------------------------------------
Gene2path <- na.omit(kegg_pathways)
#将每个通路对应的基因ID整合为一个list-------------------------------------------------------------
G2PList <- split(Gene2path$ENTREZID,Gene2path$PATH)
G2PList <- na.omit(G2PList)
matched_list <- G2PList
#找到对应hsacolname对应列序号，相当于邻接矩阵的序号
TMP <- lapply(G2PList,function(x){
  # temp <- paste("hsa",x,sep = ":")
  match(x, hsacolName)
})
#过滤列表长度为0
#TMP <- Filter(function(x) length(x) > 0, TMP)
# 根据名称去除列表中的元素全为NA的向量
names_with_na <- names(TMP)[sapply(TMP, function(x) all(is.na(x)))]
TMP <- TMP[setdiff(names(TMP), names_with_na)]
matched_list <- matched_list[setdiff(names(matched_list), names_with_na)]
#根据基因ID获取游走后的W1值
RandomW <- lapply(TMP, function(x){ W1[x] })
# 将列表中所有的NA值替换为0
RandomW <- lapply(RandomW, function(x) {replace(x, is.na(x), 0)})
RandomWSum<- lapply(RandomW, function(x) {
  sqrt(sum(x^2))
})

# 获取对应基因id的t值
sgn <- lapply(matched_list, function(x){
  as.numeric(DataF1[match(x,DataF1[,1]),2])
})
# 对列表中的每个元素,即对对应基因的t值取sgn函数
sgn <- lapply(sgn, function(x) sign(x))

#zT(gi)是训练集中基因 gi在所有样本上行标准化的表达值向量

scaled_matrix <- apply(dataFilter[,-1], 1, function(x) scale(x)) %>% t()
colnames(scaled_matrix) <- colnames(dataFilter[,-1])
scaled_matrix <- cbind(dataFilter[,1],scaled_matrix)
colnames(scaled_matrix)[1] <- "GeneId"
dataF4 <- scaled_matrix[which(scaled_matrix[,1] %in% GeneID),]


#找到对应基因ID对应行的索引
index <- lapply(G2PList,function(x){
  match(x, dataF4[,1])
})

#计算各通路活性
for(i in 2:ncol(dataF4)){
  Zt <- lapply(index, function(x){ as.numeric(dataF4[x,i])})
  result <- mapply(function(x, y, z) x * y * z, RandomW, sgn, Zt)
  pathwayActive <- lapply(result, function(x){sum(x)}) %>% Map("/", ., RandomWSum)
  if(i == 2){
    PathwaysAc <-  data.frame(Name = names(pathwayActive), Value = unlist(pathwayActive))
  }else{
    PathwaysAc <-  cbind.data.frame(PathwaysAc,unlist(pathwayActive))
      }
}
pathwayProfile <- na.omit(PathwaysAc)
colnames(pathwayProfile) <- colnames(dataF4)


# 获取邻接矩阵-------------
library(igraph)
library(dplyr)
metaAdjacency <- get.adjacency(graph_hsaMetabolism,sparse = FALSE)

# 将邻接矩阵行归一化--------------------
M1 <- apply(metaAdjacency, 1, function(x){
  x/sum(x)
})

#L1范数
L1_norm <- function(x,y){
  return(sum(abs(x-y)))
}

# 计算随机游走后的W∞值
wzhi <- function(W0){
  zhongzhi <- 1e-10
  W1 <-  0 * W0
  temp <- 0 * W0
  W2 <- W0
  while(L1_norm(W2,W1) >= zhongzhi){
    W1 <- temp
    W2 <- (1-r)* M %*% W1 + r * W0
    temp <- W2
  }
  return(W1)
}


#将case、control组计算t值------------------------------------T\N--------
tTest <- function(set){
  cancer <- set[,grep("T", colnames(set))]
  normal <- set[,grep("N", colnames(set))]
  t_test_results <- data.frame()
  t_test_p.value <- c()
  for (i in 1:nrow(set)) {
    t_test_results <-  append(t_test_results,t.test(cancer[i,],normal[i,])$statistic)
    t_test_p.value <- append(t_test_p.value,t.test(cancer[i,],normal[i,])$p.value)
    # t.test(TU,NO,var.equal = FALSE)
    t_Statistic <- unlist(t_test_results)
  }
  return(cbind(t_Statistic,t_test_p.value))
}
#将case、control组计算t值------------------------------------Tumor\Normal--------
tTest <- function(set){
  cancer <- set[,grep("Tumor", colnames(set))]
  normal <- set[,grep("Normal", colnames(set))]
  t_test_results <- data.frame()
  t_test_p.value <- c()
  for (i in 1:nrow(set)) {
    t_test_results <-  append(t_test_results,t.test(cancer[i,],normal[i,])$statistic)
    t_test_p.value <- append(t_test_p.value,t.test(cancer[i,],normal[i,])$p.value)
    # t.test(TU,NO,var.equal = FALSE)
    t_Statistic <- unlist(t_test_results)
  }
  return(cbind(t_Statistic,t_test_p.value))
}

#将case、control组计算t值------------------------------------Tumor\healthy--------
tTest <- function(set){
  cancer <- set[,grep("tumor", colnames(set))]
  normal <- set[,grep("healthy", colnames(set))]
  t_test_results <- data.frame()
  t_test_p.value <- c()
  for (i in 1:nrow(set)) {
    t_test_results <-  append(t_test_results,t.test(cancer[i,],normal[i,])$statistic)
    t_test_p.value <- append(t_test_p.value,t.test(cancer[i,],normal[i,])$p.value)
    # t.test(TU,NO,var.equal = FALSE)
    t_Statistic <- unlist(t_test_results)
  }
  return(cbind(t_Statistic,t_test_p.value))
}




#测试无效-------------------------------------------
#Map(function(x){
#temp <- paste("hsa",x,sep = ":")
#which(hsacolName %in% temp)
#},G2PList,seq_along(G2PList))

#result <- Map(function(x, y) {
#which(hsacolName %in% x)  # 找到与G2PList中的值相等的行的索引
#}, TMP, seq_along(TMP))
#matched_list <- Map(c, G2PList, TMP)
#-------------------------------------------
#train Set Zt 归一化向量为 W0---------错误
##Zt <- lapply(TMP, function(x){ W0[x]}) %>% lapply(., function(x) {replace(x, is.na(x), 0)})
#lung_tumor_columns <- grep("Lung.Cancer", colnames(dataFilter1))
#normal_lung_columns <- grep("Lung.Normal", colnames(dataFilter1))
#lung_tumor <- dataFilter[lung_tumor_columns]
#normal_lung <- dataFilter[normal_lung_columns]
# 提取t统计量和p值-------------------------------------------------------------------
#t_test_results <- data.frame()
#t_test_p.value <- c()
#rbind(t_test_results,1)
#for (i in 1:nrow(dataFilter)) {
#  t_test_results <-  append(t_test_results,t.test(lung_tumor[i,],normal_lung[i,])$statistic)
#  t_test_p.value <- append(t_test_p.value,t.test(lung_tumor[i,],normal_lung[i,])$p.value)
#  # t.test(TU,NO,var.equal = FALSE)
#}
#t_Statistic <- unlist(t_test_results)




#使用fdr方法对原始p值进行校正





