#' @title iCTC-identification of CTC
#'
#' @description this package stats::predicts whether -- is a CTC or
#' not using various preprocessing techniques and machine learning
#'
#' @param cell_samples, can contain any number of peripheral blood cells,
#' with atleast desired genes either in row or column
#'
#' @param cases, can contain any subset of (1,2,3,4,5,6,7,8,9) numbers,
#' meaning of these numbers given below:
#' 1-> Harmony_NB, 2-> Harmony_RF, 3-> Harmony_GBM
#'
#' 4-> PCA_NB, 5-> PCA_RF, 6-> PCA_GBM
#'
#' 7-> Original_NB, 8-> Original_RF, 9-> Original_GBM
#'
#' Harmony-> projects cells into a shared embedding in which cells group by
#' cell type rather than dataset-specific conditions.
#' PCA-> Principal component analysis (PCA) is a technique used to emphasize
#' variation and bring out strong patterns in a dataset.
#' Original means none of (Harmony, PCA).
#'
#' NB-> Naive Bayes(ML Model)
#' RF-> Random Forest(ML Model)
#' GBM-> Gradient Boosting Machine(ML Model)
#' @examples
#' library(harmony)
#' cell_samples=iCTC::raw_test_data$Clearcell_Polaris_sample_test
#' results<-iCTC(cell_samples=cell_samples, cases = c(4,5,6))
#' print(results)
#' @return results, will retrun table of samples and predicted values
#' corresponding cases which have given
#' row conatins cases and column with sample names
#' @export iCTC
iCTC <- function(cell_samples,cases=c(4,5,6))
{
    cases=as.integer(cases)
    if(any(cases>9) || any(cases<1))
    {
      print("please give any input in range 1 to 9")
      return()
    }
    gc()
    cases=sort(cases)
    cell_samples_name=0
    if(is.matrix(cell_samples))
    {
        cell_samples=cell_samples
    }else
    {
        cell_samples_name=names(cell_samples)
        cell_samples=matrix(cell_samples)
        rownames(cell_samples)=cell_samples_name
    }
    cell_samples=genes_filter(cell_samples)
    if(!is.matrix(cell_samples))
        return(cell_samples)
    print(cases)
    cell_samples_name=colnames(cell_samples)
    colnames(cell_samples)=paste(seq_len(length(colnames(cell_samples))),
        colnames(cell_samples),seq_len(length(colnames(cell_samples))),sep="_")
    n_models=3
    n_methods=3
    count<-list("harmony_count"=c(0),"pca_count"=c(0),"original_count"=c(0))
    for(i in seq_len(length(cases)))
    {
        for(j in seq_len(n_methods))
        {
            if(cases[[i]]<=(n_models*j) && cases[[i]]>=(n_models*(j-1)+1))
            {
                count[[j]] = count[[j]] + 1
            }
        }
    }
    methods=c("Harmony","PCA","Original")
    models=c("NB","RF","GBM")
    method_functions=list(harmony,pca,original)
    model_functions=list(NB,RF,GBM)
    for(j in seq_len(n_methods))
    {
        print(paste(methods[j]," count ", count[[j]], sep=""))
    }
    count<-list("harmony_count"=c(0),"pca_count"=c(0),"original_count"=c(0))
    all_result <- list()
    temp_method=0
    for(i in seq_len(length(cases)))
    {
        for(j in seq_len(n_methods))
        {
            if(cases[[i]]>=(n_models*(j-1)+1) &&
                cases[[i]]<=(n_models*j) && count[[j]]==0)
            {
                temp_method=j
                count[[j]]=1
                print(paste(methods[[j]], " correction running...", sep=""))
                preprocessed_data <- method_functions[[j]](cell_samples)
                print(paste(methods[[j]], " correction Done", sep=""))
             }
        }
        for(k in seq_len(n_models))
        {
            if(k==n_models)
                k=0
            if(cases[[i]] %% n_models==k)
            {
                gc()
                if(k==0)
                    k=n_models
                print(paste(models[[k]], " running...", sep=""))
                all_result <- append(all_result,
                    list(model_functions[[k]](temp_method,preprocessed_data)))
                print(paste(models[[k]], " Done", sep=""))
            }
        }
    }
    models_methods=c()
    for(i in cases)
    {
         models_methods=c(models_methods,
             paste(methods[as.integer(((i-1)/n_models)+1)],
                 models[as.integer(((i-1)%%n_models)+1)],sep="_"))
    }
    result_label=matrix(0,length(cases),length(cell_samples[1,]))
    rownames(result_label)=models_methods
    colnames(result_label)=cell_samples_name
    for(col in seq_len(length(cell_samples[1,])))
    {
        for(row in seq_len(length(cases)))
        {
            if(all_result[[row]]$te_labels[col]==0)
            {
                result_label[row,col]="CTC"
             }
             if(all_result[[row]]$te_labels[col]==1)
             {
              result_label[row,col]="Blood"
             }
        }
    }
    results=list("predicted_Labels"=result_label)
    return(results)
}

genes_filter <- function(gene_file)
{
    gc()
    common_genes = rownames(iCTC::Original_data_1$log_normalized_data_train_1)
    if(length(stats::na.omit(match(toupper(common_genes),
        toupper(rownames(gene_file))))>=0)==length(common_genes))
    {
        index=1
        print("All desired genes in row names found")
    }
    else if(length(stats::na.omit(match(toupper(common_genes),
        toupper(colnames(gene_file))))>=0)==length(common_genes))
    {
        index=2
        print("All desired genes in col names found")
    }else
    {
        index=3
        count_1=length(stats::na.omit(match(toupper(common_genes),
            toupper(rownames(gene_file))))>=0)
        count_2=length(stats::na.omit(match(toupper(common_genes),
            toupper(colnames(gene_file))))>=0)
        if(count_1>count_2 && count_2==0)
        {
            print(paste("Only ",count_1,
                " desired genes found in row names of your data",sep=""))
            print("please provide data with atleast desired genes
                either in col or row names")
            return(index)
        }
        if(count_2>count_1 && count_1==0)
        {
            print(paste("Only ",count_2,
                " desired genes found in col names of your data",sep=""))
            print("please provide data with atleast desired genes
                either in col or row names")
            return(index)
        }
        if(count_2>0 && count_1>0)
        {
            print(paste("Only ",count_1," and ",count_2,
                " desired genes found in  row and col names of your data",
                    sep=""))
            print("please provide data with atleast desired genes
                either in col or row names")
            return(index)
        }
        if(count_2==0 && count_1==0)
        {
            print("no single desired gene found in your data")
            print("please provide data with atleast desired genes
                either in col or row names")
            return(index)
         }
    }

    if(index==1)
    {
        index=apply(gene_file,2,
            function(x) sum(x>0)>=(length(gene_file[,1])*.06))
        if((length(gene_file[1,])-length(index))>0)
        {
            print(paste("only ",(length(gene_file[1,])-length(index)),
                " samples found with atlest 10 percent expressed genes
                    in column names of your data"))
        }else
        {
            print("All samples found with atlest 10 percent expressed genes
                in column names of your data")
        }
        return(gene_file[match(toupper(common_genes),
            toupper(rownames(gene_file))),index])
    }
    if(index==2)
    {
        index=apply(gene_file,1,
            function(x) sum(x>0)>=(length(gene_file[,1])*.1))
        if((length(gene_file[,1])-length(index))>0)
        {
            print(paste("only ",(length(gene_file[,1])-length(index)),
                " samples found with atlest 10 percent expressed genes
                    in column names of your data"))
        }else
        {
            print("All samples found with atlest 10 percent expressed genes
                in row names of your data ")
        }
        return(t(gene_file[index,match(toupper(common_genes),
            toupper(colnames(gene_file)))]))
    }
}
normalization <- function(data,method="log")
{
    gc()
    if (method == "log")
    {
        log_counts <- log(data + 1)
        print("Transformation Done")
        return(log_counts)
     }
     else if (method == "median")
     {
         normalized_matrix = as.matrix(data)
         total_count=colSums(data)
         med_total=stats::median(total_count)
         for(i in seq_len(ncol(normalized_matrix)))
         {
             normalized_matrix[,i] = data[,i] * (med_total/total_count[i])
         }
         print("Normalization Done")
         return(normalized_matrix)
    }else
    {
        return("error")
    }
}

pca <- function(user_data)
{
    gc()
    user_data=normalization(user_data,method = "median")
    user_data=normalization(user_data,method = "log")
    user_data=stats::predict(iCTC::pca_data_100,t(user_data))
    return(user_data)
}
harmony <- function(user_data)
{
    user_data=normalization(user_data,method = "median")
    user_data=normalization(user_data,method = "log")
    user_data=stats::predict(iCTC::pca_data_100,t(user_data))
    temp=rbind(iCTC::PCA_data_2$log_normalized_data_train_2,user_data)
    temp_1=matrix(0,ncol=4,nrow=dim(user_data)[1])
    temp_1[,1]=rownames(user_data)
    temp_1[,2]=rep("other_1",dim(user_data)[1])
    temp_1[,3]=rep("temp_1",dim(user_data)[1])
    temp_1[,4]=rep("sample_1",dim(user_data)[1])
    colnames(temp_1)=c("Type","Study","Sample_id","label")
    temp_2=matrix(0, ncol=4,
        nrow=dim(iCTC::PCA_data_2$log_normalized_data_train_2)[1])
    temp_2[,1]=rownames(iCTC::PCA_data_2$log_normalized_data_train_2)
    temp_2[,2]=rep("other_2",
        dim(iCTC::PCA_data_2$log_normalized_data_train_2)[1])
    temp_2[,3]=rep("temp_2",
        dim(iCTC::PCA_data_2$log_normalized_data_train_2)[1])
    temp_2[,4]=rep("sample_2",
        dim(iCTC::PCA_data_2$log_normalized_data_train_2)[1])
    colnames(temp_1)=c("Type","Study","Sample_id","label")
    temp_2=rbind(temp_2,temp_1)
    temp_2=data.frame(temp_2)
    harmonized_pcs <- harmony::HarmonyMatrix( data_mat = temp,
    meta_data = temp_2, vars_use  = "Study", do_pca    = FALSE)
    hyrdo_test_3=harmonized_pcs[which(unlist(temp_2$Study) %in% "other_1"),]

    return(hyrdo_test_3)
}

original <- function(user_data)
{
    gc()
    user_data=normalization(user_data,method = "median")
    user_data=normalization(user_data,method = "log")
    return(t(user_data))
}
NB <- function(data_case,te)
{
    gc()
    if(data_case==1)
    {
    test_prediction_1 =stats::predict(iCTC::nb_3,te)
    }
    if(data_case==2)
    {
    test_prediction_1 =stats::predict(iCTC::nb_2,te)
    }
    if(data_case==3)
    {
    data_2=data.frame(t(iCTC::Original_data_1$log_normalized_data_train_1),
                      "labels"=factor(c(rep(0,538),rep(1,1323))))
    test_prediction_1 =stats::predict(caret::train(data_2[,-c(dim(data_2)[2])],
        data_2[,dim(data_2)[2]],'nb',trControl=caret::trainControl(method='cv',
            number=5, returnResamp = "all")),te)
    }
    return_list <- list("te_labels"=as.numeric(as.matrix(test_prediction_1)))
    return(return_list)
}
RF <- function(data_case,te)
{
    gc()
    if(data_case==1)
    {
        test_prediction_1 =stats::predict(iCTC::rf_3,te)
    }
    if(data_case==2)
    {
        test_prediction_1 =stats::predict(iCTC::rf_2,te)
    }
    if(data_case==3)
    {
        test_prediction_1 =stats::predict(iCTC::rf_1,te)
    }
    return_list <- list("te_labels"=as.numeric(as.matrix(test_prediction_1)))
    return(return_list)
}

GBM <- function(data_case,te)
{
    gc()
    if(data_case==1)
    {
        test_prediction_1 =stats::predict(iCTC::gbm_3,te)
    }
    if(data_case==2)
    {
        test_prediction_1 =stats::predict(iCTC::gbm_2,te)
    }
    if(data_case==3)
    {
        test_prediction_1 =stats::predict(iCTC::gbm_1,te)
    }

    return_list <- list("te_labels"=as.numeric(as.matrix(test_prediction_1)))
    return(return_list)
}
