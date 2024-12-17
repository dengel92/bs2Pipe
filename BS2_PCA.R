pca.plot <- function(pca.df, cent){
  df.clean <- pca.df %>% 
    column_to_rownames("ids") %>%
    select(-contains("_rank")) %>% 
    select(-contains("_top"))
  sepcent <- map(cent, function(x){
    tmp <- t(subset(df.clean, , select = grep(x, names(df.clean))))
  }) %>% set_names(cent)
  
  ## PCA Plot
  
  pc <- map(sepcent, prcomp)%>% set_names(cent)
  ggobj <- map2(pc, sepcent, ~ggplot2::autoplot(.x, data = .y, label = F))%>% set_names(cent)
  pca <- map2(ggobj, cent, function(X, y){
    tmp <- X + ggtitle(paste0("PCA: ",y)) + geom_point(size=1) + 
      geom_label_repel(aes(label = rownames(sepcent[[y]]), colour = rownames(sepcent[[y]])), size = 6)+
      # scale_colour_brewer(type="qual", palette = 2)+
      theme_bw()+
      theme(text = element_text(size=14),
            legend.title = element_text(size=15),
            legend.text=element_text(size=15),
            axis.title.x = element_text(size=20),
            axis.title.y = element_text(size=20))
  })
  
  ## t-sne Plot
  set.seed(42)
  tsne_model <- map(sepcent, function(x){
    tmp <- Rtsne(as.matrix(x), 
                 check_duplicates=FALSE, 
                 pca=TRUE, 
                 perplexity=1, 
                 theta=0.5, 
                 dims=3, 
                 set.seed=TRUE)
  }) %>% set_names(cent)
  
  ## getting the two dimension matrix
  ## df_tsne <- as.data.frame(tsne_model$Y)
  df_tsne <- map2(tsne_model, cent, function(x, y){
    data.frame(PC1 = x$Y[,1],
               PC2 = x$Y[,2],
               group = rownames(sepcent[[y]]))
  }) %>% set_names(cent)
  
  
  ggtsne <- map(df_tsne, ~ggplot(.x, aes(x=PC1, y=PC2, group=group), label= F)) %>% set_names(cent)
  
  tsne <- map2(ggtsne, cent, function(X, y){
    tmp <- X + ggtitle(paste0("T-sne Plot: ",y)) + geom_point(size=1) + 
      geom_label_repel(aes(label = rownames(sepcent[[y]]), colour = rownames(sepcent[[y]])), size = 6)+
      # scale_colour_brewer(type="qual", palette = 2)+
      theme_bw()+
      theme(text = element_text(size=14),
            legend.title = element_text(size=15),
            legend.text=element_text(size=15),
            axis.title.x = element_text(size=20),
            axis.title.y = element_text(size=20))
  })
  
  ## RPCA
  rpc <- map(sepcent, ~tryCatch(
    PcaGrid(as.matrix(.x), 3),
    error = function(e) {
      message("Error with dataset: ", e$message)
      NULL
    }
  )) %>% set_names(cent)
  
  rpc_df <- map2(rpc, cent, ~tryCatch(
    function(x, y){
    data.frame(X = x$sd,
               Y = x$od,
               group = rownames(sepcent[[y]]))
  },
  error = function(e) {
    message("Error with dataset: ", e$message)
    NULL
  })) %>% set_names(cent)
  ggrpca <- map2(rpc_df, rpc, function(X, Y) {
    tryCatch({
      ggplot(data = X, aes(x = X, y = Y, group = group)) +
        geom_hline(yintercept = Y$cutoff.od) +
        geom_vline(xintercept = Y$cutoff.sd) +
        xlab("Score Distance") +
        ylab("Orthogonal Distance")
    }, error = function(e) {
      message("Error in ggrpca creation for RPC: ", conditionMessage(e))
      NULL
    })
  }) %>% set_names(cent)
  
  rpca <- map2(ggrpca, cent, function(X, Y) {
    tryCatch({
      tmp <- X + ggtitle(paste0("RPCA: ", Y)) +
        geom_point(size = 1) +
        geom_label_repel(aes(label = group, colour = group), size = 6) +
        theme_bw() +
        theme(
          text = element_text(size = 14),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20)
        )
      tmp
    }, error = function(e) {
      message("Error in RPCA plot for centrality ", Y, ": ", conditionMessage(e))
      NULL
    })
  }) %>% set_names(cent)
  
  ## RPCAScore
  rpcscre <- map(rpc, function(r) {
    tryCatch({
      getPrcomp(r)
    }, error = function(e) {
      message("Error in computing RPCAScore: ", conditionMessage(e))
      NULL
    })
  }) %>% set_names(cent)
  
  ggrpcascre <- map2(rpcscre, sepcent, function(X, Y) {
    tryCatch({
      ggplot2::autoplot(X, data = Y, label = F)
    }, error = function(e) {
      message("Error in RPCAScore plotting: ", conditionMessage(e))
      NULL
    })
  }) %>% set_names(cent)
  
  rpcascre <- map2(ggrpcascre, cent, function(X, y) {
    tryCatch({
      tmp <- X + ggtitle(paste0("RPCA Score: ", y)) +
        geom_point(size = 1) +
        geom_label_repel(aes(label = rownames(sepcent[[y]]), colour = rownames(sepcent[[y]])), size = 6) +
        theme_bw() +
        theme(
          text = element_text(size = 14),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20)
        )
      tmp
    }, error = function(e) {
      message("Error in RPCA Score plot for centrality ", y, ": ", conditionMessage(e))
      NULL
    })
  }) %>% set_names(cent)
  
  
  
  res <- lst(PCA = pca, Tsne = tsne, RPCA = rpca, RPCAscore = rpcascre)
  
}