formatFeatures = function(Model = Model, train_data = train_data, test_data = test_data, test = Test){
  features = unique(Model$name[2:nrow(Model)])
  for(f in features){
    tmp_f = dplyr::select(train_data, y, all_of(f))
    colnames(tmp_f) = c("y", "value")
    
 #   ggplot(tmp_f, aes(x = as.factor(y), y = value))+
    #  geom_boxplot()
    
    fg = dplyr::filter(tmp_f, y == 1)
    bg = dplyr::filter(tmp_f, y == 0)
    Model$Feature_increase[Model$name == f] = mean(fg$value) - mean(bg$value)
    Model$Feature_increase_p[Model$name == f] = wilcox.test(x = fg$value, y = bg$value)$p.value
    
    if(test == T){
      tmp_f = dplyr::select(test_data, y, all_of(f))
      colnames(tmp_f) = c("y", "value")
      fg = dplyr::filter(tmp_f, y == 1)
      bg = dplyr::filter(tmp_f, y == 0)
      Model$Feature_increase_test[Model$name == f] = mean(fg$value) - mean(bg$value)
      Model$Feature_increase_p_test[Model$name == f] = wilcox.test(x = fg$value, y = bg$value)$p.value
    }else{
      Model$Feature_increase_test[Model$name == f] = NA
      Model$Feature_increase_p_test[Model$name == f] = NA
    }
    
  }
  return(Model)
}
