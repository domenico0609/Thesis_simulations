NMVar = 3
FOutputs = 12

FST=data.frame(
  FST = c("mito Fst p1&p2", "mito Fst p1&p3", "mito Fst p1&p4", "mito Fst p2&p3", "mito Fst p2&p4", "mito Fst p3&p4", "nuc Fst p1&p2", "nuc Fst p1&p3", "nuc Fst p1&p4", "nuc Fst p2&p3", "nuc Fst p2&p4", "nuc Fst p3&p4"),
  NM0 = double(FOutputs),
  NM2 = double(FOutputs),
  NM10 = double(FOutputs))

for (n in 1:NMVar) {
  Miss = double(FOutputs)
  Nan_val = double(FOutputs)
  
  for (run in 1:runs) {
    l = nrow(large_results[[n]][[run]]$output_data)
    if (l < FOutputs & l > 0) FOutputs = nrow(large_results[[n]][[run]]$output_data)
    
    for (output in 1:FOutputs) {
      if (length(large_results[[n]][[run]]$output_data) != 0) 
        {
          value = as.double(large_results[[n]][[run]]$output_data$data[[output]])
        }
      else 
        {
          print(sprintf("missing data in n:%d, run:%d, output:%d", n, run, output))
          Miss[[output]] = Miss[output]+1
        }
      
      if (!is.nan(value)) {FST[[n+1]][[output]] = FST[[n+1]][[output]] + (value)}
      else 
      {
        print(sprintf("the NaN value was found in n:%d, run:%d, output:%d", n, run, output))
        Nan_val[[output]] = Nan_val[[output]]+1 
      }
    }
    
    FOutputs = 12
  }
  
  print(Miss)
  print(Nan_val)
  
  for (output in 1:FOutputs) {
    FST[[n+1]][[output]] = FST[[n+1]][[output]]/(runs-Miss[[output]]-Nan_val[[output]])
  }
  
}

capture.output(FST, file = "FST.txt")
