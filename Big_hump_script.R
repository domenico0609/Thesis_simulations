library(future)

plan(multisession, workers= 16)

runs = 100
NMVar = 3
FOutputs = 12

#big_humpy_0 <- slim_script_render(humpy_script_0, reps=runs)
#big_results_0 <- slim_run(big_humpy_0, parallel = TRUE, progress = TRUE)

big_humpy_2 <- slim_script_render(humpy_script_2, reps=runs)
big_results_2 <- slim_run(big_humpy_2, parallel = TRUE, progress = TRUE)

big_humpy_10 <- slim_script_render(humpy_script_10, reps=runs)
big_results_10 <- slim_run(big_humpy_10, parallel = TRUE, progress = TRUE)

large_results <- list(nm0 = big_results_0, nm2 = big_results_2, nm10 = big_results_10)

FST=data.frame(
  FST = c("mito Fst p1&p2", "mito Fst p1&p3", "mito Fst p1&p4", "mito Fst p2&p3", "mito Fst p2&p4", "mito Fst p3&p4", "nuc Fst p1&p2", "nuc Fst p1&p3", "nuc Fst p1&p4", "nuc Fst p2&p3", "nuc Fst p2&p4", "nuc Fst p3&p4"),
  NM0 = double(FOutputs),
  NM2 = double(FOutputs),
  NM10 = double(FOutputs))

for (n in 1:NMVar) {
  if (n == 1) runs = 2 else runs = 100
  print(runs)
  for (run in 1:runs) {
    for (output in 1:FOutputs) {
      value = as.double(large_results[[n]][[run]]$output_data$data[[output]])
      if (!is.nan(value)) {FST[[n+1]][[output]] = FST[[n+1]][[output]] + (value/runs)}
      else {print(sprintf("the NaN value was found in n:%d, run:%d, output:%d", n, run, output)) }
    }  
  }
}

capture.output(FST, file = "FST.txt")
