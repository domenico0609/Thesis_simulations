library(future)

plan(multisession, workers= 16)

runs = 100

big_humpy_0 <- slim_script_render(humpy_script_0, reps=runs)
big_results_0 <- slim_run(big_humpy_0, parallel = TRUE, progress = TRUE)

big_humpy_2 <- slim_script_render(humpy_script_2, reps=runs)
big_results_2 <- slim_run(big_humpy_2, parallel = TRUE, progress = TRUE)

big_humpy_10 <- slim_script_render(humpy_script_10, reps=runs)
big_results_10 <- slim_run(big_humpy_10, parallel = TRUE, progress = TRUE)

large_results <- list(nm0 = big_results_0, nm2 = big_results_2, nm10 = big_results_10)

