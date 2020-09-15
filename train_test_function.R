  ### CREATE TEST AND TRAIN SETS
  ## Function
  set.seed(1234)
  create_train_test <- function(data, size = 0.8, train = TRUE) {
    n_row = nrow(data)
    total_row = size * n_row
    train_sample <- 1: total_row
    if (train == TRUE) {
      return (data[train_sample, ])
    } else {
      return (data[-train_sample, ])
    }
  }
