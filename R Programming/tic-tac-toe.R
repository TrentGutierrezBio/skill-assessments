install.packages("tidyverse")

library(tidyverse)



if (interactive()) {
  con <- stdin()
} else {
  con <- "stdin"
}
cat("X or O?")
symbol <- readLines(con = con, n = 1)


cat("X or O?")
empty_board <- (c(NA,NA,NA,NA,NA,NA,NA,NA,NA))

boardset <- matrix(data = empty_board, byrow = TRUE, nrow = 3)

boardset[1,1]

