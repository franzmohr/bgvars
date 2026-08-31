
.make_new_plot_view <- function(y_names, x_names, title_text, lab_size) {
  
  n_y <- length(y_names)
  n_x <- length(x_names)
  
  # Define layout
  mat <- matrix(NA_integer_,
                n_y + 2 , # k_domestic + title row + column names
                n_x + 1) # n_ect + row names
  mat[1, ] <- 1
  mat[-1, 1] <- c(0, 2:(n_y + 1))
  mat[2, -1] <- (n_y + 1) + 1:n_x
  mat[-(1:2), -1] <- matrix(1:(n_y * n_x) + n_y + n_x + 1, n_y, n_x)
  graphics::layout(mat,
                   widths = c(lab_size, rep((1 - lab_size) / n_x, n_x)),
                   heights = c(.07, lab_size, rep((1 - lab_size) / n_y, n_y)))
  
  
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new() # Moves to the next subwindow
  graphics::text(0.5, 0.5, labels = title_text, cex = 1.5)
  # Fill rows
  graphics::par(mar = c(3, 0, 0, 0))
  for (j in y_names) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, labels = j, adj = 0.5)
  }
  # Fill columns
  graphics::par(mar = c(0, 0, 0, 0))
  for (j in x_names[1:n_x]) {
    graphics::plot.new(); graphics::text(0.5, 0.5, labels = j, adj = 0.5)
  }
  graphics::par(mar = c(3, 2.1, .5, 1))
}