R_path <- 'E:/Study-2024/09.RJ_2stage_resample_Oct24/4. Second_review CODE/results_c/'

library(ggplot2)
library(dplyr)
library(tidyr)

c_vec <- seq(0.1, 0.5, by = 0.1)
file_names <- paste0(R_path,"NewRes_PARA_c", c_vec, ".RData")

Res_List <- list()

for(i in 1:length(file_names)) {
  load(file_names[i]) # The name of Matrix in the .RData: Res_ALL
  Res_List[[i]] <- Res_ALL
}


Mse_Ave_c <- matrix(NA,nrow = 2,ncol =length(c_vec)) #Average MSEs and MSPEs
for (iN in 1:length(c_vec)){
  Mse_Ave_c[1,iN] <- mean(Res_List[[iN]][,1])
  Mse_Ave_c[2,iN] <- mean(Res_List[[iN]][,6])
}
results <- data.frame(
  par.c = c_vec,
  MSE =  Mse_Ave_c[1,],
  MSPE = Mse_Ave_c[2,]
)


results_long <- results %>%
  pivot_longer(cols = c(MSE, MSPE), 
               names_to = "Metric", 
               values_to = "Value")

ggplot(results_long, aes(x = par.c, y = Value, color = Metric, group = Metric)) +
  geom_line(lwd = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("MSE" = "#E69F00", "MSPE" = "#56B4E9")) +
  labs(title = "MSE and MSPE across Different c Values",
       x = "c",
       y = "Mean Value",
       color = "Metric") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(breaks = seq(0.1, 0.6, by = 0.1)) +
  scale_y_continuous(limits = c(0, max(results_long$Value) * 1.1))

