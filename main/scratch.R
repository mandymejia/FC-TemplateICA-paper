.libPaths(c(normalizePath("~/R/libs"), .libPaths()))
library(ciftiTools)

# Set Workbench binary path
ciftiTools.setOption('wb_path', '~/workbench/bin_rh_linux64/wb_command')

# Step 1: Create a basic xifti with data and medial wall mask
xii <- as.xifti(
  cortexL = rnorm(32492),
  cortexL_mwall = rep(FALSE, 32492),
  cortexR = rnorm(32492),
  cortexR_mwall = rep(FALSE, 32492)
)

# Step 2: Save PNG
png("test_xifti_plot.png", width = 1000, height = 800)
plot(xii,
     title = "Test Xifti Plot (PNG)",
     color_mode = "sequential",
     zlim = range(xii$data$cortex_left, na.rm = TRUE))
dev.off()

# Step 3: Save PDF
pdf("test_xifti_plot.pdf", width = 10, height = 8)
plot(xii,
     title = "Test Xifti Plot (PDF)",
     color_mode = "sequential",
     zlim = range(xii$data$cortex_left, na.rm = TRUE))
dev.off()



