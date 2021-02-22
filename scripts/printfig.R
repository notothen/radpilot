#### Script with functions to export high resolution figures
## 16/04/2018
## H. Christiansen

#### function for printing high resolution figures
printfig <-  function(name, plot) {
  filename <- (paste0(name,'.tiff'))
  tiff(file = filename, width = 7.25, height = 5, units = 'in', res = 500)
  par(mar = c(4.1,4.5,0.5,0.5), mgp = c(2.6,0.8,0), las = 1)
  plot(x)
  dev.off()
}

printfigPCA <-  function(name, plot) {
  filename <- (paste0(name,'.tiff'))
  tiff(file = filename, width = 5, height = 5, units = 'in', res = 500)
  par(mar = c(4.1,4.5,0.5,0.5), mgp = c(2.6,0.8,0), las = 1)
  plot(x)
  dev.off()
}

printfig3 <-  function(name, plot) {
  filename <- (paste0('/Figures/',name,'.tiff'))
  tiff(file = filename, width = 3.625, height = 6.5, units = 'in', res = 500)
  par(mar = c(4.1,4.5,0.5,0.5), mgp = c(2.6,0.8,0), las = 1)
  plot(x)
  dev.off()
}

printfigpdf <-  function(name, plot) {
  filename <- (paste0(name,'.pdf'))
  pdf(file = filename, width = 7.25, height = 5)
  par(mar = c(4.1,4.5,0.5,0.5), mgp = c(2.6,0.8,0), las = 1)
  plot(x)
  dev.off()
}

printfigpdf3 <-  function(name, plot) {
  filename <- (paste0('/Figures/',name,'.pdf'))
  pdf(file = filename, width = 3.625, height = 6.5)
  par(mar = c(4.1,4.5,0.5,0.5), mgp = c(2.6,0.8,0), las = 1)
  plot(x)
  dev.off()
}

printfigPCApdf <-  function(name, plot) {
  filename <- (paste0(name,'.pdf'))
  pdf(file = filename, width = 5, height = 5)
  par(mar = c(4.1,4.5,0.5,0.5), mgp = c(2.6,0.8,0), las = 1)
  plot(x)
  dev.off()
}
