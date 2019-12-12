width <- 64
height <- 48
d <- data.frame(row=rep(1:height, each=width),
                col=rep(1:width, height),
                x=rnorm(n=width * height))

# Construct and print the plot.
library(ggplot2)
plot <- (ggplot(data=d, aes(x=col, y=height + 1 - row, fill=x))
         + geom_raster()
         + scale_fill_gradient2()
)
pixel_size <- 10

BorderlessPlotPng <- function(plot, ...) {
  # Write a ggplot2 plot to an image file with no borders.
  #
  # Args:
  #   plot:  A ggplot2 plot object.
  #   ...:  Arguments passed to the png() function.
  require(grid)
  png(type='cairo', antialias=NULL, units='px', ...)
  print(plot
        + theme(plot.margin=unit(c(0, 0, -0.5, -0.5), 'line'),
                axis.text=element_blank(),
                axis.ticks=element_blank(),
                axis.title=element_blank(),
                legend.position='none')
        + scale_x_continuous(expand=c(0, 0))
        + scale_y_continuous(expand=c(0, 0))
        + labs(x=NULL, y=NULL)
  )
  dev.off()
}

BorderlessPlotPng(plot,
                  filename='test.png',
                  width=width * pixel_size,
                  height=height * pixel_size)
