library(ggplot2)

ggplot(data, aes(x=quarter, y=cum_sales)) +
  geom_line() +
  labs(x='Quarter', y='Cumulative Sales')



OSJ <- read.table("clipboard",header = T)
ARATH <- read.table("clipboard",header = T)

OSJ$C


plot(OSJ$Degree,type = "h")
plot(cumsum(OSJ$Degree)/sum(OSJ$Degree)*100,col="red",type="s")
points(OSJ$Degree,type = "h")

# plot(cumsum(OSJ$Degree), type="s")
# points(OSJ$Degree,cumsum(OSJ$Degree),col="red")

library(ggplot2);library(ggpubr)

cumulative <- data.frame(Degree=OSJ$Degree)
cumulative$Degree <- sort(cumulative$Degree)
cumulative$cumulative <- cumsum(OSJ$Degree/sum(OSJ$Degree))


# Fix manually the maximum value of y-axis
ymax_OSJ <- max(OSJ$Degree,na.rm = T)
OSJ_plot <- ggplot(data=cumulative,aes(x=Degree)) + 
   geom_histogram(binwidth = 0.3)+
  scale_x_log10(name = 'Degree',breaks=seq(0,5500,500))+
  geom_line(aes(x=Degree,y=cumulative*ymax_OSJ), col="red", lwd=1)+
  scale_y_continuous(name = 'Counts', breaks=seq(0,5500,500),sec.axis = sec_axis(~./ymax_OSJ*100,
                                                                   name = "Cumulative percentage [%]"))+
  ggtitle("OSJ genes undirected graph (Cumulative sum of degrees)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


cumulative2 <- data.frame(Degree=ARATH$Degree)
cumulative2$Degree <- sort(cumulative2$Degree)
cumulative2$cumulative <- cumsum(ARATH$Degree/sum(ARATH$Degree))



ymax_ARA <- max(ARATH$Degree,na.rm = T)
ARATH_plot <- ggplot(data=cumulative2,aes(x=Degree)) + 
  geom_histogram(binwidth = 0.3)+
  scale_x_log10(name = 'Degree',breaks=seq(0,9500,500))+
  geom_line(aes(x=Degree,y=cumulative*ymax_ARA), col="red", lwd=1)+
  scale_y_continuous(name = 'Counts', breaks=seq(0,9500,500),sec.axis = sec_axis(~./ymax_ARA*100,
                                                          name = "Cumulative percentage [%]"))+
  ggtitle("ARATH genes undirected graph (Cumulative sum of degrees)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



figure <- ggarrange(OSJ_plot, ARATH_plot,
                    labels = c("A", "B"),
                    ncol = 1, nrow = 2)
figure
