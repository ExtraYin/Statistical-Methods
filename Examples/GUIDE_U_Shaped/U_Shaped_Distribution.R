s1_x = rep(-10:10, each = 10)
s1_y = s1_x^2
s1_x = s1_x + rnorm(210)
s1_y = s1_y + 10 * rnorm(210)

s2_x = rep(0:20, each = 10)
s2_y = -(s2_x-10)^2+130
s2_x = s2_x + rnorm(210)
s2_y = s2_y + 10 * rnorm(210)

# plot(s1_x, s1_y, col = "red", xlim = c(-15,25), ylim = c(-10, 200))
# points(s2_x, s2_y, col = "blue")

d = data.frame(x = c(s1_x, s2_x), y = c(s1_y, s2_y), label = c(rep('red', 210), rep('blue', 210)))
plot(d$x, d$y, col = d$label)

write.csv(d, file = "U_Shaped.txt", row.names = FALSE)
