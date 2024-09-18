pct_dx_young <- (1-.561) # Undiagnosed 13-24
pct_dx_older <- (1-.295) # Undiagnosed 25-34

ci_dx_young <- 1 - c(.609, .619)

ci_dx_older <- 1 - c(.282, .300)

length_young <- length(16:24)
length_older <- length(25:29)
total_length <- length_young + length_older

age_weighted <- pct_dx_young*(length_young/total_length) + pct_dx_older*(length_older/total_length)

0.55333

total <- (1 -.174)
black <- (1 -.206)
white <- (1 -.123)
hispanic <- (1 - .208)
(black/total) * age_weighted

(c(black, hispanic, white)/total)*age_weighted

total <- (1 -.174)
black <- (1 -.206)
white <- (1 -.123)
hispanic <- (1 - .208)
(black/total) * age_weighted

(c(black, hispanic, white)/total)*age_weighted


.804/.814 * age_weighted

0.5465356428

.804*(age_weighted/.814)

.6797707*.804

(1-c(.235, .175)) * .6797707

(1-c(.240, .173)) * .6797707

(1-c(.150, .095)) * .6797707

(1-c(mean(c(0, .096, 0, .089)),
     mean(c(.447, .34, .503, .176)))) * .6797707


