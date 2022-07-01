function loss = gamma_fmin(gamma, accu_inHosp, newly_recov)

y = gamma * accu_inHosp;
logy = max(0, log(y));

loss = -sum(newly_recov(:) .* logy(:) - y(:));