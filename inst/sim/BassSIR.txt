deriv(S) <- - foi * S
deriv(I) <- foi * S - dropout * I

initial(S) <- m - I0 - A0
initial(I) <- I0

output(Inc) <- foi * S
output(FOI) <- foi
output(PrEx) <- foi_ex / foi * 100
output(Re) <- beta / dropout * S / m
output(R0) <- beta / dropout

foi <- foi_en + foi_ex
foi_en <- beta * (I / m)
foi_ex <- kappa

m <- user(1000)
I0 <- user(5)
A0 <- user(5)

kappa <- user(0.01)
beta <- user(0.3)
dropout <- user(0)
