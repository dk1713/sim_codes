sig = 3;
gaussian = @(x, w_0) exp(-( (x-0)/w_0 ).^( 2 ));

x = linspace(-10, 10, 200);

figure(10)
plot(x, gaussian(x, sig), 'b')