
u = linspace(-200, 200, 200);



for sigma = linspace(0, 1, 5)
   x = penalty_function(u, sigma);
   
   plot(x, u)
   hold on;
end

figure()
for sigma = linspace(0, 100, 5)
   x = penalty_function(u, sigma);

   plot(x, u)
   hold on;
end
