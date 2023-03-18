function [p] = spline(N, D)
%n = number of steps in interpolation
%
%Output - x and y interpolation data, as well as plot of the spline graph.

%Import data


fTemp = fopen('bessel.txt')
data = textscan(fTemp, '%f%f%f%f', 'HeaderLines', 1)
%
rho= data{1}
j0 = data{2}
j1 = data{3}
j2 = data{4}

%Create tridiagonal
L = length(rho)

%Create h vector to store h values
h = zeros(L,1)
for i=1:L-1
    h(i) = rho(i+1)-rho(i)
end
h(length(h)) = 1

A = zeros(L,L)
A(1,1) = 2*h(1)
A(1,2) = h(1)
A(L,L-1) = h(L-1)
A(L,L) = 2*(h(L-1))

for i = 2:L-1
   
    A(i,i) = 2*(h(i+1)+h(i))
    A(i,i-1) = h(i)
    A(i,i+1) = h(i)
end


%Create B matrix
B = zeros(L,1)


if D ~= 0
    B(1,1) = 6*((j1(2)-j1(1)))/h(1) - 6*D
else
    B(1,1) = 6*((j1(2)-j1(1)))/h(1) - 6*(1/2)*(j0(1)-j2(1)) 
end
B(L,1) = -6*((j1(L) - j1(L-1)))/(h(L)) + 6*(1/2)*(j0(L)-j2(L))
for i = 2:L-1
   B(i,1) = (6*(j1(i+1) -j1(i))/(h(i))) - 6*(j1(i)-j1(i-1))/h(i-1)
end

%Now solve for x vector using tridiagonal script from #1

secDir = tridiagonal(A,B);

%These give us the second derivatives of p(x).  Now plug into equations:

plot(rho,j1, '.', markersize = 20)
xlabel('x')
ylabel('J1(x)')
title('Cubic spline with D = 3')

hold on
line(xlim(), [0,0], 'LineWidth', 1, 'Color', 'k')

p = zeros(1,N)


for k = 1:L-1
    %Iterate over each pair of points
    
    px = linspace(rho(k), rho(k+1),N)
    hk = rho(k+1) - rho(k)
    p = zeros(1,N)
    for i = 1:N
         a = (secDir(k+1)-secDir(k))/(6*hk)
         b = secDir(k)/2
         c = (j1(k+1)-j1(k))/hk - (hk*secDir(k+1))/6 - (hk*secDir(k))/3
         d = j1(k)
         p(i) = d + c*(i/N)+b*(i/N).^2 + a*(i/N).^3
         
    end
    if D == 0
        plot(px, p, 'b.')
    else
        plot(px, p, 'k.')
    end
end





%plot(px1,p1, '.')