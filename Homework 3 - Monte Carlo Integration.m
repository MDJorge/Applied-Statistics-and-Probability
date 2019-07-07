%% HW3, name: Jorge Monzon Diaz, email:jmonzon000@citymail.cuny.edu
%% Problem 1
% Finding Confidence Intervals using Monte Carlo Integration

%part a
n=100000;
x=rand(n,100); %generate n 100-dimensional rvs
out=hw4h1(x); %put rvs through function
mn=mean(out);
st=std(out);
error=1.96.*st./sqrt(n-1); %find pivots
CI1a=[mn-error,mn+error] %CI for 95%= mean +- 1.96*error

%part b
n=100000;
x2=rand(n,100)*1.05; %generate n 100-dimensional rvs in [0,1.05]
out2=hw4h1(x2)./(1./(1.05.^100)); %put rvs through function and use  importance sampling
mn2=mean(out2);
st2=std(out2);
error2=1.96.*st2./sqrt(n-1); %find pivots
CI1b=[mn2-error,mn2+error] %CI for 95%= mean +- 1.96*error

%part c
inside=@(x)sum((x.^2),2)<=36; %create condition for radius of ball
%sum x_i^2 from i=1 to n <=36

n=100000;
x3=rand(n,100); %generate n 100-dimensional rvs
samples=inside(x3).*x3; %puts 100-d rvs through the condition checking function
out3=hw4h1(x3); %puts samples fitting condition through function
mn3=mean(out3);
st3=std(out3);
error3=1.96.*st3./sqrt(n-1); %find pivots
CI1c=[mn3-error,mn3+error] %CI for 95%= mean +- 1.96*error

%% functions
function h=hw4h1(x)
% x should be a Nx100 matrix

if size(x,2) ~= 100, error('wrong size'), end

h=abs( sin( 2*pi*x(:,1).*sum(x,2))).*((cos( 2*pi*x(:,2).* sum(x.^2,2))).^2);
end

function h=quest2fn(x,y,z)
    if (size(x,1) ~= 1) && (size(y,1) ~= 1), error('wrong size'), end
    if sqrt(x.^2+y.^2) > 1
        h=0;
    elseif sqrt(x.^2+y.^2) > 1/2 && sqrt(x.^2+y.^2) <= 1
        h=z.*exp(x.*y).*exp(-z.*(x.^2+y.^2));
    else
        h=z.*cos(3.*x.*y).^2.*exp(-z.*(x.^2+y.^2));
    end
end
