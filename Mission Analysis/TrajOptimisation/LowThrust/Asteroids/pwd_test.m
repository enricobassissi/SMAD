clear; clc;

n = 100;
t = linspace(0,10,n);
A = 5; f = 0.5;
y = A*sin(2*pi*f*t);
figure()
plot(t,y)

I = cumtrapz(y);
T = cumsum(t);

%% doesn't work
max_y = 0.2;
% MAX_Y = max_y*T;
for i = 1:n
    if (I(i)>=max_y)
        pwm(i) = 1;
        I
    else
        pwm(i) = 0;
    end
end
figure()
plot(t,y)
hold on
plot(t,max_y*pwm)